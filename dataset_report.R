# ====================================================================
# dataset_report.R  — unified per-dataset diagnostic report
# ====================================================================
# Replaces:  dataset_overview_v2.R, dataset-slices_overview.R
#
# Three entry points:
#
#   .prepare_dataset()  — internal workhorse (steps 1-4)
#     1. Loads deployments & observations
#     2. Finds 12-month anchor slices (or uses pre-computed ones)
#     3. Joins dataset_meta.xlsx (provider, country, camera setup)
#     4. Returns a list with single-row metrics tibble + intermediate objects
#
#   dataset_metrics()   — public: runs steps 1-4 and returns the metrics row
#
#   dataset_report()    — public: runs steps 1-4 then produces a PDF:
#        Page 1    — Dashboard: metrics + station map + effort sparkline
#        Page 2    — Effort timeline + anchor table + protocol effort table
#        Page 3    — Species × month heatmap with protocol window markers
#        Pages 4+  — Per-slice: station richness map + accumulation curve
#        Pages N+  — Species panels (top 5 + bottom 5 by detections)
#        Last      — Rejected species table
# ====================================================================

# ── packages ───────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(sf)
  library(cowplot)
  library(ggspatial)
  library(scales)
  library(grid)
  library(geosphere)
  library(gridExtra)
  library(readxl)
  library(furrr)
})

.helpers_path <- normalizePath(
  file.path(dirname(rstudioapi::getSourceEditorContext()$path), "helpers.R")
)
source(.helpers_path)

# ── private cache ──────────────────────────────────────────────────
.report_cache <- new.env(parent = emptyenv())

# daily_deployment_effort_fast() is defined in helpers.R — no local copy needed.

# ── helper: load world sf ──────────────────────────────────────────
get_world_sf <- function() {
  if (is.null(.report_cache$world_sf)) {
    .report_cache$world_sf <- if (requireNamespace("rnaturalearth", quietly = TRUE))
      rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
    else if (requireNamespace("maps", quietly = TRUE))
      st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))
    else NULL
  }
  .report_cache$world_sf
}

# protocol_calendar_dates() is defined in helpers.R — returns all protocol
# windows (SNAP_EU_CORE, SNAP_EU_BUFFER, EOW_EARLY, EOW_LATE) for a given anchor period.



# ====================================================================
# DATASET META JOIN
# ====================================================================

read_dataset_meta <- function(
    meta_path = "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/dataset_meta.xlsx"
) {
  if (!file.exists(meta_path)) {
    message("dataset_meta.xlsx not found at: ", meta_path)
    return(NULL)
  }
  dm <- readxl::read_xlsx(meta_path) |>
    filter(!is.na(dataset)) |>
    select(dataset, any_of(c("provider", "country", "target", "bait",
                              "trail", "spacing", "height", "micro_loc",
                              "note")))

  # Normalise dataset names for fuzzy matching
  dm |> mutate(dataset_join = tolower(gsub("[^a-z0-9]", "", dataset)))
}


# ====================================================================
# INTERNAL: data loading + metrics computation (steps 1-5)
# ====================================================================
.prepare_dataset <- function(
    dataset_name,
    keep_species    = spp_keep,
    anchors_tbl     = NULL,
    dataset_meta    = NULL,
    root            = "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets"
) {

  stopifnot(is.character(keep_species), length(keep_species) > 0)

  # ── load data + detect columns ──────────────────────────────────
  archive_path <- resolve_dataset_path(dataset_name, root, include_zips = TRUE)
  deploy <- read_csv_src(archive_path, "deployments\\.csv")
  obs    <- read_csv_src(archive_path, "observations\\.csv")

  cols <- set_camtrap_cols(obs, deploy)
  list2env(cols, envir = environment())

  # ── daily effort ────────────────────────────────────────────────
  daily_dep_full <- daily_deployment_effort_fast(
    deploy, start_col, end_col, drop_leap_day = FALSE
  )
  if (!nrow(daily_dep_full)) stop("No valid deployment dates for: ", dataset_name)

  # ── parse observations ──────────────────────────────────────────
  obs_ts <- obs |>
    transmute(
      deploymentID = .data[[depl_id_obs]],
      species      = .data[[sci_col]],
      timestamp    = parse_ts_safe(.data[[date_col]]),
      .empty       = {
        if ("is_empty" %in% names(obs)) {
          v <- obs$is_empty
          if (is.list(v)) map_lgl(v, \(x) identical(x, TRUE)) else as.logical(v)
        } else rep(FALSE, n())
      }
    ) |>
    filter(!is.na(timestamp))

  # ── observation span ────────────────────────────────────────────
  obs_kept <- obs_ts |>
    filter(!is.na(species), nzchar(species), species %in% keep_species)

  if (nrow(obs_kept) > 0) {
    span_tbl <- obs_kept
  } else if (nrow(obs_ts) > 0) {
    span_tbl <- obs_ts
    message("No observations for kept species; using all-obs span: ", dataset_name)
  } else {
    span_tbl <- NULL
    message("No valid timestamps; using deployment span: ", dataset_name)
  }

  if (!is.null(span_tbl)) {
    start_date <- as_date(min(span_tbl$timestamp, na.rm = TRUE))
    end_date   <- as_date(max(span_tbl$timestamp, na.rm = TRUE))
  } else {
    start_date <- min(daily_dep_full$date, na.rm = TRUE)
    end_date   <- max(daily_dep_full$date, na.rm = TRUE)
  }

  duration_days <- as.integer(end_date - start_date + 1L)
  duration_years <- round(duration_days / 365.25, 1)

  # ── restrict effort to observation span ─────────────────────────
  daily_dep <- tibble(date = seq.Date(start_date, end_date, by = "day")) |>
    left_join(daily_dep_full, by = "date") |>
    mutate(n_active_deployments = replace_na(n_active_deployments, 0L)) |>
    arrange(date)

  trap_days_total <- sum(daily_dep$n_active_deployments)
  mean_active     <- mean(daily_dep$n_active_deployments)

  r <- rle(daily_dep$n_active_deployments == 0L)
  longest_gap <- if (any(r$values)) max(r$lengths[r$values]) else 0L

  # ── temporal coverage & data quality ──────────────────────────
  days_with_effort <- sum(daily_dep$n_active_deployments > 0L)
  temporal_coverage_pct <- round(100 * days_with_effort / max(duration_days, 1L), 1)

  # Data quality: missing timestamps, missing coordinates
  n_obs_no_ts <- sum(is.na(parse_ts_safe(obs[[date_col]])))
  n_dep_no_coords <- deploy |>
    transmute(
      lat = parse_num_safe(.data[[lat_col]]),
      lon = parse_num_safe(.data[[lon_col]])
    ) |>
    filter(!is.finite(lat) | !is.finite(lon)) |>
    nrow()
  n_dep_bad_dates <- deploy |>
    transmute(
      start = parse_ts_safe(.data[[start_col]]),
      end   = parse_ts_safe(.data[[end_col]])
    ) |>
    filter(is.na(start) | is.na(end) | end <= start) |>
    nrow()

  # ── stations & geometry ─────────────────────────────────────────
  loc_coords <- deploy |>
    transmute(
      location_id = .data[[loc_col]],
      lat = parse_num_safe(.data[[lat_col]]),
      lon = parse_num_safe(.data[[lon_col]])
    ) |>
    filter(!is.na(location_id), is.finite(lat), is.finite(lon)) |>
    summarise(lat = mean(lat), lon = mean(lon), .by = location_id)

  n_locations <- n_distinct(deploy[[loc_col]])

  # Deduplicate to unique physical locations — datasets with rotating
  # deployments (e.g. BE-Leuven) assign many locationIDs to the same
  # coordinates; spatial metrics should reflect physical locations only.
  unique_coords <- loc_coords |> distinct(lat, lon)
  n_locations_unique <- nrow(unique_coords)

  # Map each locationID to its unique-coordinate row index
  loc_coords$uc_idx <- match(
    paste(loc_coords$lat, loc_coords$lon),
    paste(unique_coords$lat, unique_coords$lon)
  )

  if (n_locations_unique >= 2) {
    dmat <- geosphere::distm(unique_coords[, c("lon", "lat")],
                             fun = geosphere::distHaversine) / 1000
    trap_array_km <- max(dmat)
    diag(dmat) <- NA_real_
    median_nn_km <- median(apply(dmat, 1, min, na.rm = TRUE))
    cen <- geosphere::centroid(unique_coords[, c("lon", "lat")])
    centroid_lon <- cen[1]; centroid_lat <- cen[2]

    # Mean nearest-neighbour distance among simultaneously active cameras.
    # Uses NN (not mean pairwise) so clustered designs get their within-
    # cluster spacing, not the between-cluster spread.
    dep_loc <- deploy |>
      transmute(
        location_id = .data[[loc_col]],
        dep_start   = as_date(parse_ts_safe(.data[[start_col]])),
        dep_end     = as_date(parse_ts_safe(.data[[end_col]]))
      ) |>
      filter(!is.na(dep_start), !is.na(dep_end), dep_end > dep_start,
             location_id %in% loc_coords$location_id)

    all_days <- seq.Date(start_date, end_date, by = "day")
    n_days   <- length(all_days)
    n_uc     <- n_locations_unique

    # Activity matrix: unique coordinates × days
    # A physical location is active if ANY of its locationIDs are active.
    active_mat <- matrix(FALSE, nrow = n_uc, ncol = n_days)
    for (.r in seq_len(nrow(dep_loc))) {
      loc_idx <- loc_coords$uc_idx[match(dep_loc$location_id[.r],
                                         loc_coords$location_id)]
      if (is.na(loc_idx)) next
      d1 <- max(1L, as.integer(dep_loc$dep_start[.r] - start_date) + 1L)
      d2 <- min(n_days, as.integer(dep_loc$dep_end[.r] - start_date) + 1L)
      if (d1 <= d2) active_mat[loc_idx, d1:d2] <- TRUE
    }

    daily_mean_nn <- vapply(seq_len(n_days), function(d) {
      idx <- which(active_mat[, d])
      if (length(idx) < 2L) return(NA_real_)
      sub <- dmat[idx, idx, drop = FALSE]
      mean(apply(sub, 1, min, na.rm = TRUE))
    }, numeric(1))

    mean_cam_dist_km <- mean(daily_mean_nn, na.rm = TRUE)
    if (!is.finite(mean_cam_dist_km)) mean_cam_dist_km <- NA_real_

  } else {
    trap_array_km <- median_nn_km <- mean_cam_dist_km <- NA_real_
    centroid_lat <- if (n_locations_unique == 1) unique_coords$lat[1] else NA_real_
    centroid_lon <- if (n_locations_unique == 1) unique_coords$lon[1] else NA_real_
  }

  # ── observation counts ──────────────────────────────────────────
  obs_span <- obs_ts |>
    filter(
      timestamp >= as.POSIXct(start_date, tz = "UTC"),
      timestamp <= as.POSIXct(end_date + days(1), tz = "UTC")
    )
  n_obs_all  <- nrow(obs_span)
  n_spp_all  <- obs_span |>
    filter(!is.na(species), nzchar(species)) |>
    pull(species) |> n_distinct()

  obs_filtered <- obs_span |>
    filter(!is.na(species), nzchar(species), species %in% keep_species)
  n_obs_filt <- nrow(obs_filtered)
  n_spp_filt <- n_distinct(obs_filtered$species)

  # ── anchors ─────────────────────────────────────────────────────
  ds_anchors <- if (!is.null(anchors_tbl)) {
    anchors_tbl |> filter(grepl(dataset_name, dataset, fixed = TRUE))
  } else {
    a <- get0("anchors", envir = .GlobalEnv)
    if (!is.null(a)) {
      a |> filter(grepl(dataset_name, dataset, fixed = TRUE))
    } else {
      # Find anchors on the fly
      if (exists("find_anchors", mode = "function")) {
        find_anchors(archive_path, verbose = FALSE)
      } else {
        tibble()
      }
    }
  }

  # ── dataset metadata join ───────────────────────────────────────
  if (is.null(dataset_meta)) dataset_meta <- read_dataset_meta()
  meta_cols <- tibble(
    provider = NA_character_, country = NA_character_,
    target = NA_character_, bait = NA_character_,
    trail = NA_character_, spacing = NA_character_,
    height = NA_character_
  )
  if (!is.null(dataset_meta)) {
    join_key <- tolower(gsub("[^a-z0-9]", "", dataset_name))
    matched <- dataset_meta |> filter(dataset_join == join_key)
    if (nrow(matched) == 1) {
      meta_cols <- matched |>
        select(any_of(c("provider", "country", "target", "bait",
                         "trail", "spacing", "height")))
    }
  }

  # ── build metrics row (no climate variables — those live in env_covs
  #    and are shown in the PDF dashboard but not exported to the table) ──
  metrics_row <- bind_cols(
    tibble(
      dataset_name              = dataset_name,
      start_date                = start_date,
      end_date                  = end_date,
      duration_days             = duration_days,
      stations                  = n_locations_unique,
      mean_active_deployments   = mean_active,
      temporal_coverage_pct     = temporal_coverage_pct,
      longest_zero_effort_gap_days = as.integer(longest_gap),
      trap_days_total           = trap_days_total,
      array_diameter_km         = trap_array_km,
      median_nn_dist_km         = median_nn_km,
      mean_active_cam_dist_km   = mean_cam_dist_km,
      centroid_lat              = centroid_lat,
      centroid_lon              = centroid_lon,
      n_obs_all                 = n_obs_all,
      n_spp_all                 = n_spp_all,
      n_obs_filtered            = n_obs_filt,
      n_spp_filtered            = n_spp_filt,
      n_slices                  = nrow(ds_anchors),
      n_obs_no_timestamp        = n_obs_no_ts,
      n_dep_no_coords           = n_dep_no_coords,
      n_dep_bad_dates           = n_dep_bad_dates
    ),
    meta_cols
  )

  list(
    metrics_row     = metrics_row,
    dataset_name    = dataset_name,
    archive_path    = archive_path,
    deploy          = deploy,
    obs_ts          = obs_ts,
    obs_span        = obs_span,
    daily_dep       = daily_dep,
    daily_dep_full  = daily_dep_full,
    loc_coords      = loc_coords,
    ds_anchors      = ds_anchors,
    meta_cols       = meta_cols,
    keep_species    = keep_species,
    start_date      = start_date,
    end_date        = end_date,
    duration_days   = duration_days,
    duration_years  = duration_years,
    n_locations        = n_locations,
    n_locations_unique = n_locations_unique,
    mean_active        = mean_active,
    temporal_coverage_pct = temporal_coverage_pct,
    longest_gap     = longest_gap,
    trap_days_total = trap_days_total,
    trap_array_km   = trap_array_km,
    median_nn_km    = median_nn_km,
    mean_cam_dist_km = mean_cam_dist_km,
    centroid_lat    = centroid_lat,
    centroid_lon    = centroid_lon,
    n_obs_all       = n_obs_all,
    n_spp_all       = n_spp_all,
    n_obs_filt      = n_obs_filt,
    n_spp_filt      = n_spp_filt,
    n_obs_no_ts     = n_obs_no_ts,
    n_dep_no_coords = n_dep_no_coords,
    n_dep_bad_dates = n_dep_bad_dates,
    loc_col         = loc_col,
    depl_id_deploy  = depl_id_deploy,
    start_col       = start_col,
    end_col         = end_col,
    lat_col         = lat_col,
    lon_col         = lon_col
  )
}


# ====================================================================
# PUBLIC: metrics only (no PDF)
# ====================================================================
dataset_metrics <- function(
    dataset_name,
    keep_species    = spp_keep,
    anchors_tbl     = NULL,
    dataset_meta    = NULL,
    root            = "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets"
) {
  .prepare_dataset(dataset_name, keep_species, anchors_tbl,
                   dataset_meta, root)$metrics_row
}


# ====================================================================
# PUBLIC: metrics + PDF report
# ====================================================================
dataset_report <- function(
    dataset_name,
    keep_species    = spp_keep,
    anchors_tbl     = NULL,
    dataset_meta    = NULL,
    root            = "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets",
    output_pdf      = NULL,
    pdf_only        = FALSE
) {
  pal <- list(
    accent       = "#2166ac",
    accent2      = "#b2182b",
    anchor_fill  = alpha("#2166ac", 0.12),
    anchor_border = alpha("#2166ac", 0.40),
    core_fill    = alpha("#d6604d", 0.15),
    core_border  = "#d6604d",
    buffer_fill  = alpha("#4393c3", 0.10),
    buffer_border = "#4393c3",
    eow_e_fill   = alpha("#5aae61", 0.12),
    eow_e_border = "#5aae61",
    eow_l_fill   = alpha("#984ea3", 0.12),
    eow_l_border = "#984ea3",
    station      = "#d73027",
    effort       = "grey30",
    bg           = "grey98",
    grid         = "grey90"
  )

  prep <- .prepare_dataset(dataset_name, keep_species, anchors_tbl,
                           dataset_meta, root)
  list2env(prep, envir = environment())

  # ================================================================
  # PAGE 1 — DASHBOARD
  # ================================================================

  # Build data quality flag string
  dq_issues <- c(
    if (n_obs_no_ts > 0)    paste0(comma(n_obs_no_ts), " obs no timestamp") else NULL,
    if (n_dep_no_coords > 0) paste0(n_dep_no_coords, " dep no coords") else NULL,
    if (n_dep_bad_dates > 0) paste0(n_dep_bad_dates, " dep bad dates") else NULL
  )
  dq_label <- if (length(dq_issues)) paste(dq_issues, collapse = "; ") else "OK"

  metrics_df <- tibble(
    Metric = c(
      "Period", "Duration", "Stations",
      "Mean active cameras", "Temporal coverage",
      "Longest gap", "Trap-days",
      "Array diameter", "Median NN dist.", "Mean active NN dist.", "Centroid",
      "Species (kept / all)", "Observations (kept / all)",
      "Anchor slices", "Data quality",
      if (!is.na(meta_cols$country[1])) "Country" else NULL
    ),
    Value = c(
      paste(format(start_date, "%Y-%m-%d"), "\u2013", format(end_date, "%Y-%m-%d")),
      paste0(comma(duration_days), " days (", duration_years, " yr)"),
      if (n_locations != n_locations_unique)
        paste0(comma(n_locations_unique), " locations (", comma(n_locations), " locationIDs)")
      else comma(n_locations_unique),
      formatC(mean_active, digits = 1, format = "f"),
      paste0(temporal_coverage_pct, "% of days with \u22651 camera"),
      paste0(comma(as.integer(longest_gap)), " days"),
      comma(trap_days_total),
      paste0(formatC(trap_array_km, digits = 1, format = "f"), " km"),
      paste0(formatC(median_nn_km, digits = 2, format = "f"), " km"),
      paste0(formatC(mean_cam_dist_km, digits = 2, format = "f"), " km"),
      paste0(round(centroid_lat, 2), "\u00b0N, ", round(centroid_lon, 2), "\u00b0E"),
      paste0(n_spp_filt, " / ", n_spp_all),
      paste0(comma(n_obs_filt), " / ", comma(n_obs_all)),
      as.character(nrow(ds_anchors)),
      dq_label,
      if (!is.na(meta_cols$country[1])) meta_cols$country[1] else NULL
    )
  )

  tbl_theme <- gridExtra::ttheme_minimal(
    core = list(
      fg_params = list(hjust = 0, x = 0.02, fontsize = 10),
      bg_params = list(fill = c(pal$bg, "white"), col = NA)
    ),
    colhead = list(
      fg_params = list(hjust = 0, x = 0.02, fontsize = 10, fontface = "bold"),
      bg_params = list(fill = alpha(pal$accent, 0.15), col = NA)
    )
  )
  tbl_grob <- tableGrob(metrics_df, rows = NULL, theme = tbl_theme)
  tbl_grob$widths <- unit(c(0.42, 0.58), "npc")

  # Station map
  loc_sf <- st_as_sf(loc_coords, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  bbox <- st_bbox(loc_sf)
  bb_pad <- pad_bbox_safe(bbox, 0.15)
  world_sf <- get_world_sf()

  if (!is.null(world_sf)) {
    cen_sfc <- st_centroid(st_as_sfc(bbox))
    country_sf <- suppressWarnings(
      world_sf[st_intersects(world_sf, cen_sfc, sparse = FALSE), ]
    )
    if (nrow(country_sf) == 0) country_sf <- world_sf

    map_plot <- ggplot() +
      geom_sf(data = country_sf, fill = "grey96", color = "grey75", linewidth = 0.25) +
      geom_sf(data = loc_sf, color = pal$station, size = 1.4, alpha = 0.85, shape = 16) +
      coord_sf(xlim = c(bb_pad["xmin"], bb_pad["xmax"]),
               ylim = c(bb_pad["ymin"], bb_pad["ymax"]), expand = FALSE) +
      annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.6,
                       line_width = 0.4, height = unit(0.12, "cm")) +
      theme_void() +
      theme(plot.margin = margin(2, 2, 2, 2))

    minimap <- ggplot() +
      geom_sf(data = country_sf, fill = "white", color = "grey40", linewidth = 0.15) +
      geom_sf(data = st_as_sfc(bb_pad), fill = NA, color = pal$accent2, linewidth = 0.6) +
      theme_void() +
      theme(plot.background = element_rect(fill = "white", color = "grey70", linewidth = 0.3))
  } else {
    map_plot <- ggplot() +
      geom_point(data = loc_coords, aes(lon, lat), color = pal$station,
                 size = 1.4, alpha = 0.85) +
      theme_void()
    minimap <- ggplot() + theme_void()
  }

  map_with_inset <- ggdraw() +
    draw_plot(map_plot, 0, 0, 1, 1) +
    draw_plot(minimap, 0.68, 0.65, 0.28, 0.30)

  # Effort sparkline
  effort_spark <- ggplot(daily_dep, aes(date, n_active_deployments)) +
    geom_area(fill = alpha(pal$accent, 0.15), color = pal$accent, linewidth = 0.3) +
    scale_x_date(date_labels = "%b %Y", expand = expansion(mult = 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Active cameras") +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
      plot.margin = margin(4, 8, 4, 8)
    )

  page1 <- ggdraw() +
    draw_label(dataset_name, x = 0.04, y = 0.98, hjust = 0, vjust = 1,
               fontface = "bold", size = 16, colour = pal$accent) +
    draw_label("Dataset Overview", x = 0.04, y = 0.955, hjust = 0, vjust = 1,
               size = 10, colour = "grey50") +
    draw_grob(tbl_grob, x = 0.02, y = 0.48, width = 0.52, height = 0.45,
              hjust = 0, vjust = 0) +
    draw_plot(map_with_inset, x = 0.52, y = 0.42, width = 0.46, height = 0.55) +
    draw_line(x = c(0.04, 0.96), y = c(0.44, 0.44), color = "grey80", linewidth = 0.3) +
    draw_label("Daily effort", x = 0.04, y = 0.42, hjust = 0, vjust = 1,
               fontface = "bold", size = 10, colour = "grey40") +
    draw_plot(effort_spark, x = 0.02, y = 0.02, width = 0.96, height = 0.38)


  # ================================================================
  # PAGE 2 — EFFORT TIMELINE + ANCHORS + PROTOCOL EFFORT
  # ================================================================

  n_slices_label <- if (nrow(ds_anchors) > 0) {
    paste0(nrow(ds_anchors), " anchor slice",
           if (nrow(ds_anchors) > 1) "s" else "",
           " | CORE (red) | BUFFER (blue) | EOW-E (green) | EOW-L (purple)")
  } else {
    "No anchor slices found"
  }

  effort_plot <- ggplot(daily_dep, aes(date, n_active_deployments)) +
    geom_area(fill = alpha(pal$effort, 0.08), color = pal$effort, linewidth = 0.35)

  # Add protocol window shading (behind slice borders)
  if (nrow(ds_anchors) > 0) {
    snap_wins <- purrr::map_dfr(seq_len(nrow(ds_anchors)), \(i) {
      protocol_calendar_dates(ds_anchors$anchor_start[i], ds_anchors$anchor_end[i])
    }) |> distinct()

    win_fill <- c(SNAP_EU_CORE = pal$core_fill, SNAP_EU_BUFFER = pal$buffer_fill,
                  EOW_EARLY = pal$eow_e_fill, EOW_LATE = pal$eow_l_fill)
    for (j in seq_len(nrow(snap_wins))) {
      w <- snap_wins[j, ]
      effort_plot <- effort_plot +
        annotate("rect", xmin = w$win_start, xmax = w$win_end,
                 ymin = -Inf, ymax = Inf,
                 fill = win_fill[w$window] %||% alpha("grey50", 0.1),
                 color = NA)
    }
  }

  # Add anchor slice borders and labels
  if (nrow(ds_anchors) > 0) {
    y_label <- max(daily_dep$n_active_deployments, na.rm = TRUE) * 0.95

    for (i in seq_len(nrow(ds_anchors))) {
      effort_plot <- effort_plot +
        annotate("rect",
                 xmin = ds_anchors$anchor_start[i],
                 xmax = ds_anchors$anchor_end[i],
                 ymin = -Inf, ymax = Inf,
                 fill = pal$anchor_fill, color = pal$anchor_border,
                 linewidth = 0.4, linetype = "dashed") +
        annotate("text",
                 x = ds_anchors$anchor_start[i] +
                   (ds_anchors$anchor_end[i] - ds_anchors$anchor_start[i]) / 2,
                 y = y_label,
                 label = paste0("Slice ", ds_anchors$slice[i]),
                 size = 2.8, color = pal$accent, fontface = "italic")
    }
  }

  effort_plot <- effort_plot +
    scale_x_date(date_labels = "%b %Y", expand = expansion(mult = 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title = paste0(dataset_name, " \u2014 Effort timeline"),
      subtitle = n_slices_label,
      x = NULL, y = "Active cameras"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, color = pal$accent),
      plot.subtitle = element_text(size = 9, color = "grey50"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.margin = margin(10, 12, 10, 12)
    )

  # ── Anchor summary table ────────────────────────────────────────
  tbl_theme_sm <- ttheme_minimal(
    core = list(
      fg_params = list(fontsize = 8, hjust = 0, x = 0.02),
      bg_params = list(fill = c(pal$bg, "white"), col = NA)
    ),
    colhead = list(
      fg_params = list(fontsize = 8, fontface = "bold", hjust = 0, x = 0.02),
      bg_params = list(fill = alpha(pal$accent, 0.12), col = NA)
    )
  )

  if (nrow(ds_anchors) > 0) {
    anchor_summary <- ds_anchors |>
      transmute(
        Slice = slice,
        Period = paste(format(anchor_start, "%Y-%m-%d"), "\u2013",
                       format(anchor_end, "%Y-%m-%d")),
        `Trap-days` = comma(as.integer(trap_days_12m)),
        Score = formatC(score, digits = 3, format = "f")
      )
    anchor_grob <- tableGrob(anchor_summary, rows = NULL, theme = tbl_theme_sm)

    # ── Per-slice protocol effort table ─────────────────────────────
    dep_parsed <- deploy |>
      transmute(
        location  = .data[[loc_col]],
        dep_start = as_date(parse_ts_safe(.data[[start_col]])),
        dep_end   = as_date(parse_ts_safe(.data[[end_col]]))
      ) |>
      filter(!is.na(dep_start), !is.na(dep_end), dep_end > dep_start)

    protocol_effort_rows <- purrr::map_dfr(seq_len(nrow(ds_anchors)), \(si) {
      s <- ds_anchors[si, ]
      s_start <- as.Date(s$anchor_start)
      s_end   <- as.Date(s$anchor_end)

      # Cameras active any time during the full anchor period
      cams_full <- dep_parsed |>
        filter(!(dep_end < s_start | dep_start > s_end)) |>
        pull(location) |> unique()

      prot_wins <- protocol_calendar_dates(s_start, s_end)

      purrr::map_dfr(seq_len(nrow(prot_wins)), \(j) {
        w <- prot_wins[j, ]
        active <- dep_parsed |>
          filter(!(dep_end < w$win_start | dep_start > w$win_end))
        active_cams <- unique(active$location)

        td <- active |>
          mutate(
            eff_start = pmax(dep_start, w$win_start),
            eff_end   = pmin(dep_end, w$win_end),
            days      = as.integer(eff_end - eff_start) + 1L
          ) |>
          summarise(trap_days = sum(pmax(days, 0L)))

        n_full <- max(length(cams_full), 1L)
        tibble(
          Slice     = s$slice,
          Window    = gsub("SNAP_EU_", "", w$window),
          Period    = paste(format(w$win_start, "%b %d"), "\u2013",
                           format(w$win_end, "%b %d %Y")),
          `Trap-d`  = comma(td$trap_days),
          Sites     = length(active_cams),
          `Full`    = length(cams_full),
          `Cov.`    = paste0(round(100 * length(active_cams) / n_full), "%")
        )
      })
    })

    prot_grob <- tableGrob(protocol_effort_rows, rows = NULL, theme = tbl_theme_sm)

    page2 <- arrangeGrob(
      ggplotGrob(effort_plot),
      arrangeGrob(
        textGrob("Anchor slices", gp = gpar(fontface = "bold", cex = 0.8, col = "grey40"),
                 just = "left", x = 0.02),
        anchor_grob,
        textGrob("Protocol window effort (Sites Full = cameras in 12-mo anchor; Cov. = window / full)",
                 gp = gpar(cex = 0.65, col = "grey50"), just = "left", x = 0.02),
        prot_grob,
        ncol = 1, heights = unit(c(0.08, 0.35, 0.08, 0.49), "npc")
      ),
      heights = unit(c(0.52, 0.48), "npc"),
      padding = unit(0.5, "line")
    )
  } else {
    page2 <- ggplotGrob(effort_plot)
  }


  # ================================================================
  # PAGE 3 — SPECIES × MONTH HEATMAP
  # ================================================================

  obs_det <- obs_ts |>
    filter(
      (.empty == FALSE | is.na(.empty)),
      !is.na(species),
      !str_starts(species, regex("Undefined|Homo\\s+sapiens|Other", ignore_case = TRUE)),
      species %in% keep_species
    )

  month_min <- floor_date(start_date, "month")
  month_max <- floor_date(end_date, "month")
  all_months <- seq.Date(month_min, month_max, by = "month")

  effort_month <- monthly_deployment_effort(daily_dep)

  dets_month <- obs_det |>
    mutate(month_date = as.Date(floor_date(timestamp, "month"))) |>
    filter(
      timestamp >= as.POSIXct(start_date, tz = "UTC"),
      timestamp <= as.POSIXct(end_date + days(1), tz = "UTC")
    ) |>
    count(species, month_date, name = "n_det") |>
    complete(species, month_date = all_months, fill = list(n_det = 0L)) |>
    filter(!is.na(species))

  sp_totals <- dets_month |>
    summarise(n_tot = sum(n_det), .by = species) |>
    filter(n_tot > 0)

  if (nrow(sp_totals) > 0) {
    heatmap_data <- dets_month |>
      inner_join(sp_totals, by = "species") |>
      left_join(effort_month, by = "month_date") |>
      mutate(
        trap_nights = replace_na(trap_nights, 0L),
        det_rate = if_else(trap_nights > 0, n_det / trap_nights * 100, 0)
      )

    sp_order <- sp_totals |> arrange(n_tot) |> pull(species)
    heatmap_data <- heatmap_data |>
      mutate(species = factor(species, levels = sp_order))

    n_breaks <- max(1, length(all_months) %/% 18)
    heatmap_plot <- ggplot(heatmap_data, aes(month_date, species, fill = det_rate)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_viridis_c(
        option = "inferno", direction = -1,
        name = "Detections\nper 100\ntrap-nights",
        na.value = "grey95"
      ) +
      scale_x_date(date_labels = "%b\n%Y", expand = c(0, 0),
                   breaks = all_months[seq(1, length(all_months), by = n_breaks)]) +
      labs(
        title = paste0(dataset_name, " \u2014 Species detection rates"),
        subtitle = "Detections per 100 trap-nights by month",
        x = NULL, y = NULL
      ) +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 13, color = pal$accent),
        plot.subtitle = element_text(size = 9, color = "grey50"),
        axis.text.y = element_text(face = "italic", size = 7),
        axis.text.x = element_text(size = 7),
        panel.grid = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        plot.margin = margin(10, 10, 10, 10)
      )

    # ── Overlay protocol window boundaries ──────────────────────────
    if (nrow(ds_anchors) > 0) {
      all_prot_wins <- purrr::map_dfr(seq_len(nrow(ds_anchors)), \(i) {
        protocol_calendar_dates(ds_anchors$anchor_start[i], ds_anchors$anchor_end[i])
      }) |> distinct()

      win_border <- c(SNAP_EU_CORE = pal$core_border, SNAP_EU_BUFFER = pal$buffer_border,
                      EOW_EARLY = pal$eow_e_border, EOW_LATE = pal$eow_l_border)
      win_lty    <- c(SNAP_EU_CORE = "solid", SNAP_EU_BUFFER = "dashed",
                      EOW_EARLY = "dotted", EOW_LATE = "dotted")

      for (k in seq_len(nrow(all_prot_wins))) {
        w <- all_prot_wins[k, ]
        heatmap_plot <- heatmap_plot +
          annotate("rect",
                   xmin = w$win_start, xmax = w$win_end,
                   ymin = -Inf, ymax = Inf,
                   fill = NA,
                   color = win_border[w$window] %||% "grey50",
                   linewidth = 0.7,
                   linetype = win_lty[w$window] %||% "solid")
      }

      # Legend annotation at the bottom-right
      heatmap_plot <- heatmap_plot +
        labs(subtitle = paste0(
          "Detections per 100 trap-nights by month  |  ",
          "CORE (\u2500 red)  BUFFER (\u2504 blue)  EOW-E (\u2508 green)  EOW-L (\u2508 purple)"
        ))
    }

  } else {
    heatmap_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No detections for kept species.",
               size = 5, color = "grey50") +
      theme_void()
  }


  # ================================================================
  # PAGES 4+ — PER-SLICE: RICHNESS MAP + ACCUMULATION CURVES
  # ================================================================

  slice_pages <- list()

  for (si in seq_len(nrow(ds_anchors))) {
    s <- ds_anchors[si, ]
    s_start <- as.Date(s$anchor_start)
    s_end   <- as.Date(s$anchor_end)
    slice_lab <- paste0("Slice ", s$slice, " (", s_start, " \u2192 ", s_end, ")")

    # -- Richness map --
    deploy_ok <- deploy |>
      transmute(
        deploymentID = .data[[depl_id_deploy]],
        start = parse_ts_safe(.data[[start_col]]),
        end   = parse_ts_safe(.data[[end_col]]),
        lat   = parse_num_safe(.data[[lat_col]]),
        lon   = parse_num_safe(.data[[lon_col]])
      ) |>
      filter(!is.na(start), !is.na(end), end > start) |>
      mutate(d_start = as_date(start), d_end = as_date(end)) |>
      filter(!(d_end < s_start | d_start > s_end))

    obs_slice <- obs_ts |>
      filter(
        (.empty == FALSE | is.na(.empty)),
        !is.na(species), nzchar(species),
        species %in% keep_species,
        as_date(timestamp) >= s_start,
        as_date(timestamp) <= s_end
      )

    det_total <- nrow(obs_slice)
    det_days  <- n_distinct(as_date(obs_slice$timestamp))

    daily_dep_slice <- daily_dep_full |> filter(date >= s_start, date <= s_end)
    trap_nights_slice <- sum(daily_dep_slice$n_active_deployments, na.rm = TRUE)

    subtitle_stats <- paste0(
      "trap-days = ", comma(trap_nights_slice),
      " | detections = ", comma(det_total),
      " | det-days = ", comma(det_days)
    )

    cam_rich <- obs_slice |>
      mutate(deploymentID = as.character(deploymentID)) |>
      distinct(deploymentID, species) |>
      count(deploymentID, name = "richness")

    cams_sf <- deploy_ok |>
      mutate(deploymentID = as.character(deploymentID)) |>
      filter(is.finite(lon), is.finite(lat)) |>
      left_join(cam_rich, by = "deploymentID") |>
      mutate(richness = replace_na(richness, 0L)) |>
      st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

    if (nrow(cams_sf) > 0 && !is.null(world_sf)) {
      s2_was_on <- sf::sf_use_s2()
      sf::sf_use_s2(FALSE)
      world_crop <- tryCatch(
        sf::st_intersection(world_sf, st_as_sfc(bb_pad)),
        error = function(e) world_sf
      )
      sf::sf_use_s2(s2_was_on)

      p_richness <- ggplot() +
        geom_sf(data = world_crop, fill = "grey98", colour = "grey70", linewidth = 0.2) +
        geom_point(data = cams_sf, aes(x = lon, y = lat, colour = richness),
                   size = 2.2, alpha = 0.95) +
        scale_colour_viridis_c(name = "Species\nrichness") +
        coord_sf(xlim = c(bb_pad["xmin"], bb_pad["xmax"]),
                 ylim = c(bb_pad["ymin"], bb_pad["ymax"]), expand = FALSE) +
        annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.7) +
        labs(
          title = paste0("Station richness \u2014 ", slice_lab),
          subtitle = subtitle_stats,
          x = NULL, y = NULL
        ) +
        theme_minimal(base_size = 11) +
        theme(
          plot.title = element_text(face = "bold", size = 12, color = pal$accent),
          plot.subtitle = element_text(size = 9, color = "grey50")
        )
    } else {
      p_richness <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0("No valid coordinates — ", slice_lab),
                 size = 4, color = "grey50") +
        theme_void()
    }
    slice_pages <- c(slice_pages, list(p_richness))

    # -- Accumulation curve (chronological) --
    month_start_s <- floor_date(s_start, "month")
    month_end_s   <- floor_date(s_end, "month")
    months_sel    <- tibble(month_date = seq.Date(month_start_s, month_end_s, by = "month"))
    nM <- nrow(months_sel)

    month_labels <- if (nM <= 12 && n_distinct(year(months_sel$month_date)) == 1) {
      month.abb[month(months_sel$month_date)]
    } else {
      paste0(month.abb[month(months_sel$month_date)], " ",
             substr(year(months_sel$month_date), 3, 4))
    }

    obs_m <- obs_slice |>
      transmute(species, month_date = as.Date(floor_date(timestamp, "month"))) |>
      semi_join(months_sel, by = "month_date") |>
      group_by(month_date) |>
      summarise(species_list = list(unique(species)), .groups = "drop") |>
      right_join(months_sel, by = "month_date") |>
      mutate(species_list = coalesce(species_list, list(character()))) |>
      arrange(month_date)

    effort_month_s <- monthly_deployment_effort(daily_dep_slice) |>
      right_join(months_sel, by = "month_date") |>
      mutate(trap_nights = replace_na(trap_nights, 0L)) |>
      arrange(month_date)

    # Single chronological curve
    seen <- character(0)
    acc_data <- tibble(
      step     = seq_len(nM),
      richness = purrr::map_int(seq_len(nM), \(j) {
        seen <<- union(seen, obs_m$species_list[[j]])
        length(seen)
      }),
      new_sp   = c(obs_m$species_list[[1]] |> length(),
                   purrr::map_int(2:nM, \(j) {
                     length(setdiff(obs_m$species_list[[j]],
                                    Reduce(union, obs_m$species_list[1:(j-1)],
                                           accumulate = FALSE)))
                   })),
      trap_nights = effort_month_s$trap_nights
    )

    max_rich <- max(acc_data$richness, na.rm = TRUE)
    max_eff  <- max(acc_data$trap_nights, na.rm = TRUE)
    scale_fac <- if (is.finite(max_eff) && max_eff > 0) max_rich / max_eff else 1

    p_acc <- ggplot(acc_data, aes(step)) +
      geom_col(aes(y = new_sp), fill = alpha(pal$accent, 0.3),
               width = 0.6) +
      geom_line(aes(y = richness), colour = "black", linewidth = 1.1) +
      geom_point(aes(y = richness), colour = "black", size = 1.5) +
      geom_line(aes(y = trap_nights * scale_fac), colour = "firebrick",
                linetype = "dashed", linewidth = 0.7) +
      scale_x_continuous(breaks = seq_len(nM), labels = month_labels) +
      scale_y_continuous(
        name = "Cumulative species (bars = new spp.)",
        sec.axis = sec_axis(~ . / scale_fac, name = "Monthly trap-nights")
      ) +
      labs(
        title = paste0("Accumulation \u2014 ", slice_lab),
        subtitle = subtitle_stats,
        x = NULL
      ) +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 12, color = pal$accent),
        plot.subtitle = element_text(size = 9, color = "grey50"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
        plot.margin = margin(8, 10, 6, 10)
      )

    slice_pages <- c(slice_pages, list(p_acc))
  }


  # ================================================================
  # INDIVIDUAL SPECIES PANELS — top 5 + bottom 5 by detections
  # ================================================================

  mamm_rate <- dets_month |>
    inner_join(sp_totals, by = "species") |>
    left_join(effort_month, by = "month_date") |>
    mutate(trap_nights = replace_na(trap_nights, 0L)) |>
    mutate(
      ymax  = max(n_det, na.rm = TRUE),
      emax  = max(trap_nights, na.rm = TRUE),
      scale = if_else(emax > 0, ymax / emax, 1),
      .by = species
    ) |>
    mutate(
      species_label = paste0(species, "  (n = ", comma(n_tot), ")"),
      effort_scaled = trap_nights * scale
    )

  if (nrow(mamm_rate) > 0) {
    # Keep top 5 + bottom 5 (or all if ≤10 species)
    sp_ranked <- sp_totals |> arrange(desc(n_tot))
    n_sp <- nrow(sp_ranked)
    if (n_sp > 10) {
      focal_sp <- c(head(sp_ranked$species, 5), tail(sp_ranked$species, 5))
    } else {
      focal_sp <- sp_ranked$species
    }

    make_sp_panel <- function(d, lab) {
      sc <- d$scale[1]
      ggplot(d, aes(month_date)) +
        geom_col(aes(y = n_det), fill = alpha(pal$accent, 0.6),
                 color = NA, width = 25) +
        geom_line(aes(y = effort_scaled), linewidth = 0.5,
                  color = pal$accent2, alpha = 0.6) +
        scale_y_continuous(
          name = "Detections",
          sec.axis = sec_axis(\(x) x / sc, name = "Trap-nights"),
          expand = expansion(mult = c(0, 0.08))
        ) +
        scale_x_date(date_breaks = "2 months", date_labels = "%b %Y",
                     expand = expansion(mult = 0.01)) +
        labs(x = NULL, title = lab) +
        theme_minimal(base_size = 10) +
        theme(
          plot.title = element_text(face = "italic", size = 10, color = pal$accent),
          axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
          axis.title.y.left = element_text(size = 8, color = pal$accent),
          axis.title.y.right = element_text(size = 8, color = pal$accent2),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.margin = margin(6, 10, 4, 10)
        )
    }

    plot_list <- mamm_rate |>
      filter(species %in% focal_sp) |>
      nest(.by = c(species, species_label)) |>
      arrange(desc(map_int(data, \(d) sum(d$n_det)))) |>
      mutate(plot = pmap(list(data, species_label), make_sp_panel)) |>
      pull(plot)
  } else {
    plot_list <- list()
  }


  # ================================================================
  # REJECTED SPECIES TABLE
  # ================================================================

  rejected_tbl <- obs_span |>
    mutate(species = if_else(is.na(species) | !nzchar(species),
                             "<missing>", as.character(species))) |>
    filter(!species %in% keep_species) |>
    count(species, sort = TRUE, name = "n_rejected")

  make_rejected_pages <- function(tbl, chunk = 50L) {
    if (!nrow(tbl))
      return(list(
        arrangeGrob(
          textGrob("No rejected observations.",
                   gp = gpar(cex = 1.1, col = "grey50")),
          heights = unit(1, "npc")
        )
      ))
    idx <- split(seq_len(nrow(tbl)), ceiling(seq_len(nrow(tbl)) / chunk))
    lapply(seq_along(idx), \(i) {
      g <- tableGrob(tbl[idx[[i]], ], rows = NULL,
                     theme = ttheme_minimal(
                       core = list(
                         fg_params = list(fontsize = 8, hjust = 0, x = 0.02),
                         bg_params = list(fill = c(pal$bg, "white"), col = NA)
                       ),
                       colhead = list(
                         fg_params = list(fontsize = 9, fontface = "bold",
                                          hjust = 0, x = 0.02),
                         bg_params = list(fill = alpha(pal$accent, 0.12), col = NA)
                       )
                     ))
      arrangeGrob(
        textGrob(
          sprintf("Rejected species (page %d/%d)", i, length(idx)),
          gp = gpar(fontface = "bold", cex = 1.0, col = "grey40")
        ),
        g, heights = c(0.06, 0.94)
      )
    })
  }
  rejected_pages <- make_rejected_pages(rejected_tbl)


  # ================================================================
  # WRITE PDF
  # ================================================================

  if (is.null(output_pdf)) {
    out_dir <- dataset_output_dir(archive_path)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    output_pdf <- file.path(out_dir, paste0(dataset_name, "_report.pdf"))
  } else {
    dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  }
  output_pdf <- normalizePath(output_pdf, mustWork = FALSE)

  pdf(output_pdf, width = 8.27, height = 11.69, onefile = TRUE, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)

  # Page 1: Dashboard
  print(page1)

  # Page 2: Effort + Anchors
  grid.newpage()
  grid.draw(page2)

  # Page 3: Heatmap
  print(heatmap_plot)

  # Pages 4+: Per-slice (richness map + accumulation)
  for (pg in slice_pages) print(pg)

  # Species panels (3 per page; top 5 + bottom 5)
  if (length(plot_list)) {
    sp_header <- if (n_sp > 10) {
      paste0("Monthly detections — top 5 + bottom 5 of ", n_sp, " species")
    } else {
      "Monthly detections by species"
    }
    print(marrangeGrob(
      grobs = plot_list, nrow = 3, ncol = 1,
      top = textGrob(sp_header,
                     gp = gpar(fontface = "bold", cex = 1.0, col = pal$accent))
    ))
  }

  # Rejected species
  print(marrangeGrob(grobs = rejected_pages, nrow = 1, ncol = 1))

  message("OK: ", dataset_name)
  invisible(if (pdf_only) NULL else metrics_row)
}


# ====================================================================
# Run all datasets  (all modes parallelised)
# ====================================================================
#
# Usage — uncomment ONE mode block below, then run the whole section.
#
#   (A) Metrics only   — fast, no PDFs, saves dataset_metadata.csv
#   (B) PDF only       — regenerate report PDFs, no metrics update
#   (C) Both           — metrics + PDFs (full run)
# ====================================================================
root <- "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets"
dataset_names <- basename(list.dirs(root, recursive = FALSE, full.names = TRUE))
N_WORKERS <- min(parallel::detectCores() - 1L, length(dataset_names), 8L)

# Pre-compute shared objects so each worker doesn't repeat this work
dataset_meta_shared <- read_dataset_meta()
get_world_sf()

# Anchors (reuse from Full1.R if available)
if (!exists("anchors", inherits = TRUE)) {
  message("Computing anchors for all datasets...")
  ds_paths <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  anchors <- purrr::map_dfr(ds_paths, find_anchors)
}

# ── (A) Metrics only ─────────────────────────────────────────────────
run_fn <- \(nm) dataset_metrics(nm, keep_species = spp_keep,
                                anchors_tbl = anchors,
                                dataset_meta = dataset_meta_shared)

# ── (B) PDF only ─────────────────────────────────────────────────────
 run_fn <- \(nm) { dataset_report(nm, keep_species = spp_keep,
                                   anchors_tbl = anchors,
                                   dataset_meta = dataset_meta_shared,
                                   pdf_only = TRUE); NULL }
# ── (C) Both (metrics + PDF) ─────────────────────────────────────────
 run_fn <- \(nm) dataset_report(nm, keep_species = spp_keep,
                                 anchors_tbl = anchors,
                                 dataset_meta = dataset_meta_shared)

# ── Execute ──────────────────────────────────────────────────────────
plan(multisession, workers = N_WORKERS)
message("Running ", length(dataset_names), " datasets on ", N_WORKERS, " workers")

results <- future_map(dataset_names, \(nm) {
  source(.helpers_path, local = FALSE)
  tryCatch(run_fn(nm), error = function(e) {
    warning("FAILED: ", nm, " \u2014 ", conditionMessage(e), call. = FALSE)
    NULL
  })
}, .progress = TRUE, .options = furrr_options(
  seed = NULL,
  packages = c("tidyverse", "lubridate", "sf", "cowplot", "ggspatial",
               "scales", "grid", "geosphere", "gridExtra", "readxl")
))

plan(sequential)

# ── Collect metrics (modes A and C) ──────────────────────────────────
dataset_overview_metrics <- results |>
  purrr::compact() |>
  dplyr::bind_rows()

if (nrow(dataset_overview_metrics) > 0) {
  readr::write_csv(dataset_overview_metrics, "dataset_metadata.csv")
  saveRDS(dataset_overview_metrics, "dataset_overview_metrics.rds")
  message("Saved dataset_metadata.csv and dataset_overview_metrics.rds (",
          nrow(dataset_overview_metrics), " datasets)")
  dataset_overview_metrics
}

rm(results, dataset_meta_shared, run_fn)


# ====================================================================
# Map dataset centroids (Europe)
# ====================================================================
if (!exists("dataset_overview_metrics", envir = .GlobalEnv) &&
    file.exists("dataset_overview_metrics.rds")) {
  dataset_overview_metrics <- readRDS("dataset_overview_metrics.rds")
  message("Loaded dataset_overview_metrics.rds (", nrow(dataset_overview_metrics), " datasets)")
}

if (exists("dataset_overview_metrics", envir = .GlobalEnv)) {
  .centroids_sf <- dataset_overview_metrics |>
    transmute(dataset = dataset_name, lon = centroid_lon, lat = centroid_lat) |>
    filter(is.finite(lon), is.finite(lat)) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

  .eu_bbox <- st_bbox(c(xmin = -12, ymin = 34, xmax = 35, ymax = 72),
                      crs = st_crs(4326))

  .world_sf <- get_world_sf()
  .world_eu <- if (!is.null(.world_sf))
    suppressWarnings(st_intersection(.world_sf, st_as_sfc(.eu_bbox)))

  dataset_map <- ggplot() +
    { if (!is.null(.world_eu)) geom_sf(data = .world_eu, fill = "grey95",
                                        color = "grey75", linewidth = 0.2) } +
    geom_sf(data = .centroids_sf, size = 2.5, alpha = 0.9, color = "#2166ac") +
    ggrepel::geom_text_repel(
      data = st_drop_geometry(.centroids_sf),
      aes(x = lon, y = lat, label = dataset),
      size = 2.5, max.overlaps = 30, segment.size = 0.2,
      segment.color = "grey60", color = "grey30"
    ) +
    coord_sf(xlim = c(.eu_bbox["xmin"], .eu_bbox["xmax"]),
             ylim = c(.eu_bbox["ymin"], .eu_bbox["ymax"]), expand = FALSE) +
    labs(title = "Dataset centroids", x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 14, color = "#2166ac"))

  print(dataset_map)
  rm(.centroids_sf, .world_sf, .world_eu, .eu_bbox, dataset_map)
}

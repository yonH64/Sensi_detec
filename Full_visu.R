# ─────────────────────────────────────────────────────────
# Libraries (plotting & helpers)
# ─────────────────────────────────────────────────────────
library(tidyverse)
library(lubridate)
library(viridis)
library(glue)
library(GGally)
library(ggplot2)
library(dplyr)
library(tidyr)

# ─────────────────────────────────────────────────────────
# Metric labels for plotting
# ─────────────────────────────────────────────────────────
# Use these with as_labeller() in facet_wrap/facet_grid
metric_labels <- c(
  # Detection — values
  lambda       = "Daily detection rate (\u03bb)",
  log_rate     = "log(event rate)",
  matched_rate = "Matched-camera rate",
  spatial_cov  = "Spatial coverage",
  
  # Detection — offsets (deltas vs FULL)
  d_lambda       = "\u0394 \u03bb",
  d_rate         = "\u0394 rate",
  d_matched_rate = "\u0394 Matched rate",
  
  # Detection — standard errors
  lambda_se       = "SE(\u03bb)",
  rate_se         = "SE(rate)",
  matched_rate_se = "SE(Matched rate)",
  
  # Detection — MSE
  mse_lambda       = "MSE(\u03bb)",
  mse_rate         = "MSE(rate)",
  mse_matched_rate = "MSE(Matched rate)",
  
  # Rank preservation
  rho_lambda   = "\u03c1 \u03bb",
  rho_log_rate = "\u03c1 log(rate)",
  
  # Richness — offsets
  d_sr_raref       = "\u0394 Rarefied richness",
  prop_sr_full     = "Proportion of FULL spp detected",
  
  # Richness — SE
  sr_raref_se      = "SE(Rarefied richness)",

  # Richness — MSE
  mse_sr_raref     = "MSE(Rarefied richness)"
)


############################################################
#                   DROPPED WINDOWS CHECK                  #
############################################################

# Count by reason (a row can fail multiple reasons; we count each flag)
threshold_counts <- all_dropped_species %>%
  select(dataset, window_id, species,
         failed_min_events, failed_min_occasions_pos, failed_min_sites_pos) %>%
  pivot_longer(starts_with("failed_"),
               names_to = "reason", values_to = "failed") %>%
  filter(failed) %>%
  mutate(reason = recode(reason,
                         failed_min_events        = "Too few events (min_events)",
                         failed_min_occasions_pos = "Too few positive occasions",
                         failed_min_sites_pos     = "Too few positive sites")) %>%
  count(reason, sort = TRUE)

print(threshold_counts)



by_dataset <- all_dropped_species %>%
  select(dataset, window_id, species,
         failed_min_events, failed_min_occasions_pos, failed_min_sites_pos) %>%
  pivot_longer(starts_with("failed_"),
               names_to = "reason", values_to = "failed") %>%
  filter(failed) %>%
  mutate(reason = recode(reason,
                         failed_min_events        = "Too few events (min_events)",
                         failed_min_occasions_pos = "Too few positive occasions",
                         failed_min_sites_pos     = "Too few positive sites")) %>%
  count(dataset, reason, sort = TRUE)

print(by_dataset)



############################################################
#                  VISUALIZATION UTILITIES                 #
############################################################

# --------------------------------------------------------------------
# plot_calendar_heatmaps()
# --------------------------------------------------------------------
# Aim
#   Heat-map the window grid for detection metrics (species-level) or
#   richness metrics (dataset-level).
#
# domain = "detection" (default)
#   Loops over species (optionally merged across datasets).
#     view = "offset" -> d_lambda, d_rate, d_matched_rate
#     view = "se"     -> lambda_se, rate_se, matched_rate_se
#     view = "value"  -> lambda, rate, matched_rate
#     view = "mse"    -> mse_lambda, mse_rate, mse_matched_rate
#
# domain = "richness"
#   Dataset-level (no species column).
#     view = "offset" -> d_sr_raref, prop_sr_full
#     view = "se"     -> sr_raref_se (single facet; requires sr_raref_se column)
#     view = "value"  -> sr_raref, sr_obs
#     view = "mse"    -> mse_sr_raref (single facet)
#
# Arguments
#   df              : data frame -- all_window_species (detection) or
#                     all_window_richness (richness).
#   dropped_df      : optional; species-level dropped data (detection only).
#   domain          : "detection" | "richness".
#   merge_datasets  : TRUE = aggregate across datasets; FALSE = per-dataset.
#   agg             : "mean" | "median" for merged mode.
#   view            : which metric set to display (see above).
#   rng_offset      : half-range for diverging offset colour scale.
#   limits_se       : c(min, max) for SE colour scale.
#   limits_value    : c(min, max) for value scale (NULL = auto).
#   limits_mse      : c(min, max) for MSE scale (NULL = auto).
#   datasets        : optional character filter.
#   species         : optional character filter (detection only).
# --------------------------------------------------------------------

plot_calendar_heatmaps <- function(df,
                                   dropped_df      = NULL,
                                   domain          = c("detection", "richness"),
                                   merge_datasets  = TRUE,
                                   agg             = c("mean", "median"),
                                   view            = c("offset", "se", "value", "mse"),
                                   rng_offset      = 0.30,
                                   limits_se       = c(0, 0.02),
                                   limits_value    = NULL,
                                   limits_mse      = NULL,
                                   datasets        = NULL,
                                   species         = NULL) {

  suppressPackageStartupMessages({
    library(tidyverse); library(lubridate); library(glue)
  })

  domain <- match.arg(domain)
  agg    <- match.arg(agg)
  view   <- match.arg(view)

  # -- validate view x domain ------------------------------------------------
  valid_views <- list(
    detection = c("offset", "se", "value", "mse"),
    richness  = c("offset", "se", "value", "mse")
  )
  if (!view %in% valid_views[[domain]]) {
    stop(glue("view = '{view}' is not available for domain = '{domain}'.\n",
              "Valid views: {paste(valid_views[[domain]], collapse = ', ')}"))
  }

  # -- drop non-grid protocol windows ----------------------------------------
  is_grid <- function(x) stringr::str_detect(x, "^d\\d{3}_L\\d+$") | x == "FULL"
  df <- df |> dplyr::filter(is_grid(window_id))
  if (!is.null(dropped_df) && nrow(dropped_df))
    dropped_df <- dropped_df |> dplyr::filter(is_grid(window_id))

  # -- Snapshot Europe protocol windows ---------------------------------------
  snapshot_windows <- c("d246_L57", "d246_L64", "d239_L57", "d239_L64")
  snapshot_col     <- "#1a5fb4"

  # -- facet columns by domain x view ----------------------------------------
  facet_cols <- switch(
    paste(domain, view, sep = "."),
    detection.offset = c("d_lambda", "d_rate", "d_matched_rate"),
    detection.se     = c("lambda_se", "rate_se", "matched_rate_se"),
    detection.value  = c("lambda", "rate", "matched_rate"),
    detection.mse    = c("mse_lambda", "mse_rate", "mse_matched_rate"),
    richness.offset  = c("d_sr_raref", "prop_sr_full"),
    richness.se      = c("sr_raref_se"),
    richness.value   = c("sr_raref", "sr_obs"),
    richness.mse     = c("mse_sr_raref")
  )

  missing_cols <- setdiff(facet_cols, names(df))
  if (length(missing_cols))
    stop(glue("Column(s) not found in data: {paste(missing_cols, collapse = ', ')}.\n",
              "Run Full1.R to regenerate the data with the latest metric pipeline."))

  legend_title <- switch(view,
    offset = "Offset\n(window \u2013 full)",
    se     = "SE",
    value  = "Value",
    mse    = "MSE"
  )

  # -- smart default limits ---------------------------------------------------
  if (is.null(limits_value))
    limits_value <- if (domain == "richness") NULL else c(0, 1)
  if (is.null(limits_mse))
    limits_mse <- if (domain == "richness") NULL else c(0, 0.05)

  # -- helpers ----------------------------------------------------------------
  year0   <- as.Date("2021-01-01")
  agg_fun <- if (agg == "mean") \(x) mean(x, na.rm = TRUE)
             else               \(x) stats::median(x, na.rm = TRUE)

  if (!exists("anchors", inherits = TRUE))
    stop("`anchors` not found in the environment.")
  anchor_map <- anchors |>
    group_by(dataset) |>
    arrange(as.Date(anchor_start), .by_group = TRUE) |>
    mutate(
      dataset_plot_id = if (n() == 1L) dataset
                        else paste0(dataset, "_slice", row_number()),
      anchor_doy = yday(as.Date(anchor_start))
    ) |>
    ungroup() |>
    transmute(dataset = dataset_plot_id, anchor_doy)

  decode_anchor <- function(d) {
    d |>
      mutate(
        start_doy = readr::parse_integer(
          stringr::str_extract(window_id, "(?<=d)\\d{3}")),
        length_d  = readr::parse_integer(
          stringr::str_extract(window_id, "(?<=_L)\\d+"))
      ) |>
      left_join(anchor_map, by = "dataset") |>
      mutate(
        start_off  = ((start_doy - anchor_doy) %% 365) + 1L,
        start_date = year0 + as.integer(start_off - 1L),
        wraps      = (start_doy + length_d - 1L) > 365L
      )
  }

  tile_params <- function(df_grid) {
    ux <- sort(unique(df_grid$start_date))
    uy <- sort(unique(df_grid$length_d))
    x_step <- suppressWarnings(median(diff(ux), na.rm = TRUE))
    if (!is.finite(x_step) || x_step <= 0)
      x_step <- as.difftime(7, units = "days")
    y_step <- suppressWarnings(median(diff(uy), na.rm = TRUE))
    if (!is.finite(y_step) || y_step <= 0) y_step <- 7
    list(
      tile_w = as.numeric(x_step, units = "days") * 1.005,
      tile_h = y_step * 1.005
    )
  }

  # -- fill scale builder -----------------------------------------------------
  build_fill <- function(agg_prefix = NULL) {
    nm <- if (!is.null(agg_prefix)) glue("{agg_prefix} {view}") else legend_title
    switch(view,
      offset = scale_fill_gradient2(
        low = "#2166ac", mid = "grey98", high = "#b2182b",
        midpoint = 0, limits = c(-rng_offset, rng_offset),
        oob = scales::squish, name = nm),
      se = scale_fill_viridis_c(
        option = "magma", direction = -1,
        limits = limits_se, oob = scales::squish, name = nm),
      mse = scale_fill_viridis_c(
        option = "inferno", direction = -1,
        limits = limits_mse, oob = scales::squish, name = nm),
      scale_fill_viridis_c(
        option = "viridis",
        limits = limits_value, oob = scales::squish, name = nm)
    )
  }

  # -- shared theme -----------------------------------------------------------
  theme_heatmap <- function(base_size = 13) {
    theme_minimal(base_size = base_size) %+replace%
      theme(
        panel.grid        = element_blank(),
        panel.background  = element_rect(fill = "grey96", colour = NA),
        strip.text        = element_text(face = "bold", size = rel(0.95),
                                         hjust = 0),
        axis.text.x       = element_text(angle = 45, hjust = 1,
                                         size = rel(0.85)),
        axis.title        = element_text(size = rel(0.9)),
        plot.title        = element_text(face = "bold", size = rel(1.1),
                                         hjust = 0),
        plot.subtitle     = element_text(colour = "grey40", size = rel(0.8),
                                         hjust = 0),
        legend.key.height = unit(0.8, "cm"),
        legend.key.width  = unit(0.35, "cm"),
        plot.margin       = margin(8, 10, 8, 8)
      )
  }

  # -- facet labeller from metric_labels dict ---------------------------------
  make_labeller <- function(grid_df) {
    val_col <- if ("val_agg" %in% names(grid_df)) "val_agg" else "val"
    counts <- grid_df |>
      summarise(n_cells = sum(!is.na(.data[[val_col]])), .by = metric)
    labs_vec <- counts |>
      mutate(
        pretty = ifelse(metric %in% names(metric_labels),
                        metric_labels[metric], metric),
        label  = glue("{pretty}  ({n_cells} cells)")
      )
    as_labeller(setNames(labs_vec$label, labs_vec$metric))
  }

  # -- x / y scales -----------------------------------------------------------
  x_scale <- scale_x_date(
    breaks = seq(year0, year0 + 364, by = "1 month"),
    labels = scales::label_date("%b"),
    limits = c(year0, year0 + 364),
    expand = c(0, 0), name = "Window start"
  )
  y_scale <- scale_y_continuous(expand = c(0, 0),
                                name = "Window length (days)")

  # -- dropped-cell preprocessor (detection only, species-level) --------------
  prep_drops <- function(df_drop, sp_keep, ds_keep = NULL,
                         n_ds_tot = NULL) {
    if (is.null(df_drop) || !nrow(df_drop)) return(tibble())
    out <- df_drop |>
      filter(species %in% sp_keep) |>
      {\(d) if (!is.null(ds_keep)) filter(d, dataset %in% ds_keep) else d}() |>
      decode_anchor() |>
      distinct(species, dataset, start_off, length_d) |>
      group_by(species, start_off, length_d) |>
      summarise(n_drop = n_distinct(dataset), .groups = "drop") |>
      mutate(start_date = year0 + as.integer(start_off - 1L))
    if (!is.null(n_ds_tot))
      out <- out |> mutate(drop_class = case_when(
        n_drop >= n_ds_tot ~ "all", n_drop > 0 ~ "some", TRUE ~ "none"))
    out
  }

  # -- snapshot outline helper ------------------------------------------------
  get_snap <- function(base_df, drop_df = NULL, sp_f = NULL,
                       ds_f = NULL, use_doy = FALSE) {
    sb <- base_df |> filter(window_id %in% snapshot_windows)
    if (!is.null(ds_f)) sb <- sb |> filter(dataset %in% ds_f)
    if (use_doy) {
      sb <- sb |> distinct(window_id, start_doy, length_d) |>
        mutate(start_date = year0 + as.integer(start_doy - 1L)) |>
        select(start_date, length_d)
    } else {
      sb <- sb |> distinct(window_id, start_date, length_d) |>
        select(start_date, length_d)
    }
    sd_tbl <- tibble()
    if (!is.null(drop_df) && nrow(drop_df)) {
      dd <- drop_df |> filter(window_id %in% snapshot_windows)
      if (!is.null(ds_f)) dd <- dd |> filter(dataset %in% ds_f)
      if (!is.null(sp_f) && "species" %in% names(dd))
        dd <- dd |> filter(species == sp_f)
      if (nrow(dd)) {
        dd <- dd |> decode_anchor()
        if (use_doy) {
          sd_tbl <- dd |> distinct(window_id, start_doy, length_d) |>
            mutate(start_date = year0 + as.integer(start_doy - 1L)) |>
            select(start_date, length_d)
        } else {
          sd_tbl <- dd |> distinct(window_id, start_date, length_d) |>
            select(start_date, length_d)
        }
      }
    }
    bind_rows(sb, sd_tbl) |> distinct(start_date, length_d)
  }

  # -- core plot assembler ----------------------------------------------------
  build_heatmap <- function(grid_df, dropped_tiles, snap_tiles,
                            wrap_tiles, pars, title_str, subtitle_str,
                            agg_prefix = NULL) {
    val_col <- if ("val_agg" %in% names(grid_df)) "val_agg" else "val"
    p <- ggplot()

    # dropped tiles
    if (nrow(dropped_tiles) && "drop_class" %in% names(dropped_tiles)) {
      p <- p +
        geom_tile(data = filter(dropped_tiles, drop_class == "all"),
                  aes(start_date, length_d), fill = "grey55", colour = NA,
                  width = pars$tile_w, height = pars$tile_h) +
        geom_tile(data = filter(dropped_tiles, drop_class == "some"),
                  aes(start_date, length_d), fill = "grey80", colour = NA,
                  width = pars$tile_w, height = pars$tile_h)
    } else if (nrow(dropped_tiles)) {
      p <- p +
        geom_tile(data = dropped_tiles, aes(start_date, length_d),
                  fill = "grey75", colour = NA,
                  width = pars$tile_w, height = pars$tile_h)
    }

    # NA tiles
    na_tiles <- grid_df |>
      filter(is.na(.data[[val_col]])) |>
      distinct(metric, start_date, length_d)
    if (nrow(na_tiles))
      p <- p + geom_tile(data = na_tiles, aes(start_date, length_d),
                         fill = "#fde725", alpha = 0.5, colour = NA,
                         width = pars$tile_w, height = pars$tile_h)

    # value tiles
    p <- p +
      geom_tile(data = filter(grid_df, !is.na(.data[[val_col]])),
                aes(start_date, length_d, fill = .data[[val_col]]),
                colour = NA, width = pars$tile_w, height = pars$tile_h)

    # snapshot outlines
    if (nrow(snap_tiles))
      p <- p + geom_tile(data = snap_tiles, aes(start_date, length_d),
                         fill = NA, colour = snapshot_col, linewidth = 0.8,
                         width = pars$tile_w, height = pars$tile_h)

    # wrap marks
    if (nrow(wrap_tiles))
      p <- p + geom_point(data = wrap_tiles, aes(start_date, length_d),
                          shape = 4, size = 1, stroke = 0.2,
                          colour = "black")

    p + x_scale + y_scale + build_fill(agg_prefix) +
      facet_wrap(~ metric, ncol = 1, labeller = make_labeller(grid_df)) +
      labs(title = title_str, subtitle = subtitle_str) +
      theme_heatmap()
  }

  # ====================================================================
  #  DETECTION DOMAIN (species-level)
  # ====================================================================
  if (domain == "detection") {

    base <- df |>
      filter(window_id != "FULL") |>
      {\(d) if (!is.null(datasets)) filter(d, dataset %in% datasets) else d}() |>
      {\(d) if (!is.null(species))  filter(d, species  %in% species)  else d}() |>
      decode_anchor()

    sp_vec <- unique(base$species)
    if (!length(sp_vec)) {
      message("No rows to plot."); return(invisible(NULL))
    }

    for (sp in sp_vec) {
      base_sp <- filter(base, species == sp)

      if (!merge_datasets) {
        # -- per-dataset --
        for (ds in unique(base_sp$dataset)) {
          dat <- base_sp |>
            filter(dataset == ds) |>
            select(dataset, species, window_id, start_date, length_d,
                   wraps, all_of(facet_cols)) |>
            pivot_longer(all_of(facet_cols),
                         names_to = "metric", values_to = "val")

          dropped_ds <- prep_drops(dropped_df, sp, ds) |>
            distinct(start_date, length_d, .keep_all = TRUE) |>
            select(start_date, length_d, n_drop)

          snap_ds   <- get_snap(base_sp, dropped_df, sp, ds)
          pars      <- tile_params(dat |> distinct(start_date, length_d))
          wrap_t    <- dat |> filter(wraps) |>
            distinct(metric, start_date, length_d)
          n_cells   <- sum(!is.na(dat$val)) / length(facet_cols)

          p <- build_heatmap(
            dat, dropped_ds, snap_ds, wrap_t, pars,
            title_str    = glue("{str_to_title(view)} \u2014 {sp} \u2022 {ds}  ({round(n_cells)} cells)"),
            subtitle_str = glue("\"\u00d7\" = wrap; blue outline = Snapshot Europe")
          )
          print(p)
          if (ds != dplyr::last(unique(base_sp$dataset)))
            readline("Return for next dataset\u2026")
        }

      } else {
        # -- merged (aggregated) --
        dat_long <- base_sp |>
          select(dataset, species, window_id, start_doy, length_d,
                 all_of(facet_cols)) |>
          pivot_longer(all_of(facet_cols),
                       names_to = "metric", values_to = "val")

        grid <- dat_long |>
          mutate(start_date = year0 + as.integer(start_doy - 1L)) |>
          group_by(metric, start_doy, start_date, length_d) |>
          summarise(val_agg = agg_fun(val), .groups = "drop")

        n_ds <- n_distinct(dat_long$dataset)

        dropped_sp <- if (!is.null(dropped_df)) {
          dropped_df |>
            filter(species %in% sp,
                   dataset %in% unique(dat_long$dataset)) |>
            decode_anchor() |>
            distinct(species, dataset, start_doy, length_d) |>
            group_by(species, start_doy, length_d) |>
            summarise(n_drop = n_distinct(dataset), .groups = "drop") |>
            mutate(
              start_date = year0 + as.integer(start_doy - 1L),
              drop_class = case_when(
                n_drop >= n_ds ~ "all", n_drop > 0 ~ "some",
                TRUE ~ "none")
            ) |>
            distinct(start_date, length_d, .keep_all = TRUE) |>
            select(start_date, length_d, drop_class)
        } else tibble()

        snap_sp <- get_snap(base_sp, dropped_df, sp, use_doy = TRUE)
        pars    <- tile_params(grid |> distinct(start_date, length_d))

        p <- build_heatmap(
          grid, dropped_sp, snap_sp, tibble(), pars,
          title_str    = glue("{str_to_title(view)} \u2014 {sp}"),
          subtitle_str = glue("Merged across {n_ds} dataset",
                              "{ifelse(n_ds>1,'s','')}; ",
                              "blue outline = Snapshot Europe"),
          agg_prefix   = agg
        )
        print(p)
        if (sp != dplyr::last(sp_vec))
          readline("Return for next species\u2026")
      }
    }
  }

  # ====================================================================
  #  RICHNESS DOMAIN (dataset-level, no species column)
  # ====================================================================
  if (domain == "richness") {

    base <- df |>
      filter(window_id != "FULL") |>
      {\(d) if (!is.null(datasets)) filter(d, dataset %in% datasets) else d}() |>
      decode_anchor()

    ds_vec <- unique(base$dataset)
    if (!length(ds_vec)) {
      message("No rows to plot."); return(invisible(NULL))
    }

    if (!merge_datasets) {
      # -- per-dataset --
      for (ds in ds_vec) {
        dat <- base |>
          filter(dataset == ds) |>
          select(dataset, window_id, start_date, length_d, wraps,
                 all_of(facet_cols)) |>
          pivot_longer(all_of(facet_cols),
                       names_to = "metric", values_to = "val")

        snap_ds  <- get_snap(base |> filter(dataset == ds), ds_f = ds)
        pars     <- tile_params(dat |> distinct(start_date, length_d))
        wrap_t   <- dat |> filter(wraps) |>
          distinct(metric, start_date, length_d)
        n_cells  <- sum(!is.na(dat$val)) / max(length(facet_cols), 1)

        p <- build_heatmap(
          dat, tibble(), snap_ds, wrap_t, pars,
          title_str    = glue("{str_to_title(view)} (richness) \u2014 {ds}  ({round(n_cells)} cells)"),
          subtitle_str = glue("\"\u00d7\" = wrap; blue outline = Snapshot Europe")
        )
        print(p)
        if (ds != dplyr::last(ds_vec))
          readline("Return for next dataset\u2026")
      }

    } else {
      # -- merged (aggregated) --
      dat_long <- base |>
        select(dataset, window_id, start_doy, length_d,
               all_of(facet_cols)) |>
        pivot_longer(all_of(facet_cols),
                     names_to = "metric", values_to = "val")

      grid <- dat_long |>
        mutate(start_date = year0 + as.integer(start_doy - 1L)) |>
        group_by(metric, start_doy, start_date, length_d) |>
        summarise(val_agg = agg_fun(val), .groups = "drop")

      n_ds <- n_distinct(dat_long$dataset)

      snap_rich <- get_snap(base, use_doy = TRUE)
      pars      <- tile_params(grid |> distinct(start_date, length_d))

      p <- build_heatmap(
        grid, tibble(), snap_rich, tibble(), pars,
        title_str    = glue("{str_to_title(view)} (richness)"),
        subtitle_str = glue("Merged across {n_ds} dataset",
                            "{ifelse(n_ds>1,'s','')}; ",
                            "blue outline = Snapshot Europe"),
        agg_prefix   = agg
      )
      print(p)
    }
  }

  invisible(NULL)
}


# -- Example calls -------------------------------------------------------

# Detection: offsets, per dataset
plot_calendar_heatmaps(all_window_species,
                       dropped_df = all_dropped_species,
                       domain = "detection",
                       merge_datasets = FALSE,
                       view = "offset", rng_offset = 0.30)

# Detection: offsets, merged
plot_calendar_heatmaps(all_window_species,
                       dropped_df = all_dropped_species,
                       domain = "detection",
                       merge_datasets = TRUE,
                       view = "offset", rng_offset = 0.30)

# Detection: SEs, merged
plot_calendar_heatmaps(all_window_species,
                       dropped_df = all_dropped_species,
                       domain = "detection",
                       merge_datasets = TRUE,
                       view = "se")

# Detection: MSE, merged
plot_calendar_heatmaps(all_window_species,
                       dropped_df = all_dropped_species,
                       domain = "detection",
                       merge_datasets = TRUE,
                       view = "mse")

# Detection: per dataset, single dataset
plot_calendar_heatmaps(all_window_species,
                       dropped_df = all_dropped_species,
                       domain = "detection",
                       merge_datasets = FALSE,
                       view = "offset",
                       datasets = c("BFNP_201819"))

# Richness: offsets (d_sr_raref, prop_sr_full), merged
plot_calendar_heatmaps(all_window_richness,
                       domain = "richness",
                       merge_datasets = TRUE,
                       view = "offset", rng_offset = 3)

# Richness: SE (sr_raref_se), merged
plot_calendar_heatmaps(all_window_richness,
                       domain = "richness",
                       merge_datasets = TRUE,
                       view = "se")

# Richness: MSE (mse_sr_raref), merged
plot_calendar_heatmaps(all_window_richness,
                       domain = "richness",
                       merge_datasets = TRUE,
                       view = "mse")

# Richness: per dataset
plot_calendar_heatmaps(all_window_richness,
                       domain = "richness",
                       merge_datasets = FALSE,
                       view = "offset", rng_offset = 3)


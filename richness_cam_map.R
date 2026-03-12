richness_map <- function(
    archive_path,
    ds,                 # dataset name present in `anchors$dataset`
    anchors,            # tibble with columns: dataset, anchor_start, anchor_end
    output_file = NULL,
    date_min = NULL,    # optional hard bound on deployments (YYYY-MM-DD)
    date_max = NULL     # optional hard bound on deployments
) {
  suppressPackageStartupMessages({
    library(tidyverse); library(lubridate); library(sf)
    library(rnaturalearth); library(rnaturalearthdata); library(viridis)
  })
  
  # -------- output path
  if (is.null(output_file)) {
    output_file <- file.path(
      getwd(),
      paste0(ds, "_camera_richness_map_12m.pdf")
    )
  }
  
  # -------- load data
  deploy <- read_csv_src(archive_path, "deployments\\.csv")
  obs    <- read_csv_src(archive_path, "observations\\.csv")
  
  # -------- column detection
  cols <- set_camtrap_cols(obs, deploy)
  list2env(cols, envir = environment())

  if (any(is.na(c(sci_col, date_col, depl_id_obs, depl_id_deploy, lat_col, lon_col, start_col, end_col))))
    stop("Missing required columns (species/date/deployment/lat/lon/start/end).")
  
  # -------- parse deployment times & apply optional hard bounds
  deploy_ok <- deploy %>%
    rename(start = !!sym(start_col), end = !!sym(end_col)) %>%
    mutate(
      start = ymd_hms(start, tz = "UTC"),
      end   = ymd_hms(end,   tz = "UTC")
    ) %>%
    filter(end > start, end - start <= years(10))
  
  if (!is.null(date_min)) deploy_ok <- deploy_ok %>% mutate(start = pmax(start, ymd(date_min, tz = "UTC")))
  if (!is.null(date_max)) deploy_ok <- deploy_ok %>% mutate(end   = pmin(end,   ymd(date_max, tz = "UTC")))
  
  # -------- helper: empty events
  clean_empty <- function(df) {
    if (!"is_empty" %in% names(df)) return(rep(FALSE, nrow(df)))
    v <- df$is_empty
    if (is.list(v)) purrr::map_lgl(v, ~ identical(.x, TRUE)) else as.logical(v)
  }
  
  # -------- world once
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>% st_make_valid()
  
  # -------- bbox helpers
  bbox_from_points <- function(pts, expand = 0.08) {
    bb <- sf::st_bbox(pts)
    if (any(!is.finite(bb))) {
      xy <- sf::st_coordinates(pts)
      rx <- range(xy[,1], na.rm = TRUE); ry <- range(xy[,2], na.rm = TRUE)
      bb <- c(xmin = rx[1], ymin = ry[1], xmax = rx[2], ymax = ry[2])
    }
    dx <- max(bb["xmax"] - bb["xmin"], 1e-4)
    dy <- max(bb["ymax"] - bb["ymin"], 1e-4)
    bb["xmin"] <- bb["xmin"] - dx*expand
    bb["xmax"] <- bb["xmax"] + dx*expand
    bb["ymin"] <- bb["ymin"] - dy*expand
    bb["ymax"] <- bb["ymax"] + dy*expand
    bb
  }
  make_crop_box <- function(bb, crs = 4326) {
    bb2 <- sf::st_bbox(
      c(xmin = as.numeric(bb["xmin"]),
        ymin = as.numeric(bb["ymin"]),
        xmax = as.numeric(bb["xmax"]),
        ymax = as.numeric(bb["ymax"])),
      crs = sf::st_crs(crs)
    )
    sf::st_as_sfc(bb2)
  }
  
  # -------- all slices for this dataset
  slices <- anchors %>%
    filter(dataset == ds) %>%
    arrange(anchor_start) %>%
    mutate(
      anchor_start = as.Date(anchor_start),
      anchor_end   = as.Date(anchor_end),
      slice_id     = row_number(),
      slice_lab    = paste0("slice ", slice_id, " (", anchor_start, " \u2192 ", anchor_end, ")")
    )
  
  if (nrow(slices) == 0) stop("No slices found in `anchors` for dataset: ", ds)
  
  # -------- inner worker: one slice -> one ggplot
  plot_slice <- function(start_date, end_date, slice_lab) {
    # filter to slice (DATE only – bypasses 00:00:00 parsing issues)
    obs_det <- obs %>%
      mutate(.empty = clean_empty(.)) %>%
      filter(.empty == FALSE | is.na(.empty)) %>%
      select(-.empty) %>%
      filter(!is.na(.data[[sci_col]])) %>%
      filter(!str_detect(.data[[sci_col]], taxa_filter)) %>%
      mutate(event_date = as_date(.data[[date_col]])) %>%
      filter(event_date >= start_date, event_date <= end_date) %>%
      transmute(
        deploymentID = .data[[depl_id_obs]],
        species      = .data[[sci_col]]
      )
    
    # deployments active in the slice (any overlap)
    dep_i <- deploy_ok %>%
      mutate(d_start = as_date(start), d_end = as_date(end)) %>%
      filter(!(d_end < start_date | d_start > end_date))
    
    # richness per deployment for this slice
    cam_rich <- obs_det %>%
      distinct(deploymentID, species) %>%
      count(deploymentID, name = "richness")
    
    # camera points with richness
    cams <- dep_i %>%
      transmute(
        deploymentID = .data[[depl_id_deploy]],
        lon = suppressWarnings(as.numeric(.data[[lon_col]])),
        lat = suppressWarnings(as.numeric(.data[[lat_col]]))
      ) %>%
      filter(is.finite(lon), is.finite(lat)) %>%
      left_join(cam_rich, by = "deploymentID") %>%
      mutate(richness = tidyr::replace_na(richness, 0L)) %>%
      sf::st_as_sf(coords = c("lon","lat"), crs = 4326, remove = FALSE)
    
    # bbox + crop
    bb <- bbox_from_points(cams, expand = 0.08)
    crop_box <- make_crop_box(bb, sf::st_crs(world))
    s2_was_on <- sf::sf_use_s2()
    sf::sf_use_s2(FALSE)
    world_crop <- tryCatch(sf::st_intersection(world, crop_box), error = function(e) world)
    sf::sf_use_s2(s2_was_on)
    
    ggplot() +
      geom_sf(data = world_crop, fill = "grey98", colour = "grey70", linewidth = 0.2) +
      geom_point(data = cams, aes(x = lon, y = lat, colour = richness), size = 2.2, alpha = 0.95) +
      scale_colour_viridis_c(name = "Species richness") +
      coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"])) +
      labs(
        title = "Camera species richness",
        subtitle = paste(ds, "\u2014", slice_lab),
        x = NULL, y = NULL
      ) +
      theme_minimal(base_size = 11)
  }
  
  # -------- write ONE multi-page PDF with all slices
  pdf(output_file, width = 8.3, height = 5.8, onefile = TRUE)  # ~A5 landscape
  on.exit(dev.off(), add = TRUE)
  
  for (i in seq_len(nrow(slices))) {
    s <- slices[i, ]
    p <- plot_slice(s$anchor_start, s$anchor_end, s$slice_lab)
    print(p)
  }
  
  message("✅ Saved multi-page PDF: ", output_file)
  invisible(output_file)
}

# anchors must exist (cols: dataset, anchor_start, anchor_end)
# e.g. to make one PDF for dataset "MICA":
richness_map(
  archive_path = "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/MICA",
  ds = "MICA",
  anchors = anchors
)

richness_map(
  archive_path = "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/SI-pirc",
  ds = "SI-pirc",
  anchors = anchors
)


# ------------------------------------------------------------
# dataset_wrapper()
#
# • One dataset folder/archive in  →  one named list out
# • Relies on two objects that must already exist in memory:
#     1.  window_tbl  – template of DOY windows (from make_window_template)
#     2.  anchors     – tibble with cols: dataset, anchor_start, anchor_end
# ------------------------------------------------------------
dataset_wrapper <- function(
    archive_path,
    independence_mins  = 30,
    occ_grain_days     = 1,
    min_events         = 20,
    min_occasions_pos  = 5,
    min_sites_pos      = 5)
 {
  suppressPackageStartupMessages({ library(tidyverse); library(lubridate); library(geosphere) })
  
  ds_name <- tools::file_path_sans_ext(basename(archive_path))
  print(ds_name)

  # ---- 1. load raw tables -------------------------------------------
  deploy <- read_csv_src(archive_path, "deployments\\.csv")
  obs    <- read_csv_src(archive_path, "observations\\.csv")

  # -----2. column detection
  cols <- set_camtrap_cols(obs, deploy)
  list2env(cols, envir = environment())
  
  if (any(is.na(c(sci_col, date_col, depl_id_obs, depl_id_deploy, lat_col, lon_col, start_col, end_col))))
    stop("Missing required columns (species/date/deployment/lat/lon/start/end).")
  
  # ---- 3. crop to dataset-specific 12-month slice(s) ---------------
  anchor_rows <- anchors %>% filter(dataset == ds_name)
  
  if (nrow(anchor_rows) == 0)
    stop("Anchor start/end not found for dataset: ", ds_name)
  
  a_start = anchor_rows$anchor_start
  a_end = anchor_rows$anchor_end
  slice_id = seq_len(nrow(anchor_rows))
  
  slice_runner <- function(a_start, a_end, slice_id) {
    
    deploy_tbl <- deploy %>% 
      transmute(camera_id = .data[[depl_id_deploy]],                
                start     = ymd_hms(.data[[start_col]], tz = "UTC"),
                end       = ymd_hms(.data[[end_col]], tz = "UTC")) %>% 
      filter(end > start) %>% 
      mutate(start = pmax(start, a_start),
             end   = pmin(end,   a_end)) %>% 
      filter(end > start)
    
    deploy_tbl <- deploy_tbl %>% 
      mutate(start_virtual = as.numeric(difftime(start, a_start, units = "days")) + 1,
             end_virtual   = as.numeric(difftime(end,   a_start, units = "days")) + 1)
    
    # distance & centroid helpers
    loc_tbl <- deploy %>%
      transmute(camera_id = .data[[depl_id_deploy]],     # use depl_id_deploy here
                lat       = .data[[lat_col]],
                lon       = .data[[lon_col]]) %>%
      distinct() %>%
      drop_na()
    
    dist_mat      <- geosphere::distm(loc_tbl[, c("lon", "lat")], fun = distHaversine) / 1000
    trap_array_km <- max(dist_mat)
    latitude_deg  <- geosphere::centroid(loc_tbl[, c("lon", "lat")])[2]
    
    # Bundle once
    array_summary <- tibble(trap_array = trap_array_km,
                            latitude   = latitude_deg)
    
    obs_tbl <- obs %>% 
      transmute(camera_id = .data[[depl_id_obs]],
                species   = .data[[sci_col]],
                timestamp = ymd_hms(.data[[date_col]], tz = "UTC")) %>% 
      filter(timestamp >= a_start & timestamp <= a_end,
             !str_detect(species, taxa_filter),
             !is.na(camera_id), !is.na(species), !is.na(timestamp)) %>% 
      mutate(doy_virtual = as.numeric(difftime(timestamp, a_start, units = "days")) + 1)
    
    window_grid_ds <- window_grid %>% 
      mutate(start_virtual = start_doy,
             end_virtual   = start_doy + length_d - 1L,
             start_date    = a_start + days(start_doy - 1),
             end_date      = start_date + days(length_d - 1),
             end_date      = if_else(end_virtual > 365, end_date - days(365), end_date))
    
    core <- build_window_metrics(
      window_grid        = window_grid_ds,
      deploy_tbl         = deploy_tbl,
      obs_tbl            = obs_tbl,
      a_start            = a_start,
      a_end              = a_end,
      independence_mins  = independence_mins,
      occ_grain_days     = occ_grain_days,
      min_events         = min_events,
      min_occasions_pos  = min_occasions_pos,
      min_sites_pos      = min_sites_pos,
      trap_array_km      = trap_array_km,
      latitude_deg       = latitude_deg,
      drop_leap_day      = TRUE)
    
    suffix <- if (nrow(anchor_rows) == 1) "" else paste0("_slice", slice_id)
    
    core$window_species  <- core$window_species  %>% mutate(dataset = paste0(ds_name, suffix))
    core$dropped_species <- core$dropped_species %>% mutate(dataset = paste0(ds_name, suffix))
    core$window_camera   <- core$window_camera   %>% mutate(dataset = paste0(ds_name, suffix))
    core
  }
  
  cores <- purrr::pmap(
    list(anchor_rows$anchor_start,
         anchor_rows$anchor_end,
         seq_len(nrow(anchor_rows))),
    slice_runner)
  
  if (length(cores) == 1) {
    core <- cores[[1]]
  } else {
    core <- list(
      window_species  = purrr::map_df(cores, "window_species"),
      dropped_species = purrr::map_df(cores, "dropped_species"),
      window_camera   = purrr::map_df(cores, "window_camera")
    )
  }
}


wrapped <- list.dirs("C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets", 
                     recursive = FALSE, full.names = TRUE) %>% 
  map(dataset_wrapper)

all_window_species  <- map_df(wrapped, "window_species")
saveRDS(all_window_species, "all_window_species.rds")
all_dropped_species <- map_df(wrapped, "dropped_species")
saveRDS(all_dropped_species, "all_dropped_species.rds")

all_window_species <- readRDS("all_window_species.rds")
all_dropped_species <- readRDS("all_dropped_species.rds")

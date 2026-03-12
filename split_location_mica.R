split_location_stats <- function(archive_path, cutoff_lat = 52) {
  library(tidyverse); library(lubridate)
  
  # -- robust timestamp parser (handles "YYYY-MM-DD" and "... UTC")
  parse_ts <- function(x) {
    out <- suppressWarnings(ymd_hms(x, tz = "UTC", quiet = TRUE))
    if (any(bad <- is.na(out))) {
      out[bad] <- ymd(str_remove(x[bad], "\\s+UTC$"), tz = "UTC")
    }
    out
  }
  
  # ---- load + detect columns
  deploy <- read_csv_src(archive_path, "deployments\\.csv")
  obs    <- read_csv_src(archive_path, "observations\\.csv")
  cols <- set_camtrap_cols(obs, deploy)
  list2env(cols, envir = environment())
  
  ds <- tools::file_path_sans_ext(basename(archive_path))
  
  # ----- find a location (site) column distinct from deployment id
  pick_location_col <- function(nms, avoid) {
    cand <- nms[str_detect(nms, regex("location|site|station", TRUE))]
    cand <- setdiff(cand, avoid)
    if (length(cand)) cand[1] else stop("No site/location column found in deployments.")
  }
  location_col <- pick_location_col(names(deploy), depl_id_deploy)
  
  # ---- site latitude & region
  site_lat <- deploy %>%
    transmute(locationID = .data[[location_col]],
              lat = suppressWarnings(as.numeric(.data[[lat_col]]))) %>%
    filter(is.finite(lat)) %>%
    group_by(locationID) %>%
    summarise(site_lat = mean(lat), .groups = "drop") %>%
    mutate(region = if_else(site_lat >= cutoff_lat, "North (≥52°N)", "South (<52°N)"))
  
  # ---- parsed deployments + region
  dep <- deploy %>%
    transmute(
      deploymentID = .data[[depl_id_deploy]],
      locationID   = .data[[location_col]],
      start = parse_ts(.data[[start_col]]),
      end   = parse_ts(.data[[end_col]])
    ) %>%
    filter(!is.na(start), !is.na(end), end > start) %>%
    inner_join(site_lat, by = "locationID")
  
  # ---- anchor slices for this dataset
  slices <- anchors %>%
    filter(dataset == ds) %>%
    arrange(anchor_start) %>%
    transmute(slice_id = row_number(),
              slice_start = as.Date(anchor_start),
              slice_end   = as.Date(anchor_end))
  if (nrow(slices) == 0) stop("No anchors for dataset '", ds, "'.")
  
  # helper: longest run of zeros in a 0/1 vector
  longest_zeros <- function(x) {
    r <- rle(x == 0L)
    if (!any(r$values)) 0L else max(r$lengths[r$values])
  }
  
  slice_stats <- function(s_start, s_end) {
    # clip to slice
    dep_clip <- dep %>%
      mutate(st = pmax(as_date(start), s_start),
             en = pmin(as_date(end),   s_end)) %>%
      filter(en >= st)
    
    # sites present in slice (any overlap)
    sites_in <- dep_clip %>% distinct(region, locationID)
    
    # region-level trap-days (site-days across all deployments)
    effort_tbl <- dep_clip %>%
      mutate(date = map2(st, en, ~ seq(.x, .y, by = "day"))) %>%
      unnest(date) %>%
      filter(!(month(date) == 2 & mday(date) == 29)) %>%
      count(region, date, name = "trap_days") %>%
      group_by(region) %>%
      summarise(trap_days_total = sum(trap_days), .groups = "drop")
    
    # location-wise active days within each site's own active span
    active_days <- dep_clip %>%
      mutate(date = purrr::map2(st, en, ~ seq(.x, .y, by = "day"))) %>%
      select(region, locationID, date) %>%
      unnest(date) %>%
      filter(!(lubridate::month(date) == 2 & lubridate::mday(date) == 29)) %>%
      distinct(region, locationID, date)
    
    # site-specific bounds (first/last active day inside slice)
    bounds <- active_days %>%
      group_by(region, locationID) %>%
      summarise(t0 = min(date), t1 = max(date), .groups = "drop")
    
    # fill only between t0..t1 for each site (no padding to whole slice)
    loc_daily <- bounds %>%
      mutate(date = purrr::map2(t0, t1, ~ seq(.x, .y, by = "day"))) %>%
      unnest(date) %>%
      left_join(active_days %>% mutate(active = 1L),
                by = c("region","locationID","date")) %>%
      mutate(active = tidyr::replace_na(active, 0L)) %>%
      select(region, locationID, date, active)
    
    
    # longest zero run per site, then worst per region
    gap_by_site <- loc_daily %>%
      group_by(region, locationID) %>%
      summarise(max_gap_site = longest_zeros(active), .groups = "drop")
    
    worst <- gap_by_site %>%
      group_by(region) %>%
      slice_max(max_gap_site, n = 1, with_ties = FALSE) %>%
      rename(worst_site = locationID,
             max_gap_days_location = max_gap_site) %>%
      ungroup()
    
    sites_n <- sites_in %>% count(region, name = "n_camera_sites")
    
    sites_n %>%
      full_join(effort_tbl, by = "region") %>%
      full_join(worst,      by = "region") %>%
      mutate(slice_start = s_start, slice_end = s_end) %>%
      relocate(slice_start, slice_end, region) %>%
      arrange(region)
  }
  
  purrr::pmap_dfr(slices, ~ slice_stats(..2, ..3)) %>%
    mutate(dataset = ds) %>%
    relocate(dataset, slice_start, slice_end, region)
}


mica_loc_stats <- split_location_stats("C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/MICA")  # folder or zip
mica_loc_stats

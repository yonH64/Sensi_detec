# ─────────────────────────────────────────────────────────
# Libraries 
# ─────────────────────────────────────────────────────────
library(tidyverse)
library(lubridate)
library(viridis)
library(glue)
library(GGally)
library(ggplot2)
library(stringr)
library(geosphere)
library(mapview)
library(sf)
library(leafpop)


############################################################
#                           HELPERS                        #
############################################################
source("helpers.R")


# Protocol windows (SNAP_EU_CORE, SNAP_EU_BUFFER, EOW_EARLY, EOW_LATE) are now included
# by make_window_template() via protocol_windows() in helpers.R.
window_grid <- make_window_template(step_doy = 7, lengths_d = seq(15, 120, by = 7))

# --------------------------------------------------------------------
# compute_full_effort()
# --------------------------------------------------------------------
compute_full_effort <- function(deploy_tbl, a_start) {
  a0 <- as.Date(a_start)
  
  dep_days <- deploy_tbl %>%
    transmute(
      camera_id = camera_id,
      start_day = as.integer(as.Date(start) - a0) + 1L,
      end_day   = as.integer(as.Date(end)   - a0) + 1L
    ) %>%
    filter(!is.na(start_day), !is.na(end_day), end_day >= start_day) %>%
    mutate(
      start_day   = pmax(start_day, 1L),
      end_day     = pmin(end_day,   365L),
      effort_days = pmax(0L, end_day - start_day + 1L)
    ) %>%
    group_by(camera_id) %>%
    summarise(effort_days = sum(effort_days), .groups = "drop") %>%
    filter(effort_days > 0)
  
  tibble::tibble(
    trap_days_full = sum(dep_days$effort_days),
    n_sites_full   = nrow(dep_days)
  )
}
# --------------------------------------------------------------------

############################################################
#                           SLICES                         #
############################################################

# split_spatial_by_distance()
# --------------------------------------------------------------------
# Aim
#   Split a dataset into sub-datasets when STATIONS form spatially
#   separated clusters beyond a given threshold (km). Clustering is done on
#   unique station coordinates (averaged per station). Deployments inherit
#   their station’s cluster; observations inherit the cluster via deployment.
#
# Arguments
#   dataset_name : Name of the dataset folder/archive under `root`.
#   root         : Root directory that contains all datasets (folders or zips).
#   threshold_km : Maximum within-cluster station-to-station distance (km).
#                  Stations connected by a chain of ≤ threshold edges end up
#                  in the same connected component (cluster). Default: 20.
#   min_stations : Minimum number of STATIONS required to keep/write a cluster.
#                  Clusters with fewer stations are discarded. Default: 10.
#   export_dir   : Where to create sub-dataset folders. Defaults to `root`.
#   include_zips : Whether to resolve archives (zip/7z/tar/gz/bz2/rar). TRUE.
#   overwrite    : Overwrite existing CSVs if they exist. Default TRUE.
#
# Output
#   Character vector of created subfolder full paths (in write order).
#   If the dataset forms a single cluster or all clusters are < min_stations,
#   nothing is written and an empty character vector is returned (with message).
#
# Notes
#   - Uses helpers: resolve_dataset_path(), read_csv_src(), set_camtrap_cols(),
#     parse_num_safe(). Requires: igraph, geosphere, dplyr, readr.
#   - Deployments without valid station coordinates are dropped (warned); any
#     observations attached to those deployments are also excluded.
# --------------------------------------------------------------------
split_spatial_by_distance <- function(dataset_name,
                                      threshold_km = 20,
                                      min_stations = 10,
                                      export_dir   = NULL,
                                      include_zips = TRUE,
                                      overwrite    = TRUE) {
  suppressPackageStartupMessages({
    library(tidyverse); library(igraph); library(geosphere)
  })
  
  # ---- resolve and read -----------------------------------------------------
  root         <- ("C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets")
  archive_path <- resolve_dataset_path(dataset_name, root, include_zips = include_zips)
  ds_base      <- tools::file_path_sans_ext(basename(archive_path))
  out_root     <- if (is.null(export_dir)) root else export_dir
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  
  deploy <- read_csv_src(archive_path, "deployments\\.csv")
  obs    <- read_csv_src(archive_path, "observations\\.csv")
  cols <- set_camtrap_cols(obs, deploy)
  list2env(cols, envir = environment())
  
  # ---- station coordinates (averaged per station) ---------------------------
  parse_num <- function(x) parse_num_safe(x)
  
  station_xy <- deploy %>%
    transmute(station = .data[[loc_col]],
              lat     = parse_num(.data[[lat_col]]),
              lon     = parse_num(.data[[lon_col]])) %>%
    filter(!is.na(station)) %>%
    group_by(station) %>%
    summarise(lat = mean(lat, na.rm = TRUE),
              lon = mean(lon, na.rm = TRUE),
              .groups = "drop") %>%
    filter(is.finite(lat), is.finite(lon))
  
  n_st <- nrow(station_xy)
  if (n_st == 0) stop("No valid station coordinates (lat/lon) in deployments for: ", ds_base)
  
  # ---- graph of stations with edges ≤ threshold_km --------------------------
  dmat_km <- geosphere::distm(station_xy[, c("lon","lat")], fun = geosphere::distHaversine) / 1000
  adj     <- (dmat_km <= threshold_km) & is.finite(dmat_km)
  diag(adj) <- FALSE
  
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  comp <- igraph::components(g)
  cluster_id <- paste0("C", comp$membership)
  
  # station -> cluster
  station_cluster <- station_xy %>% mutate(cluster = cluster_id)
  
  # ---- deployments → clusters ----------------------------------------------
  dep_with_cluster <- deploy %>%
    mutate(station = .data[[loc_col]]) %>%
    left_join(station_cluster, by = "station") %>%
    mutate(has_cluster = !is.na(cluster))
  
  n_drop_dep <- sum(!dep_with_cluster$has_cluster, na.rm = TRUE)
  if (n_drop_dep > 0) {
    warning("Dropping ", n_drop_dep, " deployment(s) lacking valid station coordinates; ",
            "their observations will be excluded from outputs.")
  }
  dep_with_cluster <- dep_with_cluster %>% filter(has_cluster) %>% select(-has_cluster)
  
  # ---- observations → clusters via deployment -------------------------------
  dep_to_cluster <- dep_with_cluster %>%
    transmute(deploymentID = .data[[depl_id_deploy]], cluster)
  
  obs_with_cluster <- obs %>%
    transmute(deploymentID = .data[[depl_id_obs]],
              across(everything())) %>%
    left_join(dep_to_cluster, by = "deploymentID") %>%
    filter(!is.na(cluster))
  
  # ---- cluster sizes (stations + deployments) -------------------------------
  stations_per_cluster <- station_cluster %>% count(cluster, name = "n_stations")
  deploys_per_cluster  <- dep_with_cluster %>% count(cluster, name = "n_deployments")
  cl_info <- stations_per_cluster %>%
    left_join(deploys_per_cluster, by = "cluster") %>%
    arrange(desc(n_stations), desc(n_deployments))
  
  n_clusters <- nrow(cl_info)
  if (n_clusters <= 1) {
    message("No split needed: all stations form a single spatial cluster at ", threshold_km, " km.")
    return(character(0))
  }
  
  # ---- apply min_stations threshold -----------------------------------------
  keep_cl <- cl_info %>% filter(n_stations >= min_stations)
  drop_cl <- cl_info %>% filter(n_stations <  min_stations)
  
  if (nrow(keep_cl) == 0) {
    message("All clusters have fewer than ", min_stations, " stations; nothing written.")
    return(character(0))
  }
  if (nrow(drop_cl) > 0) {
    message("Discarding ", nrow(drop_cl), " cluster(s) with < ", min_stations, " stations: ",
            paste0(drop_cl$cluster, " (", drop_cl$n_stations, " stn) ", collapse = ", "))
  }
  
  # ---- write kept clusters ---------------------------------------------------
  created <- character(0)
  for (k in seq_len(nrow(keep_cl))) {
    cl_lab <- keep_cl$cluster[k]
    out_dir <- file.path(out_root, paste0(ds_base, "_", cl_lab))
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    dep_file <- file.path(out_dir, "deployments.csv")
    obs_file <- file.path(out_dir, "observations.csv")
    
    dep_k <- dep_with_cluster %>% filter(cluster == cl_lab) %>% select(-cluster, -station, -lat, -lon)
    obs_k <- obs_with_cluster %>% filter(cluster == cl_lab) %>% select(-cluster)
    
    if (overwrite || !file.exists(dep_file)) readr::write_csv(dep_k, dep_file)
    if (overwrite || !file.exists(obs_file)) readr::write_csv(obs_k, obs_file)
    
    message("✓ Wrote ", basename(out_dir), "/{deployments.csv,observations.csv}  (",
            keep_cl$n_stations[k], " stations; ",
            nrow(dep_k), " deployments; ",
            nrow(obs_k), " obs)")
    created <- c(created, normalizePath(out_dir, mustWork = FALSE))
  }
  
  # ---- remove original dataset folder if split produced outputs -------------
  if (length(created) > 0) {
    # only delete if original is a directory (not an archive)
    if (dir.exists(archive_path)) {
      # safety: only delete when it's directly under `root`
      ap_norm   <- normalizePath(archive_path, winslash = "/", mustWork = FALSE)
      root_norm <- normalizePath(root,         winslash = "/", mustWork = FALSE)
      if (dirname(ap_norm) == root_norm) {
        unlink(archive_path, recursive = TRUE, force = TRUE)
        message("🗑️ Deleted original dataset folder: ", archive_path)
      } else {
        message("Skipping deletion: original path is not directly under root: ", archive_path)
      }
    } else {
      message("Original dataset is an archive (not a folder); nothing to delete: ", archive_path)
    }
  }
  
  invisible(created)
}



# Default: keep clusters that have at least 10 stations
split_spatial_by_distance("NL-MICA",
                          export_dir   = "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets",
                          threshold_km = 30,
                          min_stations = 10)



# find_anchors() — moved to helpers.R
# See helpers.R for full documentation and implementation.



ds_paths <- list.dirs(
  "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets",
  recursive = FALSE, full.names = TRUE
)

anchors <- furrr::future_map_dfr(
  ds_paths,
  find_anchors,
  .options = furrr::furrr_options(seed = TRUE)  # default scheduling=1
)

############################################################
#                           PARSER                         #
############################################################

# --------------------------------------------------------------------
# build_window_metrics_fast1()
# --------------------------------------------------------------------
# Aim
#   For one 12-month “anchor” period in a dataset:
#     1) Collapse detections into independent events (per camera × species).
#     2) Intersect deployments with each sampling window from a shared grid.
#     3) For each window, build camera × occasion effort/detection matrices.
#     4) Apply minimum-data thresholds per species × window.
#     5) Compute detection-derived metrics (occupancy p̂, time-to-event P, rate)
#        and their uncertainties, and deltas to the full-year (FULL) benchmark.
#     6) Compute community metrics: observed richness (sr_obs) and incidence-
#        rarefied richness (sr_raref) to a common sampling intensity.
#
# Notes (interpretation)
#   - “FULL” is treated as a benchmark (not “truth”) for year-like reporting.
#   - Occupancy is computed for completeness; long windows can violate closure.
#
# Sampling-unit definition for rarefaction
#   - One “unit” = one camera × occasion cell with >0 effort days.
#   - Q_i = number of units where species i was detected (incidence frequency).
#   - sr_raref is the expected richness after rarefying to units_ref = min(T)
#     across windows (T = number of sampled units in the window).
#
# Arguments (key)
#   window_grid        : Tibble with window_id, start_virtual, end_virtual,
#                        start_date, end_date, length_d, start_doy (FULL optional).
#   deploy_tbl         : Deployments cropped to the anchor; columns:
#                        camera_id, start, end, start_virtual, end_virtual.
#   obs_tbl            : Observations cropped to the anchor; columns:
#                        species, timestamp (POSIXct), camera_id.
#   a_start, a_end     : POSIXct anchor boundaries (wrap logic uses a_start).
#   independence_mins  : Minimum time gap (minutes) to define independent events.
#   occ_grain_days     : Occasion width (days) for camera × occasion matrices.
#   min_events         : Minimum total events per species × window.
#   min_occasions_pos  : Minimum positive occasions per species × window.
#   min_sites_pos      : Minimum positive sites (cameras) per species × window.
#   trap_array_km      : Study-level camera array diameter (km) to attach.
#   latitude_deg       : Study centroid latitude to attach.
#   drop_leap_day      : If TRUE, drop Feb 29 detections before processing.
#
# Output
#   A list with:
#     $window_species   — per species × window metrics (and deltas to FULL).
#     $window_camera    — camera × window effort table (trap-days).
#     $dropped_species  — species × window combos failing thresholds (flags).
#     $window_richness  — per window richness metrics + deltas to FULL:
#                         sr_obs, sr_raref, units_surveyed, units_ref, etc.
#     $y_mats           — optional debug storage of y matrices for selected cases.
# --------------------------------------------------------------------
build_window_metrics_fast1 <- function(
    window_grid,
    deploy_tbl,
    obs_tbl,
    a_start,
    a_end,
    independence_mins = 30,
    occ_grain_days = 1,
    min_events = 20,          # Min events for reliable TTE rate estimation (~20-30 recommended)
    min_occasions_pos = 3,    # Min temporal spread (reduced from 5 for TTE)
    min_sites_pos = 10,       # Min cameras with detections for descriptive TTE metrics
    trap_array_km = NA_real_,
    latitude_deg  = NA_real_,
    drop_leap_day = TRUE,
    debug_cases = NULL   # tibble with columns: window_id, species; optional
) {
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(lubridate)
    library(tibble)
  })
  
  y_mats <- list()
  
  # ---- 0) Convert timestamps to integer “day within anchor” (1..365) ----
  a0 <- as.Date(a_start)
  
  dep_dt <- deploy_tbl %>%
    transmute(
      camera_id = camera_id,
      start_day = as.integer(as.Date(start) - a0) + 1L,
      end_day   = as.integer(as.Date(end)   - a0) + 1L
    ) %>%
    filter(!is.na(start_day), !is.na(end_day), end_day >= start_day)
  
  if (!nrow(dep_dt)) {
    return(list(
      window_species   = tibble(),
      window_camera    = tibble(),
      dropped_species  = tibble(),
      window_richness  = tibble(),
      y_mats           = list()
    ))
  }
  
  dep_ivl_all <- as.data.table(dep_dt)[, .(camera_id, start = start_day, end = end_day)]
  setkey(dep_ivl_all, start, end)
  
  # ---- 1) Collapse detections to independent events (camera × species) ----
  events <- obs_tbl %>%
    mutate(timestamp = with_tz(timestamp, tzone = "UTC")) %>%
    { if (drop_leap_day) filter(., !(month(timestamp) == 2 & mday(timestamp) == 29)) else . } %>%
    arrange(camera_id, species, timestamp) %>%
    group_by(camera_id, species) %>%
    mutate(
      gap = timestamp - lag(timestamp),
      new_event = is.na(gap) | gap > minutes(independence_mins)
    ) %>%
    filter(new_event) %>%
    ungroup() %>%
    mutate(day_id = as.integer(as.Date(timestamp) - a0) + 1L) %>%
    filter(!is.na(day_id), day_id >= 1L, day_id <= 365L)
  
  # If there are no events in the year, we still compute effort + richness (=0).
  # Species-level tables will be empty.
  
  # ---- 2) Build window segments (split wrapped windows into two segments) ----
  wg <- as.data.table(window_grid)
  win_start_map <- setNames(as.integer(wg$start_virtual), wg$window_id)
  
  win_seg <- rbindlist(list(
    wg[, .(window_id,
           w_start = as.integer(start_virtual),
           w_end   = as.integer(pmin(end_virtual, 365L)))],
    wg[end_virtual > 365L,
       .(window_id,
         w_start = 1L,
         w_end   = as.integer(end_virtual - 365L))]
  ), use.names = TRUE)
  
  setkey(win_seg, w_start, w_end)
  
  # ---- 3) Deployment × window overlap → effort per camera × window ----
  dep_ivl <- as.data.table(dep_dt)[, .(camera_id, start = start_day, end = end_day)]
  setkey(dep_ivl, start, end)
  
  dep_ov <- foverlaps(
    dep_ivl, win_seg,
    by.x = c("start","end"),
    by.y = c("w_start","w_end"),
    type = "any",
    nomatch = 0L
  )
  
  pick_col <- function(dt, ...) {
    cand <- c(...)
    hit <- cand[cand %in% names(dt)][1]
    if (is.na(hit)) stop("Missing expected column(s): ", paste(cand, collapse=", "),
                         "\nHave: ", paste(names(dt), collapse=", "))
    hit
  }
  
  id_col <- pick_col(dep_ov, "i.window_id", "window_id")
  ws_col <- pick_col(dep_ov, "i.w_start",   "w_start")
  we_col <- pick_col(dep_ov, "i.w_end",     "w_end")
  
  dep_ov[, effort_days := pmin(end, get(we_col)) - pmax(start, get(ws_col)) + 1L]
  
  window_camera <- dep_ov[effort_days > 0L, .(
    window_id = get(id_col),
    camera_id,
    effort_days
  )]
  
  # Wrapped windows can duplicate cameras; aggregate once.
  window_camera <- window_camera[, .(effort_days = sum(effort_days)), by = .(window_id, camera_id)]
  
  # ---- 4) Assign events to windows (vectorised) ----
  if (nrow(events)) {
    ev_ivl <- as.data.table(events)[, .(camera_id, species, timestamp, start = day_id, end = day_id)]
    setkey(ev_ivl, start, end)
    
    ev_ov <- foverlaps(
      ev_ivl, win_seg,
      by.x = c("start","end"),
      by.y = c("w_start","w_end"),
      type = "within",
      nomatch = 0L
    )
    
    window_events_raw <- ev_ov[, .(window_id, camera_id, species, timestamp, day_id = start)]
  } else {
    window_events_raw <- data.table(
      window_id = character(),
      camera_id = character(),
      species   = character(),
      timestamp = as.POSIXct(character()),
      day_id    = integer()
    )
  }
  
  # ---- 5) Loop windows; within each window loop species (fits can’t be vectorised) ----
  results_list <- vector("list", length = nrow(window_grid))
  names(results_list) <- window_grid$window_id
  
  fail_log <- vector("list", 0)
  
  setDT(window_camera);     setkey(window_camera, window_id)
  setDT(window_events_raw); setkey(window_events_raw, window_id)
  
  # Richness collectors (one row per window)
  rich_list <- vector("list", length = nrow(window_grid))
  names(rich_list) <- window_grid$window_id
  
  # Rarefaction collectors (per window)
  Q_list <- vector("list", length = nrow(window_grid))  # named integer vector (species -> Q_i)
  T_list <- integer(nrow(window_grid))                  # T = number of surveyed units
  
  for (i in seq_len(nrow(window_grid))) {
    win_id <- window_grid$window_id[i]
    L_days <- as.integer(window_grid$length_d[i])
    n_occ  <- as.integer(ceiling(L_days / occ_grain_days))
    
    cam_eff <- window_camera[J(win_id)]         # camera × window effort (trap-days)
    ev_sub  <- window_events_raw[J(win_id)]     # events pre-filtered to window segment(s)
    
    if (!nrow(cam_eff) || n_occ <= 0L) {
      # No sampling effort (or invalid window), richness is 0 and no species metrics
      rich_list[[i]] <- tibble(
        window_id        = win_id,
        window_len       = L_days,
        window_start     = window_grid$start_date[i],
        trap_days_window = 0,
        n_sites          = 0L,
        units_surveyed   = 0L,
        sr_obs           = 0L
      )
      Q_list[[i]] <- integer(0)
      T_list[i]   <- 0L
      results_list[[i]] <- tibble()
      next
    }
    
    wstart <- win_start_map[[win_id]]
    
    # ---- 5a) Build occasion segments for this window (wrap-safe) ----
    occ_tbl <- data.table(occ_id = seq_len(n_occ))
    occ_tbl[, k_start := (occ_id - 1L) * occ_grain_days]
    occ_tbl[, k_end   := pmin(occ_id * occ_grain_days - 1L, L_days - 1L)]
    occ_tbl[, o_start := ((wstart - 1L + k_start) %% 365L) + 1L]
    occ_tbl[, o_end   := ((wstart - 1L + k_end)   %% 365L) + 1L]
    
    occ_seg <- rbindlist(list(
      occ_tbl[o_start <= o_end, .(occ_id, o_start, o_end)],
      occ_tbl[o_start >  o_end, .(occ_id, o_start, o_end = 365L)],
      occ_tbl[o_start >  o_end, .(occ_id, o_start = 1L, o_end)]
    ), use.names = TRUE)
    setkey(occ_seg, o_start, o_end)
    
    # ---- 5b) Deployment × occasion overlap → effort per camera × occasion ----
    # Speed: restrict deployments to cameras that actually have window-level effort.
    dep_ivl_w <- dep_ivl_all[camera_id %in% cam_eff$camera_id]
    setkey(dep_ivl_w, start, end)
    
    ov_occ <- foverlaps(
      dep_ivl_w, occ_seg,
      by.x = c("start","end"),
      by.y = c("o_start","o_end"),
      type = "any",
      nomatch = 0L
    )
    
    if (!nrow(ov_occ)) {
      rich_list[[i]] <- tibble(
        window_id        = win_id,
        window_len       = L_days,
        window_start     = window_grid$start_date[i],
        trap_days_window = as.numeric(sum(cam_eff$effort_days)),
        n_sites          = as.integer(uniqueN(cam_eff$camera_id)),
        units_surveyed   = 0L,
        sr_obs           = 0L
      )
      Q_list[[i]] <- integer(0)
      T_list[i]   <- 0L
      results_list[[i]] <- tibble()
      next
    }
    
    ov_occ[, eff_days := pmin(end, o_end) - pmax(start, o_start) + 1L]
    cam_eff_occ <- ov_occ[eff_days > 0L, .(effort_days = sum(eff_days)), by = .(camera_id, occ_id)]
    
    if (!nrow(cam_eff_occ)) {
      rich_list[[i]] <- tibble(
        window_id        = win_id,
        window_len       = L_days,
        window_start     = window_grid$start_date[i],
        trap_days_window = as.numeric(sum(cam_eff$effort_days)),
        n_sites          = as.integer(uniqueN(cam_eff$camera_id)),
        units_surveyed   = 0L,
        sr_obs           = 0L
      )
      Q_list[[i]] <- integer(0)
      T_list[i]   <- 0L
      results_list[[i]] <- tibble()
      next
    }
    
    # ---- 5c) Build effort matrix (camera × occasion) ----
    cam_ids <- sort(unique(cam_eff_occ$camera_id))
    
    eff_mat <- matrix(
      0L,
      nrow = length(cam_ids),
      ncol = n_occ,
      dimnames = list(cam_ids, paste0("oc", seq_len(n_occ)))
    )
    
    ridx <- match(cam_eff_occ$camera_id, cam_ids)
    cidx <- cam_eff_occ$occ_id
    ok   <- !is.na(ridx) & !is.na(cidx)
    if (any(ok)) eff_mat[cbind(ridx[ok], cidx[ok])] <- cam_eff_occ$effort_days[ok]
    
    # Totals for this window
    trap_days <- sum(eff_mat, na.rm = TRUE)
    n_sites   <- length(cam_ids)
    
    # Rarefaction “sample size” = number of surveyed camera × occasion cells
    T_units <- sum(eff_mat > 0L, na.rm = TRUE)
    T_list[i] <- as.integer(T_units)
    
    # ---- 5d) Prepare events for this window (can be empty) ----
    if (nrow(ev_sub)) {
      # Convert absolute day_id into “offset within the window”
      ev_sub[, offset := ((day_id - wstart) %% 365L)]
      ev_sub <- ev_sub[offset < L_days]
      
      # Assign occasions
      ev_sub[, occ_id := (offset %/% occ_grain_days) + 1L]
      ev_sub <- ev_sub[occ_id >= 1L & occ_id <= n_occ]
      
      # Keep only events from cameras that were actually sampled in this window
      ev_sub <- ev_sub[camera_id %in% cam_ids]
      
      # Keep only events in surveyed camera × occasion cells (effort > 0)
      ev_sub[, unit_id := paste(camera_id, occ_id, sep = "__")]
      sampled_units <- paste(cam_eff_occ$camera_id, cam_eff_occ$occ_id, sep = "__")
      ev_sub <- ev_sub[unit_id %in% sampled_units]
    }
    
    # ---- 5e) Species richness (observed + incidence frequencies for rarefaction) ----
    sr_obs <- if (nrow(ev_sub)) uniqueN(ev_sub$species) else 0L
    
    if (nrow(ev_sub) && T_units > 0L) {
      # Q_i = number of sampled units where species i was detected
      Q_dt <- ev_sub[, .(Q = uniqueN(unit_id)), by = species]
      Q_list[[i]] <- setNames(as.integer(Q_dt$Q), Q_dt$species)
    } else {
      Q_list[[i]] <- integer(0)
    }
    
    # Store richness row now (rarefied richness is added after looping all windows)
    rich_list[[i]] <- tibble(
      window_id        = win_id,
      window_len       = L_days,
      window_start     = window_grid$start_date[i],
      trap_days_window = as.numeric(trap_days),
      n_sites          = as.integer(n_sites),
      units_surveyed   = as.integer(T_units),
      sr_obs           = as.integer(sr_obs)
    )
    
    # If there are no detections, skip species-level calculations but keep richness/effort
    if (!nrow(ev_sub)) {
      results_list[[i]] <- tibble()
      next
    }
    
    # ---- 5f) Species-level metrics (thresholded) ----
    sp_vec <- unique(ev_sub$species)
    sp_out <- vector("list", length(sp_vec))
    names(sp_out) <- sp_vec
    
    for (sp in sp_vec) {
      ev_sp <- ev_sub[species == sp]
      n_events_total <- nrow(ev_sp)
      n_sites_pos    <- uniqueN(ev_sp$camera_id)
      n_occ_pos      <- uniqueN(ev_sp$occ_id)
      
      if (n_events_total >= min_events &&
          n_occ_pos      >= min_occasions_pos &&
          n_sites_pos    >= min_sites_pos) {
        
        # camera × occasion detection matrix (1 = detected, 0 = not detected, NA = not surveyed)
        mat <- matrix(
          0L,
          nrow = length(cam_ids),
          ncol = n_occ,
          dimnames = list(cam_ids, paste0("oc", seq_len(n_occ)))
        )
        
        ridx <- match(ev_sp$camera_id, cam_ids)
        cidx <- ev_sp$occ_id
        ok   <- !is.na(ridx) & !is.na(cidx)
        if (any(ok)) mat[cbind(ridx[ok], cidx[ok])] <- 1L
        
        # Cells with zero effort are unsurveyed → NA (important for occupancy framework)
        mat[eff_mat == 0L] <- NA_integer_
        
        if (!is.null(debug_cases) &&
            nrow(debug_cases) &&
            any(debug_cases$window_id == win_id & debug_cases$species == sp)) {
          key <- paste(win_id, sp, sep = "__")
          y_mats[[key]] <- mat
        }
        
        # ---- Camera-level detection summaries ----
        # For each camera: detection (0/1) and event count
        cam_detections <- rowSums(mat == 1L, na.rm = TRUE) > 0  # logical: any detection?
        cam_event_counts <- sapply(cam_ids, function(cid) {
          sum(ev_sp$camera_id == cid, na.rm = TRUE)
        })
        
        # Get effort per camera for this window (sum across all occasions)
        cam_effort <- rowSums(eff_mat, na.rm = TRUE)
        
        # ---- Spatial detection coverage ----
        spatial_cov <- n_sites_pos / n_sites
        spatial_cov_se <- sqrt(spatial_cov * (1 - spatial_cov) / n_sites)
        
        # ---- Time-to-event (TTE) proxy: first detection time per camera ----
        first_time <- apply(mat, 1, function(x) {
          w <- which(x == 1L)
          if (length(w)) (w[1] - 1) * occ_grain_days + 0.001 else NA_real_
        })
        
        censored <- is.na(first_time)
        n_events_cam <- sum(!censored)
        exposure <- sum(first_time[!censored]) + sum(censored) * L_days
        
        lambda <- if (exposure > 0) n_events_cam / exposure else 0
        se_lambda <- if (exposure > 0 && n_events_cam > 0) sqrt(n_events_cam) / exposure else NA_real_
        
        # ---- Event rate (events per trap-day) ----
        # Camera-level rates (empirical variance approach)
        cam_rates <- cam_event_counts / cam_effort
        cam_rates[!is.finite(cam_rates)] <- NA_real_
        
        rate <- if (trap_days > 0) n_events_total / trap_days else NA_real_
        # SE from camera-level variability (accounts for spatial heterogeneity)
        rate_se <- if (n_sites >= 2 && sum(!is.na(cam_rates)) >= 2) {
          sd(cam_rates, na.rm = TRUE) / sqrt(sum(!is.na(cam_rates)))
        } else {
          NA_real_
        }
        
        # log-rate for deltas (small epsilon avoids log(0))
        log_rate <- log(rate + 1e-6)
        log_rate_se <- if (!is.na(rate_se) && !is.na(rate)) rate_se / (rate + 1e-6) else NA_real_
        
        # ---- Camera-level detail for matched-camera analysis ----
        # cam_det_ids: cameras that detected this species in this window
        # cam_rate_tbl: per-camera event rate (events/effort) for paired comparisons
        cam_det_ids_sp <- cam_ids[cam_detections]
        cam_rate_tbl_sp <- tibble(
          camera_id    = cam_ids,
          cam_events   = as.integer(cam_event_counts),
          cam_effort   = as.numeric(cam_effort),
          cam_rate     = cam_rates
        ) |> filter(cam_effort > 0)

        sp_out[[sp]] <- tibble(
          window_id        = win_id,
          window_len       = L_days,
          window_start     = window_grid$start_date[i],
          species          = sp,
          n_events_total   = as.integer(n_events_total),
          n_occasions_pos  = as.integer(n_occ_pos),
          n_sites_pos      = as.integer(n_sites_pos),
          trap_days_window = as.numeric(trap_days),
          n_sites          = as.integer(n_sites),
          trap_array       = trap_array_km,
          latitude         = latitude_deg,
          spatial_cov      = spatial_cov,
          spatial_cov_se   = spatial_cov_se,
          lambda           = lambda,
          lambda_se        = se_lambda,
          rate             = rate,
          rate_se          = rate_se,
          log_rate         = log_rate,
          log_rate_se      = log_rate_se,
          cam_det_ids      = list(cam_det_ids_sp),
          cam_rate_tbl     = list(cam_rate_tbl_sp)
        )
        
      } else {
        fail_log[[length(fail_log) + 1]] <- tibble(
          window_id                = win_id,
          window_len               = L_days,
          window_start             = window_grid$start_date[i],
          species                  = sp,
          n_events_total           = as.integer(n_events_total),
          n_occasions_pos          = as.integer(n_occ_pos),
          n_sites_pos              = as.integer(n_sites_pos),
          failed_min_events        = n_events_total  < min_events,
          failed_min_occasions_pos = n_occ_pos      < min_occasions_pos,
          failed_min_sites_pos     = n_sites_pos    < min_sites_pos
        )
      }
    }
    
    results_list[[i]] <- bind_rows(sp_out)
  }
  
  window_species  <- bind_rows(results_list)
  dropped_species <- bind_rows(fail_log)
  window_richness <- bind_rows(rich_list)
  
  # ---- 6) Incidence-based rarefaction to common number of sampled units ----
  # E[S_m] = sum_i [1 - C(T - Q_i, m) / C(T, m)]
  sr_raref_incidence <- function(Q, T, m) {
    if (m <= 0 || T <= 0 || length(Q) == 0) return(0)
    Q <- pmin(pmax(Q, 0L), T)
    term <- exp(lchoose(T - Q, m) - lchoose(T, m))
    sum(1 - term)
  }
  
  # Var[E(S_m)] — analytical variance of the incidence rarefaction estimator
  # (Colwell et al. 2012, eq. 5):
  #   Var = sum_i p_i(1 - p_i)
  #       + 2 * sum_{i<j} [ C(T - Q_i - Q_j, m) / C(T, m)  -  p_i^c * p_j^c ]
  # where p_i^c = C(T - Q_i, m) / C(T, m)  (prob species i NOT detected in m units)
  #       p_i   = 1 - p_i^c                 (prob species i detected)
  # Returns NA if Q_i + Q_j > T for any pair (would give lchoose of negative = -Inf)
  sr_raref_var <- function(Q, T, m) {
    if (m <= 0 || T <= 0 || length(Q) == 0) return(0)
    Q  <- pmin(pmax(Q, 0L), T)
    S  <- length(Q)
    lc_T_m  <- lchoose(T, m)                           # log C(T, m)
    pc <- exp(lchoose(T - Q, m) - lc_T_m)              # p_i^c for each species
    p  <- 1 - pc                                        # p_i
    
    # Diagonal (independent) terms
    var_diag <- sum(p * pc)
    
    # Off-diagonal (covariance) terms — O(S^2); fast for typical S < 200
    if (S < 2L) return(var_diag)
    
    # Pairwise C(T - Q_i - Q_j, m) / C(T, m)
    # If T - Q_i - Q_j < m the species pair never co-occurs → term = 0
    Qi <- Q[rep(seq_len(S), times = S)]
    Qj <- Q[rep(seq_len(S), each  = S)]
    keep <- Qi < Qj                                     # upper triangle only
    Qi <- Qi[keep]; Qj <- Qj[keep]
    pc_i <- pc[rep(seq_len(S), times = S)][keep]
    pc_j <- pc[rep(seq_len(S), each  = S)][keep]
    
    denom_arg <- T - Qi - Qj
    joint <- ifelse(denom_arg >= m,
                    exp(lchoose(denom_arg, m) - lc_T_m),
                    0)
    var_offdiag <- 2 * sum(joint - pc_i * pc_j)
    
    max(0, var_diag + var_offdiag)   # clamp: numerical noise can give tiny negatives
  }
  
  units_ref <- suppressWarnings(min(T_list[T_list > 0L], na.rm = TRUE))
  if (!is.finite(units_ref)) units_ref <- NA_integer_
  
  sr_raref_vec     <- rep(NA_real_, nrow(window_grid))
  sr_raref_var_vec <- rep(NA_real_, nrow(window_grid))
  
  if (!is.na(units_ref) && units_ref > 0L) {
    for (j in seq_along(Q_list)) {
      Tj <- T_list[j]
      if (Tj >= units_ref && Tj > 0L) {
        sr_raref_vec[j]     <- sr_raref_incidence(Q_list[[j]], Tj, units_ref)
        sr_raref_var_vec[j] <- sr_raref_var(Q_list[[j]], Tj, units_ref)
      } else if (Tj > 0L) {
        sr_raref_vec[j]     <- NA_real_
        sr_raref_var_vec[j] <- NA_real_
      } else {
        sr_raref_vec[j]     <- 0
        sr_raref_var_vec[j] <- 0
      }
    }
  }
  
  window_richness <- window_richness %>%
    mutate(
      units_ref    = as.integer(units_ref),
      sr_raref     = as.numeric(sr_raref_vec),
      sr_raref_var = as.numeric(sr_raref_var_vec),  # sampling variance of sr_raref
      sr_raref_se  = sqrt(as.numeric(sr_raref_var_vec))  # SE = sqrt(variance)
    )
  
  # ---- 7) Richness deltas to FULL benchmark ----
  # ---- Richness deltas relative to FULL (12-month benchmark) ----
  if ("FULL" %in% window_richness$window_id) {
    
    sr_full_val <- window_richness %>%
      filter(window_id == "FULL") %>%
      slice(1) %>%
      pull(sr_obs)
    
    sr_raref_full_val <- window_richness %>%
      filter(window_id == "FULL") %>%
      slice(1) %>%
      pull(sr_raref)
    
    # Species sets for Jaccard and proportion detected
    spp_full <- if (nrow(window_species)) {
      window_species %>%
        filter(window_id == "FULL") %>%
        pull(species) %>%
        unique()
    } else {
      character(0)
    }
    
    # Per-window species sets (excluding FULL)
    spp_by_window <- if (nrow(window_species)) {
      window_species %>%
        filter(window_id != "FULL") %>%
        group_by(window_id) %>%
        summarise(spp_list = list(unique(species)), .groups = "drop")
    } else {
      tibble(window_id = character(0), spp_list = list())
    }
    
    window_richness <- window_richness %>%
      filter(window_id != "FULL") %>%
      left_join(spp_by_window, by = "window_id") %>%
      mutate(
        # Observed richness vs FULL
        sr_full  = as.integer(sr_full_val),
        d_sr     = as.integer(sr_obs - sr_full_val),
        d_sr_rel = if (sr_full_val > 0) {
          as.numeric(sr_obs - sr_full_val) / as.numeric(sr_full_val)
        } else {
          NA_real_
        },
        
        # Rarefied richness vs FULL (same units_ref)
        sr_raref_full  = as.numeric(sr_raref_full_val),
        d_sr_raref     = as.numeric(sr_raref - sr_raref_full_val),
        d_sr_raref_rel = if (!is.na(sr_raref_full_val) && sr_raref_full_val > 0) {
          (as.numeric(sr_raref) - as.numeric(sr_raref_full_val)) / as.numeric(sr_raref_full_val)
        } else {
          NA_real_
        },
        
        # MSE criterion for rarefied richness: bias² + sampling variance
        # Mirrors the species-level MSE (mse_lambda, mse_rate).
        # sr_raref_var is the analytical incidence-rarefaction variance (Colwell 2012).
        mse_sr_raref = (as.numeric(sr_raref) - as.numeric(sr_raref_full_val))^2 +
          ifelse(!is.na(sr_raref_var), sr_raref_var, NA_real_),
        
        # Proportion of FULL species detected: |A ∩ B| / |B|
        # How many of the FULL species were recovered in this window?
        prop_sr_full = purrr::map_dbl(spp_list, function(spp_win) {
          if (length(spp_full) == 0) return(NA_real_)
          length(intersect(spp_win, spp_full)) / length(spp_full)
        })
      ) %>%
      select(-spp_list)
    
  } else {
    window_richness <- window_richness %>%
      mutate(
        sr_full = NA_integer_, d_sr = NA_integer_, d_sr_rel = NA_real_,
        sr_raref_full = NA_real_, d_sr_raref = NA_real_, d_sr_raref_rel = NA_real_,
        mse_sr_raref  = NA_real_,
        prop_sr_full  = NA_real_
      )
  }
  
  
  # ---- 8) Species deltas to FULL benchmark and MSE criterion ----
  if (nrow(window_species) && "FULL" %in% window_species$window_id) {
    truth_tbl <- window_species %>%
      filter(window_id == "FULL") %>%
      select(
        species,
        spatial_cov_full    = spatial_cov,
        spatial_cov_se_full = spatial_cov_se,
        lambda_full         = lambda,
        lambda_se_full      = lambda_se,
        rate_full           = rate,
        rate_se_full        = rate_se,
        cam_det_ids_full    = cam_det_ids,
        cam_rate_tbl_full   = cam_rate_tbl
      )
    
    window_species <- window_species %>%
      left_join(truth_tbl, by = "species") %>%
      mutate(
        # Deltas (bias) — lambda and rate are window-length-independent
        d_lambda      = lambda   - lambda_full,
        d_rate        = rate     - rate_full,
        
        # MSE criterion: (bias)² + Var
        mse_lambda       = (lambda - lambda_full)^2 + ifelse(!is.na(lambda_se), lambda_se^2, NA_real_),
        mse_rate         = (rate - rate_full)^2 + ifelse(!is.na(rate_se), rate_se^2, NA_real_)
      )
    
    # ---- 8a) Matched-camera metrics ----
    # For each non-FULL window × species, restrict to cameras active in BOTH
    # the sub-window and FULL, then recompute spatial_cov, rate, and lambda
    # on the shared camera set. This removes the denominator confound from
    # spatial_cov and enables camera-paired deviations for all metrics.
    
    compute_matched <- function(cam_rate_tbl_win, cam_det_ids_win,
                                cam_rate_tbl_full, cam_det_ids_full,
                                win_len) {
      # Fallback for missing data
      na_row <- tibble(
        n_matched            = NA_integer_,
        spatial_cov_m        = NA_real_,
        spatial_cov_m_se     = NA_real_,
        matched_rate         = NA_real_,
        rate_m               = NA_real_,
        rate_m_se            = NA_real_,
        log_rate_m           = NA_real_,
        log_rate_m_se        = NA_real_,
        spatial_cov_full_m   = NA_real_,
        spatial_cov_full_m_se = NA_real_,
        matched_rate_full    = NA_real_,
        rate_full_m          = NA_real_,
        rate_full_m_se       = NA_real_,
        log_rate_full_m      = NA_real_,
        log_rate_full_m_se   = NA_real_
      )
      if (is.null(cam_rate_tbl_win) || is.null(cam_rate_tbl_full)) return(na_row)
      
      win_tbl  <- cam_rate_tbl_win
      full_tbl <- cam_rate_tbl_full
      if (!nrow(win_tbl) || !nrow(full_tbl)) return(na_row)
      
      shared_cams <- intersect(win_tbl$camera_id, full_tbl$camera_id)
      n_matched   <- length(shared_cams)
      if (n_matched < 2L) return(na_row)
      
      win_m  <- win_tbl  |> filter(camera_id %in% shared_cams)
      full_m <- full_tbl |> filter(camera_id %in% shared_cams)
      
      # Matched spatial_cov (same denominator = n_matched)
      n_det_win  <- sum(shared_cams %in% cam_det_ids_win)
      n_det_full <- sum(shared_cams %in% cam_det_ids_full)
      scov_w  <- n_det_win  / n_matched
      scov_f  <- n_det_full / n_matched
      
      # Matched rate = total events / total effort on shared cameras
      rate_w <- sum(win_m$cam_events)  / sum(win_m$cam_effort)
      rate_f <- sum(full_m$cam_events) / sum(full_m$cam_effort)
      
      # SEs from camera-level rate variability on matched cameras
      cam_rates_w <- win_m$cam_rate
      cam_rates_f <- full_m$cam_rate
      n_finite_w <- sum(is.finite(cam_rates_w))
      n_finite_f <- sum(is.finite(cam_rates_f))
      rate_se_w <- if (n_finite_w >= 2) sd(cam_rates_w, na.rm = TRUE) / sqrt(n_finite_w) else NA_real_
      rate_se_f <- if (n_finite_f >= 2) sd(cam_rates_f, na.rm = TRUE) / sqrt(n_finite_f) else NA_real_
      
      # Matched detection rate from matched spatial_cov:
      # P(detect at camera in L days) ≈ scov → rate = -log(1 - scov) / L
      matched_rate_w <- -log(1 - pmin(scov_w, 1 - 1e-6)) / win_len
      matched_rate_f <- -log(1 - pmin(scov_f, 1 - 1e-6)) / 365
      
      tibble(
        n_matched            = as.integer(n_matched),
        spatial_cov_m        = scov_w,
        spatial_cov_m_se     = sqrt(scov_w * (1 - scov_w) / n_matched),
        matched_rate         = matched_rate_w,
        rate_m               = rate_w,
        rate_m_se            = rate_se_w,
        log_rate_m           = log(rate_w + 1e-6),
        log_rate_m_se        = if (!is.na(rate_se_w)) rate_se_w / (rate_w + 1e-6) else NA_real_,
        spatial_cov_full_m   = scov_f,
        spatial_cov_full_m_se = sqrt(scov_f * (1 - scov_f) / n_matched),
        matched_rate_full    = matched_rate_f,
        rate_full_m          = rate_f,
        rate_full_m_se       = rate_se_f,
        log_rate_full_m      = log(rate_f + 1e-6),
        log_rate_full_m_se   = if (!is.na(rate_se_f)) rate_se_f / (rate_f + 1e-6) else NA_real_
      )
    }
    
    # Only compute for non-FULL windows
    matched_idx <- which(window_species$window_id != "FULL")
    
    if (length(matched_idx) > 0) {
      matched_results <- purrr::pmap_dfr(
        list(
          window_species$cam_rate_tbl[matched_idx],
          window_species$cam_det_ids[matched_idx],
          window_species$cam_rate_tbl_full[matched_idx],
          window_species$cam_det_ids_full[matched_idx],
          window_species$window_len[matched_idx]
        ),
        compute_matched
      )
      
      # Bind matched columns to the non-FULL rows, keep FULL rows with NAs
      full_idx <- which(window_species$window_id == "FULL")
      na_matched <- tibble(
        n_matched = NA_integer_, spatial_cov_m = NA_real_,
        spatial_cov_m_se = NA_real_, matched_rate = NA_real_,
        rate_m = NA_real_, rate_m_se = NA_real_,
        log_rate_m = NA_real_, log_rate_m_se = NA_real_,
        spatial_cov_full_m = NA_real_, spatial_cov_full_m_se = NA_real_,
        matched_rate_full = NA_real_, rate_full_m = NA_real_,
        rate_full_m_se = NA_real_,
        log_rate_full_m = NA_real_, log_rate_full_m_se = NA_real_
      )
      all_matched <- bind_rows(
        na_matched[rep(1, length(full_idx)), ],
        matched_results
      )
      # Reorder to match window_species row order
      row_order <- order(c(full_idx, matched_idx))
      all_matched <- all_matched[row_order, ]
      
      window_species <- bind_cols(window_species, all_matched)
      
      # Matched-camera deviations and MSE
      window_species <- window_species %>%
        mutate(
          d_matched_rate  = matched_rate - matched_rate_full,
          d_rate_m        = rate_m       - rate_full_m,
          # matched_rate SE via delta method: se = se_scov / ((1 - scov) * L)
          matched_rate_se = ifelse(!is.na(spatial_cov_m_se) & spatial_cov_m < (1 - 1e-6),
                                  spatial_cov_m_se / ((1 - spatial_cov_m) * window_len),
                                  NA_real_),
          mse_matched_rate = d_matched_rate^2 +
            ifelse(!is.na(matched_rate_se), matched_rate_se^2, NA_real_),
          mse_rate_m = d_rate_m^2 +
            ifelse(!is.na(rate_m_se), rate_m_se^2, NA_real_)
        )
    }
    
    # ---- 8b) Rank preservation (Spearman ρ) against FULL ----
    full_vals <- window_species %>%
      filter(window_id == "FULL") %>%
      select(species,
             lambda_full_rp      = lambda,
             log_rate_full_rp    = log_rate,
             spatial_cov_full_rp = spatial_cov)
    
    safe_spearman <- function(x, y) {
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) >= 3L) cor(x[ok], y[ok], method = "spearman") else NA_real_
    }
    
    rank_pres <- window_species %>%
      filter(window_id != "FULL") %>%
      inner_join(full_vals, by = "species") %>%
      group_by(window_id) %>%
      summarise(
        n_shared_spp    = n(),
        rho_lambda      = safe_spearman(lambda, lambda_full_rp),
        rho_log_rate    = safe_spearman(log_rate, log_rate_full_rp),
        rho_spatial_cov = safe_spearman(spatial_cov, spatial_cov_full_rp),
        .groups = "drop"
      )
    
    window_richness <- window_richness %>%
      left_join(rank_pres, by = "window_id")
    
    # Clean up: remove FULL rows and drop _full helper columns and list columns
    window_species <- window_species %>%
      filter(window_id != "FULL") %>%
      select(
        -spatial_cov_full, -spatial_cov_se_full,
        -lambda_full, -lambda_se_full,
        -rate_full, -rate_se_full,
        -cam_det_ids_full, -cam_rate_tbl_full
      )
  } else if (nrow(window_species)) {
    window_species <- window_species %>%
      mutate(
        d_lambda         = NA_real_,
        d_rate           = NA_real_,
        mse_lambda       = NA_real_,
        mse_rate         = NA_real_,
        n_matched        = NA_integer_,
        spatial_cov_m    = NA_real_,
        spatial_cov_m_se = NA_real_,
        matched_rate     = NA_real_,
        rate_m           = NA_real_,
        rate_m_se        = NA_real_,
        log_rate_m       = NA_real_,
        log_rate_m_se    = NA_real_,
        spatial_cov_full_m    = NA_real_,
        spatial_cov_full_m_se = NA_real_,
        matched_rate_full = NA_real_,
        rate_full_m      = NA_real_,
        rate_full_m_se   = NA_real_,
        log_rate_full_m  = NA_real_,
        log_rate_full_m_se = NA_real_,
        d_matched_rate   = NA_real_,
        d_rate_m         = NA_real_,
        matched_rate_se  = NA_real_,
        mse_matched_rate = NA_real_,
        mse_rate_m       = NA_real_
      )
    window_richness <- window_richness %>%
      mutate(
        n_shared_spp    = NA_integer_,
        rho_lambda      = NA_real_,
        rho_log_rate    = NA_real_,
        rho_spatial_cov = NA_real_
      )
  }
  
  list(
    window_species   = window_species,
    window_camera    = as_tibble(window_camera),
    dropped_species  = if (nrow(dropped_species)) dropped_species else tibble(),
    window_richness  = window_richness,
    y_mats           = y_mats
  )
}



############################################################
#                          WRAPPER                         #
############################################################
# --------------------------------------------------------------------
# dataset_wrapper1()
# --------------------------------------------------------------------
# Aim
#   End-to-end processing for one dataset archive/folder:
#     1) Load deployments/observations and standardise columns.
#     2) Define one or more 12-month anchors (from `anchors` or auto-slicing).
#     3) Build the window grid in “virtual days” aligned to each anchor.
#     4) Call build_window_metrics_fast1() for each anchor slice.
#     5) Tag outputs with dataset (and “…_slice#” if multiple anchors).
#
# Output
#   A list with elements:
#     $window_species, $dropped_species, $window_camera, $window_richness, $y_mats
#     plus:
#     $full_effort_tbl — one row per dataset slice with trap_days_full, n_sites_full
#     $skipped (optional) — if no anchors were found
# --------------------------------------------------------------------
dataset_wrapper1 <- function(
    archive_path,
    independence_mins  = 30,
    occ_grain_days     = 1,
    min_events         = 20,
    min_occasions_pos  = 5,
    min_sites_pos      = 5,
    keep_species       = spp_keep,
    debug_cases        = NULL
) {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
    library(geosphere)
  })
  
  ds_name <- tools::file_path_sans_ext(basename(archive_path))
  message("— ", ds_name, " —")
  
  # ---- load tables ----
  deploy <- read_csv_src(archive_path, "deployments\\.csv")
  obs    <- read_csv_src(archive_path, "observations\\.csv")
  cols <- set_camtrap_cols(obs, deploy)
  list2env(cols, envir = environment())
  
  # ---- anchors ----
  anchor_rows <- tibble()
  if (exists("anchors", inherits = TRUE)) {
    anchor_rows <- anchors %>%
      mutate(dataset_lc = tolower(dataset)) %>%
      filter(dataset_lc == tolower(ds_name)) %>%
      select(-dataset_lc)
  }
  
  if (nrow(anchor_rows) == 0) {
    anchor_rows <- find_anchors(
      archive_path,
      max_slices   = 20,
      keep_species = keep_species,
      verbose      = FALSE
    )
  }
  
  if (nrow(anchor_rows) == 0) {
    message("Skipping (no anchors): ", ds_name)
    return(list(
      window_species   = tibble(),
      dropped_species  = tibble(),
      window_camera    = tibble(),
      window_richness  = tibble(),
      y_mats           = list(),
      full_effort_tbl  = tibble(),
      skipped          = tibble(dataset = ds_name, reason = "no_anchors")
    ))
  }
  
  # ---- per-slice runner ----
  slice_runner <- function(a_start, a_end, slice_id) {
    suffix <- if (nrow(anchor_rows) == 1) "" else paste0("_slice", slice_id)
    ds_tag <- paste0(ds_name, suffix)
    
    deploy_tbl <- deploy %>%
      transmute(
        camera_id = .data[[depl_id_deploy]],
        start     = parse_ts_safe(.data[[start_col]]),
        end       = parse_ts_safe(.data[[end_col]])
      ) %>%
      filter(end > start) %>%
      mutate(start = pmax(start, a_start),
             end   = pmin(end,   a_end)) %>%
      filter(end > start) %>%
      mutate(
        start_virtual = as.numeric(difftime(start, a_start, units = "days")) + 1,
        end_virtual   = as.numeric(difftime(end,   a_start, units = "days")) + 1
      )
    
    full_eff <- compute_full_effort(deploy_tbl, a_start)
    
    # Study geometry (array size + latitude)
    loc_tbl <- deploy %>%
      transmute(
        camera_id = .data[[depl_id_deploy]],
        lat = .data[[lat_col]],
        lon = .data[[lon_col]]
      ) %>%
      distinct() %>%
      drop_na()
    
    dist_mat      <- geosphere::distm(loc_tbl[, c("lon", "lat")], fun = distHaversine) / 1000
    trap_array_km <- max(dist_mat)
    latitude_deg  <- geosphere::centroid(loc_tbl[, c("lon", "lat")])[2]
    
    # ---- observations (cropped to anchor + keep list) ----
    obs_tbl <- obs %>%
      transmute(
        camera_id = .data[[depl_id_obs]],
        species   = stringr::str_squish(stringr::str_trim(as.character(.data[[sci_col]]))),
        timestamp = parse_ts_safe(.data[[date_col]])
      ) %>%
      filter(!is.na(camera_id), !is.na(species), nzchar(species), !is.na(timestamp),
             timestamp >= a_start, timestamp <= a_end) %>%
      { if (!is.null(keep_species)) filter(., species %in% keep_species) else . }
    
    # ---- build window grid in virtual-day space for this anchor ----
    anchor_date_for_doy <- function(doy, a_start) {
      doy <- as.integer(doy)
      y0  <- lubridate::year(a_start)
      d0  <- lubridate::yday(a_start)
      base_year <- ifelse(doy >= d0, y0, y0 + 1)
      as.Date(lubridate::make_date(base_year, 1, 1) + lubridate::days(doy - 1L))
    }
    to_virtual <- function(doy, a_start) {
      ((as.integer(doy) - lubridate::yday(a_start)) %% 365L) + 1L
    }
    
    window_grid_ds <- window_grid %>%
      mutate(
        start_virtual = to_virtual(start_doy, a_start),
        end_virtual   = start_virtual + length_d - 1L,
        start_date    = anchor_date_for_doy(start_doy, a_start),
        end_doy_cal   = ((start_doy - 1L + length_d - 1L) %% 365L) + 1L,
        end_date      = anchor_date_for_doy(end_doy_cal, a_start)
      ) %>%
      select(-end_doy_cal)
    
    if (!is.null(debug_cases) && nrow(debug_cases)) {
      window_grid_ds <- window_grid_ds %>%
        filter(window_id %in% unique(debug_cases$window_id))
    }
    
    core <- build_window_metrics_fast1(
      window_grid       = window_grid_ds,
      deploy_tbl        = deploy_tbl,
      obs_tbl           = obs_tbl,
      a_start           = a_start,
      a_end             = a_end,
      independence_mins = independence_mins,
      occ_grain_days    = occ_grain_days,
      min_events        = min_events,
      min_occasions_pos = min_occasions_pos,
      min_sites_pos     = min_sites_pos,
      trap_array_km     = trap_array_km,
      latitude_deg      = latitude_deg,
      drop_leap_day     = TRUE,
      debug_cases       = debug_cases
    )
    
    # Attach dataset-level constants (works even for empty tibbles)
    core$window_species  <- core$window_species  %>%
      mutate(trap_days_full = full_eff$trap_days_full,
             n_sites_full   = full_eff$n_sites_full,
             dataset        = ds_tag)
    
    core$dropped_species <- core$dropped_species %>%
      mutate(trap_days_full = full_eff$trap_days_full,
             n_sites_full   = full_eff$n_sites_full,
             dataset        = ds_tag)
    
    core$window_camera   <- core$window_camera   %>%
      mutate(trap_days_full = full_eff$trap_days_full,
             n_sites_full   = full_eff$n_sites_full,
             dataset        = ds_tag)
    
    core$window_richness <- core$window_richness %>%
      mutate(trap_days_full = full_eff$trap_days_full,
             n_sites_full   = full_eff$n_sites_full,
             dataset        = ds_tag)
    
    # One-row summary per slice (robust even if window_species is empty)
    core$full_effort_tbl <- tibble(
      dataset        = ds_tag,
      trap_days_full = full_eff$trap_days_full,
      n_sites_full   = full_eff$n_sites_full
    )
    
    if (!is.null(core$y_mats) && length(core$y_mats)) {
      names(core$y_mats) <- paste0(ds_tag, "__", names(core$y_mats))
    }
    
    core
  }
  
  cores <- purrr::pmap(
    list(anchor_rows$anchor_start,
         anchor_rows$anchor_end,
         as.integer(anchor_rows$slice)),
    slice_runner
  )
  
  res <- if (length(cores) == 1) {
    cores[[1]]
  } else {
    list(
      window_species   = purrr::map_df(cores, "window_species"),
      dropped_species  = purrr::map_df(cores, "dropped_species"),
      window_camera    = purrr::map_df(cores, "window_camera"),
      window_richness  = purrr::map_df(cores, "window_richness"),
      y_mats           = purrr::flatten(purrr::map(cores, ~ .x$y_mats %||% list())),
      full_effort_tbl  = purrr::map_df(cores, "full_effort_tbl")
    )
  }
  
  # If single-slice, full_effort_tbl is already present; if not, ensure it exists
  if (is.null(res$full_effort_tbl)) {
    res$full_effort_tbl <- tibble()
  }
  
  res
}





# ────────────────────────────────────────────────────────────────
# USAGE
# ────────────────────────────────────────────────────────────────

t_start <- proc.time()
message(sprintf("[%s] Starting dataset_wrapper1 across %d datasets …",
                format(Sys.time(), "%H:%M:%S"), length(ds_paths)))

wrapped <- furrr::future_map(
  ds_paths,
  dataset_wrapper1,
  .options = furrr::furrr_options(seed = TRUE)  # default scheduling=1: tasks distributed evenly
)

t_elapsed <- proc.time() - t_start
message(sprintf("[%s] Finished. Elapsed: %.1f min (%.0f s)",
                format(Sys.time(), "%H:%M:%S"),
                t_elapsed["elapsed"] / 60,
                t_elapsed["elapsed"]))

saveRDS(wrapped, "wrapped_cam.rds")

wrapped <- readRDS("wrapped_cam.rds")
all_window_species  <- map_df(wrapped, "window_species")
all_window_richness  <- map_df(wrapped, "window_richness")
all_dropped_species <- map_df(wrapped, "dropped_species")

saveRDS(all_window_species, "all_window_species.rds")
saveRDS(all_window_richness, "all_window_richness.rds")
saveRDS(all_dropped_species, "all_dropped_species.rds")



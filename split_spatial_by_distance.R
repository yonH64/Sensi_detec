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
  
  # ---- cluster sizes (stations + deployments) --------------------------------
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
  
  invisible(created)
}


# Default: keep clusters that have at least 10 stations
split_spatial_by_distance("NL-MICA",
                          export_dir   = "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/_splits",
                          threshold_km = 30,
                          min_stations = 10)

split_spatial_by_distance("SE-Tim1",
                          export_dir   = "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/_splits",
                          threshold_km = 15,
                          min_stations = 18)

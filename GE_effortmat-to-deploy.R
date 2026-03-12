fix_wide_deployments_to_camtrapdp <- function(dataset_name,
                                              root,
                                              include_zips = TRUE,
                                              date_order = c("auto","dmy","mdy"),
                                              tz_out = "UTC",
                                              backup = TRUE,
                                              overwrite = TRUE,
                                              verbose = TRUE) {
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(stringr); library(lubridate); library(readr)
  })
  date_order <- match.arg(date_order)
  
  # ---- locate dataset path (must be a folder if overwriting) ----------
  archive_path <- resolve_dataset_path(dataset_name, root, include_zips = include_zips)
  if (!dir.exists(archive_path)) {
    stop("Dataset '", dataset_name, "' is an archive. Can't overwrite inside archives. Unzip first.")
  }
  
  # Find actual file paths (for overwrite)
  files <- list.files(archive_path, recursive = TRUE, full.names = TRUE)
  depl_path <- files[str_detect(basename(files), regex("^deployments\\.csv$", ignore_case = TRUE))][1]
  obs_path  <- files[str_detect(basename(files), regex("^observations\\.csv$", ignore_case = TRUE))][1]
  if (is.na(depl_path)) stop("No deployments.csv found for dataset: ", dataset_name)
  
  # ---- read (suppress "New names" messages) --------------------------
  deploy_raw <- suppressMessages(read_csv_src(archive_path, "deployments\\.csv"))
  if (all(c("deploymentID","locationID","deploymentStart","deploymentEnd") %in% names(deploy_raw))) {
    if (verbose) message("ℹ️ ", dataset_name, ": deployments already camtrapDP-like; skipping.")
    return(invisible(FALSE))
  }
  
  names(deploy_raw) <- trimws(names(deploy_raw))
  deploy_raw <- deploy_raw %>% select(-any_of(c("...1","X1")))
  
  # ---- identify camera/location id column ----------------------------
  id_candidates <- c("cameraID","camera_id","locationID","location_id","location.id",
                     "station","station_id","site","site_id")
  id_col <- intersect(id_candidates, names(deploy_raw))[1]
  if (is.na(id_col)) id_col <- names(deploy_raw)[1]
  
  # ---- parse date out of colname robustly ----------------------------
  parse_date_from_colname <- function(nm, date_order_resolved) {
    x <- nm
    x <- str_remove(x, "\\.\\.\\.\\d+$")     # remove readr suffix ...2
    x <- str_remove(x, "^X")                # remove leading X
    x <- str_replace_all(x, "[^0-9]", "_")   # normalize separators to _
    x <- str_replace_all(x, "_+", "_")
    x <- str_replace_all(x, "^_|_$", "")
    parts <- unlist(str_split(x, "_"))
    if (length(parts) < 3) return(NA_Date_)
    parts <- parts[1:3]
    
    a <- parts[1]; b <- parts[2]; c <- parts[3]
    
    # YMD if first token looks like year
    if (nchar(a) == 4) {
      yy <- a; mm <- b; dd <- c
    } else {
      yy <- c
      if (nchar(yy) == 2) yy <- paste0("20", yy)  # adjust if you have 1990s etc.
      if (date_order_resolved == "dmy") { dd <- a; mm <- b } else { mm <- a; dd <- b }
    }
    
    dd <- suppressWarnings(as.integer(dd))
    mm <- suppressWarnings(as.integer(mm))
    yy <- suppressWarnings(as.integer(yy))
    if (anyNA(c(dd, mm, yy))) return(NA_Date_)
    if (mm < 1 || mm > 12 || dd < 1 || dd > 31 || yy < 1900 || yy > 2100) return(NA_Date_)
    
    as.Date(sprintf("%04d-%02d-%02d", yy, mm, dd))
  }
  
  resolve_date_order <- function(nms) {
    cleaned <- nms |>
      str_remove("\\.\\.\\.\\d+$") |>
      str_remove("^X") |>
      str_replace_all("[^0-9]", "_") |>
      str_replace_all("_+", "_") |>
      str_replace_all("^_|_$", "")
    
    spl <- str_split_fixed(cleaned, "_", 3)
    p1 <- spl[,1]; p2 <- spl[,2]
    is_year_first <- nchar(p1) == 4
    p1n <- suppressWarnings(as.integer(p1[!is_year_first]))
    p2n <- suppressWarnings(as.integer(p2[!is_year_first]))
    
    if (length(p1n) == 0 || length(p2n) == 0) return("dmy")
    if (any(p1n > 12, na.rm = TRUE)) return("dmy")
    if (any(p2n > 12, na.rm = TRUE)) return("mdy")
    "dmy"
  }
  
  # Candidate columns = everything except id_col; keep only those that parse to valid dates
  candidate_cols <- setdiff(names(deploy_raw), id_col)
  
  date_order_resolved <- if (date_order == "auto") {
    out <- resolve_date_order(candidate_cols)
    if (verbose) message("ℹ️ ", dataset_name, ": date_order=auto resolved to '", out, "'.")
    out
  } else date_order
  
  # IMPORTANT FIX: preserve Date class (no vapply here)
  parsed_dates <- do.call(c, lapply(candidate_cols, parse_date_from_colname, date_order_resolved = date_order_resolved))
  # parsed_dates is a Date vector; names align with candidate_cols by position
  
  date_cols_orig <- candidate_cols[!is.na(parsed_dates)]
  if (length(date_cols_orig) < 5) {
    stop(dataset_name, ": could not identify enough date columns after normalization/parsing. ",
         "Check deployments.csv header format.")
  }
  
  # ---- rename date columns to canonical dYYYY_MM_DD (+ __dupN if needed)
  canon_base  <- paste0("d", format(parsed_dates[match(date_cols_orig, candidate_cols)], "%Y_%m_%d"))
  dup_n       <- ave(canon_base, canon_base, FUN = seq_along)
  canon_names <- ifelse(dup_n == 1, canon_base, paste0(canon_base, "__dup", dup_n))
  
  # IMPORTANT: rename() wants new = old
  rename_map  <- stats::setNames(date_cols_orig, canon_names)
  
  deploy_norm <- dplyr::rename(deploy_raw, !!!rename_map)
  date_cols_new <- canon_names
  
  # ---- build per-camera daily series (dataset-specific date ranges) ---
  long_daily <- deploy_norm %>%
    mutate(cam_id = str_trim(as.character(.data[[id_col]])),
           cam_id = str_replace(cam_id, "\\.0+$", "")) %>%  # helps if ids came as 201.0
    select(cam_id, all_of(date_cols_new)) %>%
    pivot_longer(cols = all_of(date_cols_new), names_to = "date_str", values_to = "effort") %>%
    mutate(
      date_base = str_remove(date_str, "__dup\\d+$"),
      date = as.Date(str_remove(date_base, "^d"), format = "%Y_%m_%d"),
      effort_num = suppressWarnings(as.numeric(effort)),
      effort_chr = toupper(trimws(as.character(effort))),
      has_effort = case_when(
        !is.na(effort_num) ~ effort_num > 0,
        effort_chr %in% c("TRUE","T","YES","Y","1") ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    filter(!is.na(cam_id), nzchar(cam_id), !is.na(date)) %>%
    group_by(cam_id, date) %>%
    summarise(has_effort = any(has_effort), .groups = "drop") %>%
    group_by(cam_id) %>%
    complete(date = seq(min(date), max(date), by = "day"),
             fill = list(has_effort = FALSE)) %>%
    arrange(cam_id, date) %>%
    ungroup()
  
  # ---- split into bouts: every no-effort day breaks deployment --------
  bouts <- long_daily %>%
    group_by(cam_id) %>%
    mutate(
      start_bout = has_effort & !lag(has_effort, default = FALSE),
      bout_id = cumsum(start_bout)
    ) %>%
    filter(has_effort) %>%
    group_by(cam_id, bout_id) %>%
    summarise(
      start_date = min(date),
      end_date   = max(date),
      .groups = "drop"
    ) %>%
    mutate(
      locationID   = cam_id,   # camera id as location ID (your requirement)
      locationName = cam_id,
      cameraID     = cam_id,
      deploymentID = paste0(cam_id, "_", sprintf("%02d", bout_id), "_", format(start_date, "%Y%m%d")),
      deploymentStart = as.POSIXct(start_date, tz = tz_out),
      deploymentEnd   = as.POSIXct(end_date + days(1), tz = tz_out) - seconds(1)
    )
  
  if (!nrow(bouts)) stop(dataset_name, ": no effort days >0/TRUE found; nothing to convert.")
  
  # ---- join coordinates from observations.csv by location_id ----------
  coords <- tibble(locationID = character(), latitude = numeric(), longitude = numeric())
  if (!is.na(obs_path)) {
    obs <- suppressMessages(read_csv_src(archive_path, "observations\\.csv"))
    names(obs) <- trimws(names(obs))
    
    # you asked specifically for 'location_id' matching
    if (!("location_id" %in% names(obs))) {
      if (verbose) message("⚠️ ", dataset_name, ": observations.csv has no 'location_id' column; leaving lat/lon NA.")
    } else {
      # flexible lat/lon detection in observations
      lat_col <- names(obs)[str_detect(names(obs), regex("(^lat$|latitude)", ignore_case = TRUE))][1]
      lon_col <- names(obs)[str_detect(names(obs), regex("(^lon$|longitude|lng|long)", ignore_case = TRUE))][1]
      
      if (is.na(lat_col) || is.na(lon_col)) {
        if (verbose) message("⚠️ ", dataset_name, ": observations.csv missing lat/lon columns; leaving lat/lon NA.")
      } else {
        coords <- obs %>%
          transmute(
            locationID = str_trim(as.character(.data[["location_id"]])),
            locationID = str_replace(locationID, "\\.0+$", ""),
            latitude   = parse_num_safe(.data[[lat_col]]),
            longitude  = parse_num_safe(.data[[lon_col]])
          ) %>%
          filter(!is.na(locationID), nzchar(locationID), is.finite(latitude), is.finite(longitude)) %>%
          group_by(locationID) %>%
          summarise(latitude = mean(latitude, na.rm = TRUE),
                    longitude = mean(longitude, na.rm = TRUE),
                    .groups = "drop")
      }
    }
  } else if (verbose) {
    message("⚠️ ", dataset_name, ": observations.csv not found; leaving lat/lon NA.")
  }
  
  bouts2 <- bouts %>% left_join(coords, by = "locationID")
  
  # ---- format timestamps to ISO --------------------------------------
  fmt_iso <- function(x) {
    if (toupper(tz_out) == "UTC") return(format(x, "%Y-%m-%dT%H:%M:%SZ"))
    s <- format(x, "%Y-%m-%dT%H:%M:%S%z")
    sub("([+-]\\d{2})(\\d{2})$", "\\1:\\2", s)
  }
  
  deployments_camtrapdp <- bouts2 %>%
    mutate(
      deploymentStart = fmt_iso(deploymentStart),
      deploymentEnd   = fmt_iso(deploymentEnd)
    ) %>%
    transmute(
      deploymentID,
      locationID,
      locationName,
      latitude = latitude,
      longitude = longitude,
      coordinateUncertainty = NA_real_,
      deploymentStart,
      deploymentEnd,
      setupBy = NA_character_,
      cameraID,
      cameraModel = NA_character_,
      cameraDelay = NA_real_,
      cameraHeight = NA_real_,
      cameraDepth = NA_real_,
      cameraTilt = NA_real_,
      cameraHeading = NA_real_,
      detectionDistance = NA_real_,
      timestampIssues = NA_character_,
      baitUse = NA_character_,
      featureType = NA_character_,
      habitat = NA_character_,
      deploymentGroups = NA_character_,
      deploymentTags = NA_character_,
      deploymentComments = NA_character_
    )
  
  # ---- backup + overwrite --------------------------------------------
  if (backup) {
    # WITHOUT "deployments" in its name (your requirement)
    bkp <- file.path(dirname(depl_path), paste0("backup_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
    file.copy(depl_path, bkp, overwrite = TRUE)
    if (verbose) message("🗂️ Backup written: ", bkp)
  }
  
  if (overwrite) {
    write_csv(deployments_camtrapdp, depl_path, na = "")
    if (verbose) message("✅ Overwrote deployments.csv for ", dataset_name, " (", nrow(deployments_camtrapdp), " deployments)")
    return(invisible(TRUE))
  } else {
    return(invisible(deployments_camtrapdp))
  }
}


root <- "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets"
fix_wide_deployments_to_camtrapdp("GE-Berchtesgaden_NP", root = root, date_order = "auto")
fix_wide_deployments_to_camtrapdp("GE-Black_Forest_NP", root = root)
fix_wide_deployments_to_camtrapdp("GE-Eifel_NP", root = root)
fix_wide_deployments_to_camtrapdp("GE-Hainich_NP", root = root)
fix_wide_deployments_to_camtrapdp("GE-Harz_NP", root = root)
fix_wide_deployments_to_camtrapdp("GE-Hunsrueck_Hochwald_NP", root = root)
fix_wide_deployments_to_camtrapdp("GE-Kellerwald_Edersee_NP", root = root)
fix_wide_deployments_to_camtrapdp("GE-Mueritz_NP", root = root)

raw <- readxl::read_excel("C:/Users/yonah/Downloads/CAM_all.xlsx") 
  ############################################################
#                 EXCEL  ->  CamtrapDP CSVs                #
############################################################

# ─────────────────────────────────────────────────────────
# Libraries
# ─────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(readr)
  library(camtrapdp)
})

############################################################
#                           HELPERS                        #
############################################################

# --------------------------------------------------------------------
# iso_dt()
# --------------------------------------------------------------------
# Aim
#   Format POSIXct datetimes consistently for CSV export.
#
# Arguments
#   x  : POSIXct vector
#
# Output
#   Character vector "YYYY-mm-dd HH:MM:SS"
# --------------------------------------------------------------------
iso_dt <- function(x) {
  format(as.POSIXct(x, tz = "UTC"), "%Y-%m-%d %H:%M:%S")
}

# --------------------------------------------------------------------
# ensure_camtrap_cols()
# --------------------------------------------------------------------
# Aim
#   Ensure a tibble has exactly the CamtrapDP columns (add missing as NA,
#   drop extras, and order columns).
#
# Arguments
#   df       : input tibble
#   colnames : target column names (character vector)
#
# Output
#   Tibble with exactly 'colnames' in that order
# --------------------------------------------------------------------
ensure_camtrap_cols <- function(df, colnames) {
  missing <- setdiff(colnames, names(df))
  if (length(missing)) {
    for (nm in missing) df[[nm]] <- NA
  }
  df %>% select(any_of(colnames))
}

############################################################
#                           MAIN                           #
############################################################

# --------------------------------------------------------------------
# excel_to_camtrapdp_csv()
# --------------------------------------------------------------------
# Aim
#   Convert your Excel event table to CamtrapDP-compatible deployments
#   and observations CSVs (no media.csv).
#
# Arguments
#   xlsx_path : path to the Excel file
#   out_dir   : where to write deployments.csv and observations.csv
#   sheet     : sheet name or index
#   tz        : timezone to assume for 'date' (default "UTC")
#
# Output
#   Writes:
#     out_dir/deployments.csv
#     out_dir/observations.csv
#   Returns a list(deployments=..., observations=...)
# --------------------------------------------------------------------
excel_to_camtrapdp_csv <- function(xlsx_path,
                                   out_dir = dirname(xlsx_path),
                                   sheet = 1,
                                   tz = "UTC") {
  
  # ---- CamtrapDP column templates (keeps you aligned with the spec)
  tmpl <- camtrapdp::example_dataset()
  dep_cols <- names(camtrapdp::deployments(tmpl))
  obs_cols <- names(camtrapdp::observations(tmpl))
  
  # ---- read Excel
  raw <- readxl::read_excel(xlsx_path) %>%
    rename_with(~ str_replace_all(tolower(.x), "\\s+", "_"))
  
  # ---- basic checks (fail early, loudly)
  need <- c("n_cam", "date", "seqnum", "obs", "count", "position")
  miss <- setdiff(need, names(raw))
  if (length(miss)) stop("Missing required columns in Excel: ", paste(miss, collapse = ", "))
  
  # ---- parse datetime + coordinates
  dat <- raw %>%
    mutate(
      # Excel often gives POSIXct already; force tz handling
      date = as.POSIXct(date, tz = tz),
      
      obs = str_to_lower(str_trim(as.character(obs))),
      position = str_trim(as.character(position))
    ) %>%
    separate(
      position,
      into = c("latitude", "longitude"),
      sep = "\\s*,\\s*",
      convert = TRUE,
      remove = FALSE
    )
  
  # ---- species mapping (EDIT THIS if needed)
  # Rules:
  #  - values set to NA will be dropped by default (non-species or non-target)
  #  - you can change those NA to a taxon name if you prefer keeping them
  species_lut <- c(
    blaireau      = "Meles meles",
    chat          = "Felis catus",
    chevre        = "Capra hircus",
    chevreuil     = "Capreolus capreolus",
    chien         = "Canis lupus familiaris",
    ecureuil      = "Sciurus vulgaris",
    equide        = "Equus caballus",
    genette       = "Genetta genetta",
    herisson      = "Erinaceus europaeus",
    ragondin      = "Myocastor coypus",
    renard        = "Vulpes vulpes",
    sanglier      = "Sus scrofa",
    
    # drop these by default (not species-level or not useful for your workflow)
    humain        = NA_character_,
    vehicule      = NA_character_,
    vide          = NA_character_,
    oiseau        = NA_character_,
    micromammifere= NA_character_,
    lagomorphe    = NA_character_,
    mustelide     = NA_character_
  )
  
  dat <- dat %>%
    mutate(
      scientificName = dplyr::recode(obs, !!!species_lut, .default = obs),
      count = suppressWarnings(as.integer(count)),
      count = if_else(is.na(count) | count < 1L, 1L, count)
    )
  
  # ---- show you what is NOT covered by the mapping
  unknown_labels <- dat %>%
    distinct(obs, scientificName) %>%
    filter(scientificName == obs) %>%
    pull(obs) %>%
    sort()
  
  if (length(unknown_labels)) {
    message("Unmapped 'obs' labels found (kept as-is in scientificName): ",
            paste(unknown_labels, collapse = ", "))
  }
  
  # ---- drop rows that became NA after mapping (vide/humain/etc by default)
  dat_keep <- dat %>% filter(!is.na(scientificName))
  
  # ==========================================================
  # DEPLOYMENTS
  # ==========================================================
  # Assumption:
  #   1 deployment per camera ID (n_cam) using min/max 'date' seen.
  # If your cameras were moved, tell me and we’ll split deployments
  # when 'position' changes (or when long gaps occur).
  deployments <- dat %>%
    group_by(n_cam) %>%
    summarise(
      deploymentStart = min(date, na.rm = TRUE),
      deploymentEnd   = max(date, na.rm = TRUE),
      latitude  = dplyr::first(latitude),
      longitude = dplyr::first(longitude),
      .groups = "drop"
    ) %>%
    mutate(
      deploymentID = paste0("cam_", n_cam),
      locationID   = deploymentID,
      cameraID     = as.character(n_cam),
      
      # CamtrapDP expects datetimes in character form in CSVs
      deploymentStart = iso_dt(deploymentStart),
      deploymentEnd   = iso_dt(deploymentEnd)
    )
  
  # ==========================================================
  # OBSERVATIONS
  # ==========================================================
  # We collapse multiple rows that share (camera, seqnum, scientificName)
  # into ONE event-level observation:
  observations <- dat_keep %>%
    mutate(deploymentID = paste0("cam_", n_cam),
           eventID = as.character(seqnum)) %>%
    group_by(deploymentID, eventID, scientificName) %>%
    summarise(
      eventStart = min(date, na.rm = TRUE),
      eventEnd   = max(date, na.rm = TRUE),
      count      = max(count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      observationID = paste0("obs_", row_number()),
      
      # Strong defaults (edit if you want)
      observationLevel = "event",
      observationType  = "animal",
      classificationMethod = "human",
      
      eventStart = iso_dt(eventStart),
      eventEnd   = iso_dt(eventEnd)
    )
  
  # ---- conform to CamtrapDP column sets
  deployments_out  <- ensure_camtrap_cols(deployments, dep_cols)
  observations_out <- ensure_camtrap_cols(observations, obs_cols)
  
  # ---- write files
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  readr::write_csv(deployments_out,  file.path(out_dir, "deployments.csv"),  na = "")
  readr::write_csv(observations_out, file.path(out_dir, "observations.csv"), na = "")
  
  message("✓ Wrote deployments.csv (", nrow(deployments_out), " rows)")
  message("✓ Wrote observations.csv (", nrow(observations_out), " rows)")
  message("Date range (from observations): ",
          min(dat_keep$date, na.rm = TRUE), " → ", max(dat_keep$date, na.rm = TRUE))
  
  invisible(list(deployments = deployments_out, observations = observations_out))
}

# ─────────────────────────────────────────────────────────
# RUN
# ─────────────────────────────────────────────────────────
# Example:
excel_to_camtrapdp_csv("C:/Users/yonah/Downloads/CAM_all.xlsx", 
                       out_dir = "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/FR-montpellier")

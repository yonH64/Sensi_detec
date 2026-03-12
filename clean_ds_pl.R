# ============================================================
# Convert deployments1.csv + observations1.csv to CamtrapDP-like
# schema matching deployments0.csv / observations0.csv
# and OVERWRITE the existing *1.csv files.
#
# NOTE: ignores the missing 75_22 deployment for now.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(lubridate)
  library(tidyr)
  library(purrr)
})

# ---- paths (edit if needed) ---------------------------------
dir_in <- "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/PL-kampinos_NP"  # folder where the CSVs live
deploy_in <- file.path(dir_in, "deployments.csv")
obs_in    <- file.path(dir_in, "observations.csv")

# Optional safety backup (comment out if you truly want no backups)
file.copy(deploy_in, paste0(deploy_in, ".bak"), overwrite = TRUE)
file.copy(obs_in,    paste0(obs_in,    ".bak"), overwrite = TRUE)

# ---- helpers ------------------------------------------------
make_id <- function(prefix, n) sprintf("%s_%06d", prefix, seq_len(n))

clean_location_id <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\\\", "/")          # windows paths -> /
  x <- str_split_fixed(x, "/", 2)[, 1]          # keep first token
  x <- str_replace_all(x, "-", "_")             # 105-22 -> 105_22
  str_trim(x)
}

parse_ts <- function(x) {
  out <- suppressWarnings(ymd_hms(x, tz = "UTC", quiet = TRUE))
  if (all(is.na(out))) out <- suppressWarnings(ymd(x, tz = "UTC", quiet = TRUE))
  out
}

pick_col <- function(df, candidates, required = TRUE) {
  nms <- names(df)
  # exact match first
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0) return(hit[1])
  
  # case-insensitive match
  low <- tolower(nms)
  cand_low <- tolower(candidates)
  idx <- match(cand_low, low)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(nms[idx[1]])
  
  if (required) {
    stop("Could not find any of these columns: ",
         paste(candidates, collapse = ", "),
         "\nAvailable columns are: ",
         paste(nms, collapse = ", "))
  } else {
    return(NA_character_)
  }
}

as_num <- function(x) suppressWarnings(as.numeric(x))

# ---- target schemas (must match *0.csv) ---------------------
deploy_cols <- c(
  "deploymentID", "locationID", "locationName", "latitude", "longitude",
  "coordinateUncertainty", "deploymentStart", "deploymentEnd", "setupBy",
  "cameraID", "cameraModel", "cameraDelay", "cameraHeight", "cameraDepth",
  "cameraTilt", "cameraHeading", "detectionDistance", "timestampIssues",
  "baitUse", "featureType", "habitat", "deploymentGroups",
  "deploymentTags", "deploymentComments"
)

obs_cols <- c(
  "observationID", "deploymentID", "mediaID", "eventID", "eventStart",
  "eventEnd", "observationLevel", "observationType", "cameraSetupType",
  "scientificName", "count", "lifeStage", "sex", "behavior",
  "individualID", "individualPositionRadius", "individualPositionAngle",
  "individualSpeed", "bboxX", "bboxY", "bboxWidth", "bboxHeight",
  "classificationMethod", "classifiedBy", "classificationTimestamp",
  "classificationProbability", "observationTags", "observationComments"
)

# ---- read ---------------------------------------------------
dep1 <- read_csv(deploy_in, show_col_types = FALSE)
obs1 <- read_csv(obs_in,    show_col_types = FALSE)

message("Deployments1 columns: ", paste(names(dep1), collapse = ", "))
message("Observations1 columns: ", paste(names(obs1), collapse = ", "))

# ---- auto-detect column names in deployments1 ----------------
dep_loc_col   <- pick_col(dep1, c("name","Name","locationID","locationId","location","site","Site","station","Station"))
dep_start_col <- pick_col(dep1, c("Start","start","deploymentStart","deployment_start","from","begin","Begin","date_start","start_date"))
dep_end_col   <- pick_col(dep1, c("End","end","deploymentEnd","deployment_end","to","finish","Finish","date_end","end_date"))

# latitude candidates (common exports)
dep_lat_col <- pick_col(dep1, c("Y","y","lat","Lat","latitude","Latitude","LAT","northing","Northing"))
# longitude candidates (common exports) — include Y.1 fallback if that’s what you have
dep_lon_col <- pick_col(dep1, c("X","x","lon","Lon","longitude","Longitude","LON","easting","Easting","Y.1","y.1"))

# ---- DEPLOYMENTS: reshape to CamtrapDP-like ------------------
dep_out <- dep1 %>%
  transmute(
    deploymentID = make_id("dep", n()),
    locationID   = clean_location_id(.data[[dep_loc_col]]),
    locationName = clean_location_id(.data[[dep_loc_col]]),
    latitude     = as_num(.data[[dep_lat_col]]),
    longitude    = as_num(.data[[dep_lon_col]]),
    coordinateUncertainty = NA_real_,
    deploymentStart = parse_ts(.data[[dep_start_col]]),
    deploymentEnd   = parse_ts(.data[[dep_end_col]]),
    setupBy = NA_character_,
    cameraID = NA_character_,
    cameraModel = NA_character_,
    cameraDelay = NA_character_,
    cameraHeight = NA_character_,
    cameraDepth = NA_character_,
    cameraTilt = NA_character_,
    cameraHeading = NA_character_,
    detectionDistance = NA_character_,
    timestampIssues = NA_character_,
    baitUse = NA_character_,
    featureType = NA_character_,
    habitat = NA_character_,
    deploymentGroups = NA_character_,
    deploymentTags = NA_character_,
    deploymentComments = NA_character_
  ) %>%
  mutate(across(setdiff(deploy_cols, names(.)), ~ NA)) %>%
  select(all_of(deploy_cols))

# soft warning if lon == lat everywhere (your known issue)
if (all(is.finite(dep_out$latitude)) && all(is.finite(dep_out$longitude))) {
  if (isTRUE(all.equal(dep_out$latitude, dep_out$longitude, tolerance = 0))) {
    message("⚠️ longitude appears identical to latitude for all rows. ",
            "If this is wrong, replace the longitude column in deployments1 with true X/longitude values.")
  }
}

# ---- auto-detect column names in observations1 ----------------
obs_path_col <- pick_col(obs1, c("RelativePath","relativePath","relative_path","path","Path","folder","Folder"))
obs_time_col <- pick_col(obs1, c("DateTime","datetime","date_time","Timestamp","timestamp","time","Time","Date"))
obs_spp_col  <- pick_col(obs1, c("Species","species","taxon","Taxon","label","Label"))
obs_n_col    <- pick_col(obs1, c("individuals","Individuals","count","Count","n","N"), required = FALSE)

# robust species mapping (lowercase keys)
species_map <- c(
  "human"         = "Homo sapiens",
  "dog-no_leash"  = "Canis lupus familiaris",
  "dog-leash"     = "Canis lupus familiaris",
  "fox"           = "Vulpes vulpes",
  "roe deer"      = "Capreolus capreolus",
  "red deer"      = "Cervus elaphus",
  "wolf"          = "Canis lupus",
  "badger"        = "Meles meles",
  "wild boar"     = "Sus scrofa",
  "moose"         = "Alces alces",
  "lynx"          = "Lynx lynx",
  "raccoon dog"   = "Nyctereutes procyonoides"
)

# optional comment columns (keep whatever exists)
comment_candidates <- c("com to record","Com to record","com_to_record",
                        "Com to species","com to species","com_to_species",
                        "Season","season")
comment_present <- names(obs1)[tolower(names(obs1)) %in% tolower(comment_candidates)]

row_comment <- function(...) {
  vals <- c(...)
  vals <- vals[!is.na(vals) & nzchar(trimws(vals))]
  if (!length(vals)) return(NA_character_)
  str_squish(paste(vals, collapse = " | "))
}

# ---- OBSERVATIONS: reshape + map species ---------------------
obs_core <- obs1 %>%
  mutate(
    locationID = clean_location_id(.data[[obs_path_col]]),
    eventStart = parse_ts(.data[[obs_time_col]]),
    eventEnd   = eventStart,
    .sp_raw = str_to_lower(str_squish(as.character(.data[[obs_spp_col]]))),
    scientificName = recode(.sp_raw, !!!species_map, .default = NA_character_),
    count = if (!is.na(obs_n_col)) suppressWarnings(as.integer(.data[[obs_n_col]])) else NA_integer_,
    observationID = make_id("obs", n()),
    eventID       = make_id("evt", n()),
    observationLevel = "event",
    observationType  = "unknown",
    cameraSetupType  = NA_character_,
    mediaID          = NA_character_,
    lifeStage = NA_character_,
    sex       = NA_character_,
    behavior  = NA_character_,
    individualID = NA_character_,
    individualPositionRadius = NA_real_,
    individualPositionAngle  = NA_real_,
    individualSpeed = NA_real_,
    bboxX = NA_real_, bboxY = NA_real_, bboxWidth = NA_real_, bboxHeight = NA_real_,
    classificationMethod = "human",
    classifiedBy = NA_character_,
    classificationTimestamp = NA_character_,
    classificationProbability = NA_real_,
    observationTags = NA_character_
  ) %>%
  {
    if (length(comment_present)) {
      comments_df <- select(., all_of(comment_present))
      mutate(., observationComments = pmap_chr(comments_df, row_comment))
    } else {
      mutate(., observationComments = NA_character_)
    }
  } %>%
  select(-any_of(c(".sp_raw", comment_present)))

# ---- assign deploymentID by interval match -------------------
dep_for_join <- dep_out %>%
  mutate(deploymentStart = as.POSIXct(deploymentStart, tz = "UTC"),
         deploymentEnd   = as.POSIXct(deploymentEnd,   tz = "UTC"))

obs_for_join <- obs_core %>%
  mutate(eventStart = as.POSIXct(eventStart, tz = "UTC"))

cand <- obs_for_join %>%
  select(observationID, locationID, eventStart) %>%
  inner_join(dep_for_join %>% select(deploymentID, locationID, deploymentStart, deploymentEnd),
             by = "locationID") %>%
  filter(!is.na(eventStart),
         !is.na(deploymentStart), !is.na(deploymentEnd),
         eventStart >= deploymentStart,
         eventStart <= deploymentEnd) %>%
  group_by(observationID) %>%
  slice(1) %>%   # if multiple deployments match, keep first
  ungroup() %>%
  select(observationID, deploymentID)

obs_out <- obs_for_join %>%
  left_join(cand, by = "observationID") %>%
  mutate(across(setdiff(obs_cols, names(.)), ~ NA)) %>%
  select(all_of(obs_cols))

message("Obs rows with no matched deploymentID (expected until missing deployments fixed): ",
        sum(is.na(obs_out$deploymentID)))

# ---- overwrite files ----------------------------------------
write_csv(dep_out, deploy_in, na = "")
write_csv(obs_out, obs_in,    na = "")

message("✅ Overwrote: ", deploy_in)
message("✅ Overwrote: ", obs_in)

library(data.table)
library(lubridate)
library(stringr)

# ============================================================
# FILES (EDIT THIS)
#   - If you put this script inside the dataset folder, use "."
# ============================================================
dataset_dir <- "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/IT-Alps"   # e.g. "C:/Users/yonah/.../Datasets/IT-Alps"
deploy_file <- file.path(dataset_dir, "deployments.csv")
obs_file    <- file.path(dataset_dir, "observations.csv")

stopifnot(file.exists(deploy_file), file.exists(obs_file))

# ============================================================
# BACKUPS (timestamped, so you don't overwrite older backups)
# ============================================================
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
deploy_bak <- sub("\\.csv$", paste0("_backup_", stamp, ".csv"), deploy_file)
obs_bak    <- sub("\\.csv$", paste0("_backup_", stamp, ".csv"), obs_file)

file.copy(deploy_file, deploy_bak, overwrite = TRUE)
file.copy(obs_file,    obs_bak,    overwrite = TRUE)

# ============================================================
# HELPERS (copied from your "Full script sampling window.txt")
#   - Avoid sourcing the whole file because it contains runnable pipeline code
# ============================================================
parse_ts_safe <- function(x, tz = "UTC") {
  suppressPackageStartupMessages({ library(lubridate); library(stringr) })
  
  if (inherits(x, "POSIXct")) return(x)
  if (inherits(x, "Date"))    return(as.POSIXct(x, tz = tz))
  
  x_chr <- trimws(as.character(x))
  x_chr[x_chr == ""] <- NA_character_
  out <- rep(as.POSIXct(NA, tz = tz), length(x_chr))
  
  # Excel serial dates
  x_num <- suppressWarnings(as.numeric(x_chr))
  looks_serial <- !is.na(x_num) & x_num > 20000 & x_num < 80000 & !str_detect(x_chr, "[-/:A-Za-z]")
  out[looks_serial] <- as.POSIXct(x_num[looks_serial] * 86400, origin = "1899-12-30", tz = tz)
  
  # ISO normalisation
  x_norm <- x_chr
  x_norm <- str_replace(x_norm, "T", " ")
  x_norm <- str_replace(x_norm, "([+-]\\d{2}):(\\d{2})$", "\\1\\2")
  
  # DMY first
  is_dmy <- is.na(out) & !looks_serial & str_detect(x_norm, "^\\d{1,2}[\\./-]\\d{1,2}[\\./-]\\d{2,4}")
  out[is_dmy] <- suppressWarnings(parse_date_time(
    x_norm[is_dmy],
    orders = c("dmy HMSz","dmy HMz","dmy z","dmy HMS","dmy HM","dmy"),
    tz = tz, quiet = TRUE
  ))
  
  # Everything else: YMD
  idx <- is.na(out)
  out[idx] <- suppressWarnings(parse_date_time(
    x_norm[idx],
    orders = c("ymd HMSz","ymd HMz","ymd z","ymd HMS","ymd HM","ymd"),
    tz = tz, quiet = TRUE
  ))
  
  out
}

parse_num_safe <- function(x) {
  if (is.numeric(x)) return(x)
  x <- stringr::str_replace_all(as.character(x), ",", ".")
  x <- stringr::str_replace_all(x, "[^0-9\\.-]", "")
  suppressWarnings(as.numeric(x))
}

# ============================================================
# READ
# ============================================================
dep <- fread(deploy_file)
obs <- fread(obs_file)

stopifnot(all(c("location","lat","lon","ct_model","start","end") %in% names(dep)))
stopifnot(all(c("location","timestamp","species") %in% names(obs)))

# ============================================================
# 1) OBS: rename columns that break autodetection / empty handling
# ============================================================
# dataset_overview_pdf checks for is_empty (not "empty")
if ("empty" %in% names(obs) && !"is_empty" %in% names(obs)) {
  setnames(obs, "empty", "is_empty")
}

# avoid set_camtrap_cols picking camera_trapper as the camera_id column
if ("camera_trapper" %in% names(obs) && !"camera_trapper_flag" %in% names(obs)) {
  setnames(obs, "camera_trapper", "camera_trapper_flag")
}

# ============================================================
# 2) OBS: parse timestamp safely, drop nonsense, then rewrite ISO in SAME column
# ============================================================
obs[, location := str_trim(as.character(location))]

ts <- parse_ts_safe(obs$timestamp, tz = "UTC")

# hard guard against nonsense future/past years (prevents the "2031" mess)
yr <- year(ts)
bad_year <- !is.na(ts) & (yr < 1970 | yr > year(Sys.Date()) + 1)
ts[bad_year] <- as.POSIXct(NA, tz = "UTC")

# if unparsed or bad -> set timestamp to NA (do NOT create extra date columns)
obs[is.na(ts), timestamp := NA_character_]
obs[!is.na(ts), timestamp := format(ts[!is.na(ts)], "%Y-%m-%d %H:%M:%S")]

max_obs_date <- as.Date(max(ts, na.rm = TRUE))
if (!is.finite(max_obs_date)) stop("No parseable observation timestamps after cleaning.")

# ============================================================
# 3) DEP: fill lat/lon for each location (copy/paste within grid cell)
# ============================================================
dep[, location := str_trim(as.character(location))]
dep[, lat := parse_num_safe(lat)]
dep[, lon := parse_num_safe(lon)]

dep[, `:=`(
  lat = ifelse(is.na(lat), mean(lat, na.rm = TRUE), lat),
  lon = ifelse(is.na(lon), mean(lon, na.rm = TRUE), lon)
), by = location]

dep[!is.finite(lat), lat := NA_real_]
dep[!is.finite(lon), lon := NA_real_]

# ============================================================
# 4) DEP: parse start/end, replace "ongoing" with (max obs date + 1 day), rewrite ISO
# ============================================================
start_ts <- parse_ts_safe(dep$start, tz = "UTC")
end_ts   <- parse_ts_safe(dep$end,   tz = "UTC")   # "ongoing" -> NA naturally

# if start/end were date-only, enforce start at 00:00:00, end at 23:59:59
start_has_time <- str_detect(str_trim(as.character(dep$start)), ":")
end_has_time   <- str_detect(str_trim(as.character(dep$end)),   ":")

start_ts[!start_has_time & !is.na(start_ts)] <- as.POSIXct(as.Date(start_ts[!start_has_time & !is.na(start_ts)]), tz = "UTC")
end_ts[!end_has_time   & !is.na(end_ts)]     <- as.POSIXct(as.Date(end_ts[!end_has_time   & !is.na(end_ts)]),   tz = "UTC") +
  hours(23) + minutes(59) + seconds(59)

# replace missing end (ongoing) with last obs date + 1 day (end of that day)
fill_end <- as.POSIXct(max_obs_date + days(1), tz = "UTC") + hours(23) + minutes(59) + seconds(59)
end_ts[is.na(end_ts) & !is.na(start_ts)] <- fill_end

# write ISO back into the SAME columns start/end (so ymd_hms() works downstream)
dep[!is.na(start_ts), start := format(start_ts[!is.na(start_ts)], "%Y-%m-%d %H:%M:%S")]
dep[!is.na(end_ts),   end   := format(end_ts[!is.na(end_ts)],     "%Y-%m-%d %H:%M:%S")]

# ============================================================
# 5) DEP: create deploymentID = location_camera_id
#    Some cameras rotate between grid cells (e.g. C60: cell 12 then cell 23),
#    so camera_id alone is not a unique deployment identifier.
#    location_camera_id is unique per deployment row.
# ============================================================
dep[, deploymentID := paste0(location, "_", camera_id)]
stopifnot(anyDuplicated(dep$deploymentID) == 0L)

# ============================================================
# 6) OBS: assign camera_id and deploymentID by joining on location + timestamp overlap
# ============================================================
obs[, camera_id    := NA_character_]
obs[, deploymentID := NA_character_]
obs[, row_id := .I]

dep_j <- dep[, .(
  location,
  start_ts = parse_ts_safe(start, tz = "UTC"),
  end_ts   = parse_ts_safe(end,   tz = "UTC"),
  camera_id_dep    = as.character(camera_id),
  deploymentID_dep = as.character(deploymentID)
)]
dep_j <- dep_j[!is.na(start_ts) & !is.na(end_ts) & end_ts >= start_ts]

setkey(dep_j, location, start_ts, end_ts)

obs_j <- obs[, .(
  row_id,
  location = location,
  ts = parse_ts_safe(timestamp, tz = "UTC")
)]

m <- dep_j[obs_j,
           on = .(location, start_ts <= ts, end_ts >= ts),
           mult = "last",          # if overlap, take latest start
           nomatch = NA,
           .(row_id = i.row_id,
             camera_id    = camera_id_dep,
             deploymentID = deploymentID_dep)
]

obs[m, on = "row_id", `:=`(camera_id = i.camera_id, deploymentID = i.deploymentID)]

# cleanup helpers
obs[, row_id := NULL]

# put deploymentID first, then camera_id — pipeline picks deploymentID for joins
setcolorder(obs, c("deploymentID", "camera_id",
                    setdiff(names(obs), c("deploymentID", "camera_id"))))
setcolorder(dep, c("deploymentID", setdiff(names(dep), "deploymentID")))

# ============================================================
# WRITE (OVERWRITE)
# ============================================================
fwrite(dep, deploy_file)
fwrite(obs, obs_file)

cat("✅ Backups written:\n  ", deploy_bak, "\n  ", obs_bak, "\n", sep="")
cat("✅ Overwrote:\n  ", deploy_file, "\n  ", obs_file, "\n", sep="")
cat("Deploy rows:", nrow(dep), " | Obs rows:", nrow(obs), "\n")
cat("Unique deploymentIDs in deploy:", uniqueN(dep$deploymentID), "\n")
cat("Obs timestamp NA after cleaning:", sum(is.na(obs$timestamp)), "\n")
cat("Obs camera_id still NA:", sum(is.na(obs$camera_id)), "\n")
cat("Obs deploymentID still NA:", sum(is.na(obs$deploymentID)), "\n")
cat("Obs deploymentIDs not in deploy:",
    sum(!obs$deploymentID %in% dep$deploymentID, na.rm = TRUE), "\n")

# quick sanity checks for your pipeline:
cat("ymd_hms(obs$timestamp) NA count:",
    sum(is.na(lubridate::ymd_hms(obs$timestamp, tz = 'UTC'))), "\n")
cat("ymd_hms(dep$start) NA count:",
    sum(is.na(lubridate::ymd_hms(dep$start, tz = 'UTC'))), "\n")
cat("ymd_hms(dep$end) NA count:",
    sum(is.na(lubridate::ymd_hms(dep$end, tz = 'UTC'))), "\n")


# helpers.R
# ─────────────────────────────────────────────────────────────────────────────
# Shared utility functions for the Sensitivity of Detection project.
# Source this file at the top of any script that needs these helpers.
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(tidyverse))

# ── Tiny helper used internally ───────────────────────────────────────────────
trySuppressWarnings <- function(expr) {
  suppressWarnings(try(expr, silent = TRUE))
}


# ─────────────────────────────────────────────────────────────────────────────
# resolve_dataset_path()
# Find a dataset by name under `root` — returns either a folder path or an
# archive path (.zip / .7z / .tar etc.).
# ─────────────────────────────────────────────────────────────────────────────
resolve_dataset_path <- function(name, root, include_zips = TRUE) {
  exts <- c("zip", "7z", "tar", "gz", "bz2", "rar")
  # 1) exact dir match (case-insensitive)
  dirs    <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dir_hit <- dirs[tolower(basename(dirs)) == tolower(name)]
  if (length(dir_hit) == 1) return(dir_hit)
  if (length(dir_hit) > 1) stop("Multiple directories named '", name, "' under root: ", root)
  # 2) exact archive match
  if (isTRUE(include_zips)) {
    pats    <- paste0("\\.(", paste(exts, collapse = "|"), ")$")
    files   <- list.files(root, pattern = pats, full.names = TRUE, ignore.case = TRUE)
    zip_hit <- files[tolower(tools::file_path_sans_ext(basename(files))) == tolower(name)]
    if (length(zip_hit) == 1) return(zip_hit)
    if (length(zip_hit) > 1) stop("Multiple archives named '", name, "' under root: ", root)
  }
  stop("Dataset not found under root: ", name)
}


# ─────────────────────────────────────────────────────────────────────────────
# base_dataset_name() / ds_name_from_path()
# Thin utilities used across pipeline scripts.
# ─────────────────────────────────────────────────────────────────────────────
base_dataset_name <- function(x) stringr::str_remove(x, "_slice\\d+$")
ds_name_from_path <- function(p) tools::file_path_sans_ext(basename(p))


# ─────────────────────────────────────────────────────────────────────────────
# dataset_output_dir()
# Return the directory that should receive output files for a dataset.
# ─────────────────────────────────────────────────────────────────────────────
dataset_output_dir <- function(archive_path) {
  if (dir.exists(archive_path)) archive_path else dirname(archive_path)
}


# ─────────────────────────────────────────────────────────────────────────────
# read_csv_src()
# Robustly read a CSV from an archive or folder, with character fallback on
# parse problems.
# ─────────────────────────────────────────────────────────────────────────────
read_csv_src <- function(archive_path, pat_end,
                         guess_max = 100000,
                         fallback_to_character = TRUE,
                         verbose = TRUE) {
  re  <- stringr::regex(paste0(pat_end, "$"), ignore_case = TRUE)
  zip <- tools::file_ext(archive_path) %in% c("zip", "7z", "tar", "gz", "bz2", "rar")

  if (zip) {
    lst    <- unzip(archive_path, list = TRUE)$Name
    target <- lst[stringr::str_detect(basename(lst), re)][1]
    if (is.na(target)) stop("No ", pat_end, " found in archive.")
    open_con   <- function() unz(archive_path, target)
    file_label <- basename(target)
  } else {
    files  <- list.files(archive_path, recursive = TRUE, full.names = TRUE)
    target <- files[stringr::str_detect(basename(files), re)][1]
    if (is.na(target)) stop("No ", pat_end, " found in folder.")
    open_con   <- function() target
    file_label <- basename(target)
  }

  strict <- trySuppressWarnings(
    readr::read_csv(open_con(),
                    guess_max      = guess_max,
                    show_col_types = FALSE,
                    progress       = FALSE,
                    locale         = readr::locale(encoding = "UTF-8"))
  )

  if (!inherits(strict, "try-error")) {
    probs <- readr::problems(strict)
    if (nrow(probs) == 0 || !fallback_to_character) return(strict)
    if (verbose) message("Parsing issues in ", file_label,
                         " (", nrow(probs), "); re-reading as all-character.")
  } else if (verbose) {
    message("Strict read failed for ", file_label, "; re-reading as all-character.")
  }

  soft <- readr::read_csv(
    open_con(),
    col_types = readr::cols(.default = readr::col_character()),
    na        = c("", "NA", "N/A"),
    trim_ws   = TRUE,
    progress  = FALSE,
    locale    = readr::locale(encoding = "UTF-8")
  )
  attr(soft, "readr_fallback") <- TRUE
  soft
}


# ─────────────────────────────────────────────────────────────────────────────
# set_camtrap_cols()
# Auto-detect canonical column names from camtrap-style observations and
# deployments tables.
#
# Returns a named list:
#   sci_col, date_col, depl_id_obs, depl_id_deploy,
#   loc_col, lat_col, lon_col, start_col, end_col
#
# Usage at call sites (explicit unpacking):
#   cols <- set_camtrap_cols(obs, deploy)
#   list2env(cols, envir = environment())
# ─────────────────────────────────────────────────────────────────────────────
set_camtrap_cols <- function(obs, deploy) {
  n_obs <- names(obs)
  n_dep <- names(deploy)

  pick_first <- function(pool, exact = NULL, re = NULL) {
    hit <- intersect(exact, pool)
    if (length(hit)) return(hit[1])
    hit <- pool[stringr::str_detect(pool, stringr::regex(re, ignore_case = TRUE))]
    if (length(hit)) return(hit[1])
    NA_character_
  }

  cols <- list(
    sci_col        = pick_first(n_obs, exact = c("scientificName", "species_latin", "species")),
    date_col       = pick_first(n_obs, exact = c("date_recorded", "timestamp", "eventStart")),
    depl_id_obs    = pick_first(n_obs, exact = c("deploymentID", "cameraID"),
                                re = "deploy|deployment|camera|station|device|site|trap"),
    depl_id_deploy = pick_first(n_dep, exact = c("deploymentID"),
                                re = "deploy|deployment|camera|station|device|site|trap"),
    loc_col        = pick_first(n_dep, re = "location|locationID"),
    lat_col        = pick_first(n_dep, exact = c("latitude", "lat", "deploy_lat"), re = "lat"),
    lon_col        = pick_first(n_dep, exact = c("longitude", "long", "lng", "deploy_lon"),
                                re = "lon|lng"),
    start_col      = pick_first(n_dep, re = "start$|^start(date|_time)?$"),
    end_col        = pick_first(n_dep, re = "end$|^end(date|_time)?$")
  )

  missing <- names(cols)[is.na(unlist(cols)) | !nzchar(unlist(cols))]
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

  cols
}


# ─────────────────────────────────────────────────────────────────────────────
# parse_ts_safe()
# Robust timestamp parser: handles POSIXct/Date pass-through, ISO 8601,
# timezone offsets like +01:00, DMY formats, and Excel serial dates.
# ─────────────────────────────────────────────────────────────────────────────
parse_ts_safe <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXct")) return(x)
  if (inherits(x, "Date"))    return(as.POSIXct(x, tz = tz))

  x_chr <- trimws(as.character(x))
  x_chr[x_chr == ""] <- NA_character_
  out <- rep(as.POSIXct(NA, tz = tz), length(x_chr))

  # Excel serial dates (numeric-looking, no letters or separators)
  x_num       <- suppressWarnings(as.numeric(x_chr))
  looks_serial <- !is.na(x_num) & x_num > 20000 & x_num < 80000 &
                  !stringr::str_detect(x_chr, "[-/:A-Za-z]")
  out[looks_serial] <- as.POSIXct(x_num[looks_serial] * 86400,
                                   origin = "1899-12-30", tz = tz)

  # ISO 8601: normalise "T" separator and "+HH:MM" offsets
  x_norm <- x_chr
  x_norm <- stringr::str_replace(x_norm, "T", " ")
  x_norm <- stringr::str_replace(x_norm, "([+-]\\d{2}):(\\d{2})$", "\\1\\2")

  # DMY-first (day first)
  is_dmy <- is.na(out) & !looks_serial &
            stringr::str_detect(x_norm, "^\\d{1,2}[\\./-]\\d{1,2}[\\./-]\\d{2,4}")
  out[is_dmy] <- suppressWarnings(lubridate::parse_date_time(
    x_norm[is_dmy],
    orders = c("dmy HMSz", "dmy HMz", "dmy z", "dmy HMS", "dmy HM", "dmy"),
    tz = tz, quiet = TRUE
  ))

  # Everything else: YMD
  idx      <- is.na(out)
  out[idx] <- suppressWarnings(lubridate::parse_date_time(
    x_norm[idx],
    orders = c("ymd HMSz", "ymd HMz", "ymd z", "ymd HMS", "ymd HM", "ymd"),
    tz = tz, quiet = TRUE
  ))

  out
}


# ─────────────────────────────────────────────────────────────────────────────
# parse_num_safe()
# NA-tolerant numeric coercion with comma→dot normalisation.
# ─────────────────────────────────────────────────────────────────────────────
parse_num_safe <- function(x) {
  if (is.numeric(x)) return(x)
  x <- stringr::str_replace_all(x, ",", ".")
  x <- stringr::str_replace_all(x, "[^0-9\\.-]", "")
  suppressWarnings(as.numeric(x))
}


# ─────────────────────────────────────────────────────────────────────────────
# pad_bbox_safe()
# Pad an sf bbox by fraction f; ensures a minimum extent of min_deg in each
# direction to avoid zero-width bboxes.
# ─────────────────────────────────────────────────────────────────────────────
pad_bbox_safe <- function(bb, f = 0.12, min_deg = 0.02) {
  dx <- as.numeric(bb["xmax"] - bb["xmin"])
  dy <- as.numeric(bb["ymax"] - bb["ymin"])
  if (dx == 0) { bb["xmin"] <- bb["xmin"] - min_deg / 2; bb["xmax"] <- bb["xmax"] + min_deg / 2 }
  if (dy == 0) { bb["ymin"] <- bb["ymin"] - min_deg / 2; bb["ymax"] <- bb["ymax"] + min_deg / 2 }
  dx <- as.numeric(bb["xmax"] - bb["xmin"])
  dy <- as.numeric(bb["ymax"] - bb["ymin"])
  bb["xmin"] <- bb["xmin"] - f * dx
  bb["xmax"] <- bb["xmax"] + f * dx
  bb["ymin"] <- bb["ymin"] - f * dy
  bb["ymax"] <- bb["ymax"] + f * dy
  bb
}


# ─────────────────────────────────────────────────────────────────────────────
# daily_deployment_effort_fast()
# Count active deployments per day using a vectorised cumsum approach.
# O(n_deployments + n_days) — much faster than the map2+unnest version for
# large or long-running deployments.
# ─────────────────────────────────────────────────────────────────────────────
daily_deployment_effort_fast <- function(deploy, start_col, end_col,
                                         drop_leap_day = FALSE) {
  dep_tbl <- deploy |>
    transmute(
      start = parse_ts_safe(.data[[start_col]]),
      end   = parse_ts_safe(.data[[end_col]])
    ) |>
    filter(!is.na(start), !is.na(end), end > start) |>
    mutate(start_d = lubridate::as_date(start),
           end_d   = lubridate::as_date(end))

  if (!nrow(dep_tbl))
    return(tibble(date = as.Date(character()), n_active_deployments = integer()))

  # +1 on start day, -1 on day after end
  events <- c(
    stats::setNames(rep(+1L, nrow(dep_tbl)), as.character(dep_tbl$start_d)),
    stats::setNames(rep(-1L, nrow(dep_tbl)), as.character(dep_tbl$end_d + 1L))
  )
  event_df <- tibble(
    date  = as.Date(names(events)),
    delta = unname(events)
  ) |>
    summarise(delta = sum(delta), .by = date) |>
    arrange(date)

  full_dates <- tibble(date = seq(min(dep_tbl$start_d), max(dep_tbl$end_d), by = "day"))
  if (drop_leap_day)
    full_dates <- full_dates |>
      filter(!(lubridate::month(date) == 2L & lubridate::mday(date) == 29L))

  full_dates |>
    left_join(event_df, by = "date") |>
    mutate(
      delta                = replace_na(delta, 0L),
      n_active_deployments = as.integer(cumsum(delta))
    ) |>
    select(date, n_active_deployments)
}


# ─────────────────────────────────────────────────────────────────────────────
# daily_deployment_effort()
# Original (slower) map2+unnest version. Kept for backward compatibility.
# Prefer daily_deployment_effort_fast() for new code.
# ─────────────────────────────────────────────────────────────────────────────
daily_deployment_effort <- function(deploy, start_col, end_col,
                                    drop_leap_day = FALSE) {
  dep_tbl <- deploy |>
    transmute(
      start = parse_ts_safe(.data[[start_col]]),
      end   = parse_ts_safe(.data[[end_col]])
    ) |>
    filter(!is.na(start), !is.na(end), end > start)

  if (!nrow(dep_tbl))
    return(tibble(date = as.Date(character()), n_active_deployments = integer()))

  start_date <- lubridate::as_date(min(dep_tbl$start))
  end_date   <- lubridate::as_date(max(dep_tbl$end))
  full_dates <- tibble(date = seq(start_date, end_date, by = "day"))
  if (drop_leap_day)
    full_dates <- full_dates |>
      filter(!(lubridate::month(date) == 2 & lubridate::mday(date) == 29))

  dep_tbl |>
    mutate(date_seq = purrr::map2(lubridate::as_date(start),
                                  lubridate::as_date(end),
                                  \(s, e) seq(s, e, by = "day"))) |>
    select(date_seq) |>
    tidyr::unnest(date_seq) |>
    rename(date = date_seq) |>
    (\(d) if (drop_leap_day)
       filter(d, !(lubridate::month(date) == 2 & lubridate::mday(date) == 29))
     else d)() |>
    count(date, name = "n_active_deployments") |>
    right_join(full_dates, by = "date") |>
    mutate(n_active_deployments = replace_na(n_active_deployments, 0L)) |>
    arrange(date)
}


# ─────────────────────────────────────────────────────────────────────────────
# monthly_deployment_effort()
# Aggregate daily effort (from daily_deployment_effort*) to months.
# ─────────────────────────────────────────────────────────────────────────────
monthly_deployment_effort <- function(daily_df) {
  daily_df |>
    mutate(month_date = lubridate::floor_date(date, "month")) |>
    summarise(trap_nights = sum(n_active_deployments), .by = month_date) |>
    arrange(month_date)
}


# ─────────────────────────────────────────────────────────────────────────────
# spp_keep
# European terrestrial mammals accepted by the detection-sensitivity pipeline.
# ─────────────────────────────────────────────────────────────────────────────
spp_keep <- c(
  # ── Lagomorpha
  "Oryctolagus cuniculus",
  "Lepus europaeus", "Lepus timidus", "Lepus granatensis",
  "Lepus castroviejoi", "Lepus corsicanus", "Lepus tolai",
  # ── Eulipotyphla / Rodentia
  "Erinaceus europaeus", "Erinaceus roumanicus",
  "Talpa europaea", "Talpa occidentalis", "Talpa romana", "Talpa caeca",
  "Sciurus vulgaris",
  # ── Primates
  "Macaca sylvanus",
  # ── Carnivora
  "Canis lupus", "Canis aureus", "Vulpes vulpes", "Nyctereutes procyonoides",
  "Ursus arctos",
  "Lynx lynx", "Felis silvestris", "Lynx pardinus",
  "Meles meles",
  "Martes martes", "Martes foina",
  "Mustela nivalis", "Mustela erminea", "Mustela putorius", "Mustela lutreola",
  "Neogale vison", "Mustela vison",
  "Vormela peregusna",
  "Gulo gulo",
  "Lutra lutra",
  "Genetta genetta",
  "Herpestes ichneumon",
  "Procyon lotor",
  "Myocastor coypus",
  "Castor fiber",
  # ── Artiodactyla
  "Sus scrofa",
  "Alces alces",
  "Capreolus capreolus",
  "Cervus elaphus", "Cervus nippon",
  "Dama dama",
  "Rangifer tarandus",
  "Bison bonasus",
  "Rupicapra rupicapra", "Rupicapra pyrenaica",
  "Capra ibex", "Capra pyrenaica", "Capra aegagrus",
  "Ovis orientalis", "Ovis gmelini", "Ovis musimon",
  "Bos taurus", "Ovis aries", "Capra hircus",
  "Ammotragus lervia", "Muntiacus reevesi", "Hydropotes inermis"
)

# ─────────────────────────────────────────────────────────────────────────────
# protocol_windows()
# Canonical definitions of named protocol windows:
#   SNAP_EU_CORE   — Snapshot Europe core (Sep 1 – Oct 31, 61d)
#   SNAP_EU_BUFFER — Snapshot Europe buffer (Aug 18 – Nov 14, 89d)
#   EOW_EARLY      — EOW/ENETWILD early half (60d ending before CORE center)
#   EOW_LATE       — EOW/ENETWILD late half  (60d starting at CORE center)
# Returns a tibble with DOY-based columns matching the window_grid format,
# plus start_md / end_md for calendar-date conversion.
# ─────────────────────────────────────────────────────────────────────────────
protocol_windows <- function() {
  doy <- function(m, d) lubridate::yday(as.Date(sprintf("2001-%02d-%02d", m, d)))

  core_start <- doy(9, 1)    # Sep 1  = 244
  core_end   <- doy(10, 31)  # Oct 31 = 304
  core_len   <- core_end - core_start + 1L  # 61

  buf_start <- doy(8, 18)    # Aug 18 = 230
  buf_end   <- doy(11, 14)   # Nov 14 = 318
  buf_len   <- buf_end - buf_start + 1L     # 89

  # EOW/ENETWILD: two 60-day windows around the CORE center (DOY 274 = Oct 1)
  core_mid    <- core_start + floor(core_len / 2)  # DOY 274
  eow_len     <- 60L
  eow_e_start <- core_mid - eow_len                # DOY 214 (Aug 2)
  eow_e_end   <- core_mid - 1L                     # DOY 273 (Sep 30)
  eow_l_start <- core_mid                          # DOY 274 (Oct 1)
  eow_l_end   <- core_mid + eow_len - 1L           # DOY 333 (Nov 29)

  ref  <- as.Date("2001-01-01")
  ids  <- c("SNAP_EU_CORE", "SNAP_EU_BUFFER", "EOW_EARLY", "EOW_LATE")
  doys <- c(core_start,  buf_start,  eow_e_start, eow_l_start)
  doye <- c(core_end,    buf_end,    eow_e_end,   eow_l_end)
  lens <- c(core_len,    buf_len,    eow_len,     eow_len)

  tibble::tibble(
    window_id = ids,
    start_doy = doys,
    length_d  = lens,
    end_doy   = doye,
    start_md  = format(ref + doys - 1L, "%m-%d"),
    end_md    = format(ref + doye - 1L, "%m-%d")
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# protocol_calendar_dates()
# Convert protocol windows to calendar dates for a given anchor period.
# Returns tibble: window, win_start, win_end (Date).
# ─────────────────────────────────────────────────────────────────────────────
protocol_calendar_dates <- function(anchor_start, anchor_end) {
  pw    <- protocol_windows()
  years <- unique(c(lubridate::year(anchor_start), lubridate::year(anchor_end)))

  purrr::map_dfr(years, \(y) {
    pw |>
      dplyr::transmute(
        window    = window_id,
        win_start = as.Date(paste0(y, "-", start_md)),
        win_end   = as.Date(paste0(y, "-", end_md))
      )
  }) |>
    dplyr::filter(win_end >= anchor_start, win_start <= anchor_end) |>
    dplyr::distinct()
}

# ─────────────────────────────────────────────────────────────────────────────
# make_window_template()
# Build a shared calendar of rolling windows by start DOY and length.
# Appends protocol windows (SNAP_EU_CORE, SNAP_EU_BUFFER, EOW_EARLY, EOW_LATE) by default.
# Writes `window_grid` to .GlobalEnv and returns it invisibly.
# ─────────────────────────────────────────────────────────────────────────────
make_window_template <- function(step_doy = 7,
                                 lengths_d = seq(15, 120, by = 7),
                                 include_protocols = TRUE) {
  start_seq  <- seq(1, 365, by = step_doy)
  window_tbl <- tidyr::crossing(start_doy = start_seq, length_d = lengths_d) |>
    mutate(
      end_doy   = start_doy + length_d - 1L,
      window_id = sprintf("d%03d_L%d", start_doy, length_d),
      .before   = 1
    )

  window_tbl <- dplyr::bind_rows(
    tibble(window_id = "FULL", start_doy = 1L, length_d = 365L, end_doy = 365L),
    window_tbl
  )

  if (include_protocols) {
    pw <- protocol_windows() |> dplyr::select(window_id, start_doy, length_d, end_doy)
    window_tbl <- dplyr::bind_rows(window_tbl, pw) |>
      dplyr::distinct(window_id, .keep_all = TRUE)
  }

  assign("window_grid", window_tbl, envir = .GlobalEnv)
  invisible(window_tbl)
}


make_window_template(step_doy = 7, lengths_d = seq(15, 120, by = 7))


# ====================================================================
# find_anchors()
# ====================================================================
# Identify optimal 12-month "anchor" windows in a dataset by scoring
# rolling 365-day windows on effort, detection density, and gap quality.
# Returns a tibble with one row per kept slice.
#
# Dependencies:  read_csv_src(), daily_deployment_effort(), parse_ts_safe()
#                (all defined above in helpers.R)
# External pkgs: zoo (for rolling window calculations)
# ====================================================================

find_anchors <- function(
    archive_path,
    max_slices            = 20,
    min_gap_days          = 365,
    min_effort_per_slice  = 2000,
    min_balance_ratio     = 0.20,
    drop_feb29            = TRUE,
    keep_species      = spp_keep,
    independence_mins = 30,
    # allow searching slightly outside obs span (kept species)
    obs_buffer_days = 30,
    # quality controls (defaults)
    min_detect_days = 250,  # drop slices with fewer detect-days (truncated anchors)
    max_zero_frac   = 0.50,
    max_zero_gap    = 60,
    w_zero_frac     = 0.6,
    w_long_gap      = 0.8,
    
    # minimal logging
    verbose = TRUE
) {
  suppressPackageStartupMessages({ library(tidyverse); library(lubridate); library(zoo); library(stringr) })
  vmsg <- function(...) if (isTRUE(verbose)) message(...)
  
  ds <- tools::file_path_sans_ext(basename(archive_path))
  vmsg("\n--- ", ds, " ---")
  
  deploy <- read_csv_src(archive_path, "deployments\\.csv")
  obs <- tryCatch(read_csv_src(archive_path, "observations\\.csv"),
                  error = function(e) NULL)
  
  # ---- minimal column detection ----
  pick_first <- function(pool, exact = NULL, re = NULL) {
    if (is.null(exact)) exact <- character(0)
    hit <- intersect(exact, pool); if (length(hit)) return(hit[1])
    if (!is.null(re)) {
      hit <- pool[str_detect(pool, regex(re, ignore_case = TRUE))]
      if (length(hit)) return(hit[1])
    }
    NA_character_
  }
  
  # deployments: start/end
  n_dep <- names(deploy)
  start_col <- pick_first(
    n_dep,
    exact = c("deploymentStart","start","start_date","startDate","start_datetime"),
    re    = "deployment\\s*start|deploy.*start|(^|_)start($|_)|start(date|_time)?$"
  )
  end_col <- pick_first(
    n_dep,
    exact = c("deploymentEnd","end","end_date","endDate","end_datetime"),
    re    = "deployment\\s*end|deploy.*end|(^|_)end($|_)|end(date|_time)?$"
  )
  if (!nzchar(start_col) || !nzchar(end_col)) {
    stop("Could not detect deployment start/end columns in deployments table for: ", ds)
  }
  
  # ---- empty-result template (shared columns for early returns) ----
  empty_result <- function() {
    out <- tibble(dataset=character(), slice=character(),
                  anchor_start=as.Date(character()), anchor_end=as.Date(character()),
                  trap_days_12m=numeric(), events_12m=numeric(), detect_days_12m=numeric(),
                  p_E=numeric(), zero_frac_12m=numeric(), max_zero_gap_12m=integer(),
                  score=numeric())
    attr(out, "gates_relaxed") <- FALSE
    out
  }
  
  # ---- daily effort (trap-days) over FULL deployment span ----
  daily_effort <- daily_deployment_effort(
    deploy, start_col, end_col, drop_leap_day = drop_feb29
  ) %>% transmute(date, trap_days = as.numeric(n_active_deployments))
  
  if (!nrow(daily_effort) || sum(daily_effort$trap_days, na.rm = TRUE) == 0) {
    vmsg("No valid effort dates (empty/zero). Returning empty.")
    return(empty_result())
  }
  
  dep_min <- min(daily_effort$date)
  dep_max <- max(daily_effort$date)
  vmsg("Deployment span (effort-derived): ", dep_min, " -> ", dep_max)
  
  # ---- daily detections (independent events; defines OBS SPAN using kept species) ----
  daily_detects <- tibble(date = as.Date(character()), detects = integer())
  obs_min <- NA_Date_
  obs_max <- NA_Date_
  
  if (!is.null(obs)) {
    n_obs <- names(obs)
    sci_col  <- pick_first(n_obs, exact = c("scientificName","species_latin","species"),
                           re = "scientificname|species")
    date_col <- pick_first(n_obs, exact = c("date_recorded","timestamp","eventStart"),
                           re = "date|time|timestamp|eventstart")
    cam_col  <- pick_first(n_obs, re = "camera|station|device|deploy|deployment|site|trap")
    
    if (nzchar(sci_col) && nzchar(date_col)) {
      obs_f <- obs %>%
        transmute(
          species   = str_trim(as.character(.data[[sci_col]])),
          ts        = parse_ts_safe(.data[[date_col]]),
          camera_id = if (nzchar(cam_col)) as.character(.data[[cam_col]]) else "__all__"
        ) %>%
        filter(!is.na(ts), !is.na(species), nzchar(species)) %>%
        { if (!is.null(keep_species)) filter(., species %in% keep_species) else . } %>%
        mutate(date = as_date(ts)) %>%
        { if (drop_feb29) filter(., !(month(date)==2 & mday(date)==29)) else . }
      
      if (nrow(obs_f)) {
        obs_min <- min(obs_f$date, na.rm = TRUE)
        obs_max <- max(obs_f$date, na.rm = TRUE)
        
        # ---- collapse to independent events (camera x species x time gap) ----
        if (independence_mins > 0 && nzchar(cam_col)) {
          obs_f <- obs_f %>%
            arrange(camera_id, species, ts) %>%
            group_by(camera_id, species) %>%
            mutate(
              gap      = ts - lag(ts),
              new_event = is.na(gap) | gap > minutes(independence_mins)
            ) %>%
            filter(new_event) %>%
            ungroup() %>%
            select(-gap, -new_event)
          
          vmsg("Collapsed to independent events (", independence_mins,
               " min gap, camera x species): ", nrow(obs_f), " events.")
        } else if (!nzchar(cam_col)) {
          vmsg("WARNING: No camera/deployment ID column found in observations. ",
               "Cannot collapse to independent events; using raw observation counts.")
        }
        
        daily_detects <- obs_f %>%
          count(date, name = "detects") %>%
          mutate(detects = as.integer(detects))
        
        vmsg("Obs span (kept spp): ", obs_min, " -> ", obs_max,
             " | total events=", sum(daily_detects$detects, na.rm = TRUE),
             " | det-days=", nrow(daily_detects))
      } else {
        vmsg("Obs span unavailable after keep_species filter (no rows). Using full deployment span.")
      }
    } else {
      vmsg("Obs columns not detected (sci/date). Using full deployment span.")
    }
  } else {
    vmsg("observations.csv not readable. Using full deployment span.")
  }
  
  # ---- define search span (deployment span, clamped to buffered obs span if available) ----
  span_start <- dep_min
  span_end   <- dep_max
  
  if (!is.na(obs_min) && !is.na(obs_max)) {
    buf <- as.integer(obs_buffer_days)
    span_start <- max(span_start, obs_min - days(buf))
    span_end   <- min(span_end,   obs_max + days(buf))
    vmsg("Search span (obs +/-", buf, "d, clamped): ", span_start, " -> ", span_end)
  } else {
    vmsg("Search span: full deployment span (no obs restriction).")
  }
  
  # calendar
  full_days <- tibble(date = seq(span_start, span_end, by = "day"))
  if (drop_feb29) full_days <- full_days %>% filter(!(month(date)==2 & mday(date)==29))
  
  if (nrow(full_days) < 365) {
    vmsg("Span < 365 days after restriction/Feb29. Returning empty.")
    return(empty_result())
  }
  
  daily <- full_days %>%
    left_join(daily_effort,  by = "date") %>%
    left_join(daily_detects, by = "date") %>%
    mutate(trap_days = replace_na(trap_days, 0),
           detects   = replace_na(detects,   0L)) %>%
    arrange(date)
  
  # ---- K_max explanation (this is the key "why only K") ----
  eff_nonzero <- daily %>% filter(trap_days > 0)
  if (nrow(eff_nonzero)) {
    eff_span_days <- as.integer(max(eff_nonzero$date) - min(eff_nonzero$date)) + 1L
    K_max_raw <- floor(eff_span_days / 365)
    K_max <- max(1L, min(max_slices, K_max_raw))
    vmsg("Effort-active span within search window: ", eff_span_days, " days",
         " -> floor(span/365)=", K_max_raw,
         " -> K_max=", K_max)
    if (K_max == 1L) {
      vmsg("Why only K=1: you have < 730 effort-active days in the search span, so >=2 full 365d slices cannot fit.")
      vmsg("Tip: increase obs_buffer_days or relax obs-span restriction if you truly have >1 year of usable observations.")
    } else if (K_max_raw > max_slices) {
      vmsg("Note: K is capped by max_slices=", max_slices, " (raw would be ", K_max_raw, ").")
    }
  } else {
    K_max <- 1L
    vmsg("No nonzero effort days in search span -> K_max=1.")
  }
  
  # ---- rolling window metrics (365-row windows) ----
  roll365_sum  <- function(x) zoo::rollsum(x, 365, align = "left", fill = NA)
  roll365_mean <- function(x) zoo::rollmean(x, 365, align = "left", fill = NA)
  roll365_max_gap <- function(is_zero) {
    zoo::rollapply(is_zero, 365, align = "left", fill = NA, FUN = function(v) {
      r <- rle(v)
      if (any(r$values, na.rm = TRUE)) max(r$lengths[r$values], na.rm = TRUE) else 0L
    })
  }
  
  is_zero <- daily$trap_days == 0
  
  d2 <- daily %>%
    mutate(
      anchor_end = dplyr::lead(date, 364),
      E          = roll365_sum(trap_days),
      D          = roll365_sum(detects),                    # independent events
      Dd         = roll365_sum(as.integer(detects > 0)),    # tie-break: detection-days
      Z          = roll365_mean(as.numeric(is_zero)),
      G          = roll365_max_gap(is_zero)
    ) %>%
    filter(!is.na(E), !is.na(Z), !is.na(G), !is.na(anchor_end))
  
  totalE <- sum(daily$trap_days, na.rm = TRUE)
  
  p95E <- suppressWarnings(as.numeric(quantile(d2$E, probs = 0.95, na.rm = TRUE, type = 7)))
  if (!is.finite(p95E) || p95E <= 0) p95E <- max(d2$E, na.rm = TRUE)
  
  d2 <- d2 %>%
    mutate(
      E_norm = pmin(1, if (p95E > 0) E / p95E else 0),
      G_norm = pmin(1, G / 365),
      score  = E_norm - w_zero_frac * Z - w_long_gap * G_norm
    )
  
  # contiguity filter (track whether relaxed)
  gates_relaxed <- FALSE
  n_before <- nrow(d2)
  d2_ok <- d2 %>% filter(Z <= max_zero_frac, G <= max_zero_gap)
  vmsg("Contiguity filter kept ", nrow(d2_ok), "/", n_before, " anchors",
       " (Z<=", max_zero_frac, ", G<=", max_zero_gap, ").")
  if (!nrow(d2_ok)) {
    vmsg("WARNING: No anchors met contiguity thresholds -> relaxing (using all anchors).")
    d2_ok <- d2
    gates_relaxed <- TRUE
  }
  
  # ---- selection with tie-break: score -> detection-days -> midpoint distance ----
  mid_date <- span_start + floor(as.numeric(span_end - span_start) / 2)
  
  choose_windows <- function(df, k, gap_days, mid_date) {
    picked <- tibble()
    pool <- df %>%
      arrange(desc(score), desc(Dd), abs(as.numeric(date - mid_date)))
    
    while (nrow(picked) < k && nrow(pool) > 0) {
      top <- slice(pool, 1)
      picked <- bind_rows(picked, top)
      pool <- pool %>%
        filter(abs(as.numeric(difftime(date, top$date, units = "days"))) > gap_days)
    }
    picked
  }
  
  vmsg("Trying K = 1..", K_max, " (min_gap_days=", min_gap_days, ").")
  
  last_keep <- choose_windows(d2_ok, 1, min_gap_days, mid_date)
  
  for (K in seq_len(K_max)) {
    cand <- choose_windows(d2_ok, K, min_gap_days, mid_date)
    if (!nrow(cand)) break
    
    if (K == 1L) {
      pass <- TRUE
      reason <- "K=1 always accepted."
    } else {
      E_min <- min(cand$E)
      E_max <- max(cand$E)
      
      # Gate (a): absolute effort floor
      pass_floor <- E_min >= min_effort_per_slice
      # Gate (b): balance — weakest slice vs strongest
      balance     <- if (E_max > 0) E_min / E_max else 1
      pass_balance <- balance >= min_balance_ratio
      
      pass <- pass_floor && pass_balance
      
      if (pass) {
        total_selected_E <- sum(cand$E, na.rm = TRUE)
        coverage_pct <- round(100 * total_selected_E / totalE, 1)
        reason <- paste0("PASS: min E=", round(E_min), 
                         ", balance=", round(balance, 2),
                         " (coverage: ", coverage_pct, "% of all effort)")
      } else {
        parts <- character(0)
        if (!pass_floor)
          parts <- c(parts, paste0("min E=", round(E_min),
                                   " < floor ", min_effort_per_slice))
        if (!pass_balance)
          parts <- c(parts, paste0("balance=", round(balance, 2),
                                   " < ", min_balance_ratio))
        reason <- paste0("FAIL: ", paste(parts, collapse = "; "))
      }
    }
    
    vmsg("K=", K, ": picked ", nrow(cand), " | ", reason)
    if (pass) last_keep <- cand
  }
  
  vmsg("Final: K=", nrow(last_keep),
       " | anchors: ", paste(as.character(last_keep$date), collapse = ", "))
  
  result <- last_keep %>%
    transmute(
      dataset          = ds,
      anchor_start     = date,
      anchor_end       = anchor_end,
      trap_days_12m    = E,
      events_12m       = as.numeric(D),
      detect_days_12m  = as.numeric(Dd),
      p_E              = if (totalE > 0) E / totalE else NA_real_,
      zero_frac_12m    = Z,
      max_zero_gap_12m = as.integer(G),
      score            = score
    ) %>%
    arrange(anchor_start) %>%
    mutate(slice = as.character(row_number()), .after = dataset)
  
  # ---- drop truncated slices (< min_detect_days of actual observation days) ----
  if (min_detect_days > 0 && nrow(result)) {
    n_before <- nrow(result)
    result <- result %>% filter(detect_days_12m >= min_detect_days)
    n_dropped <- n_before - nrow(result)
    if (n_dropped > 0) {
      vmsg("Dropped ", n_dropped, " slice(s) with detect_days < ", min_detect_days)
      # Re-number remaining slices sequentially
      result <- result %>% mutate(slice = as.character(row_number()))
    }
  }
  
  attr(result, "gates_relaxed") <- gates_relaxed
  result
}

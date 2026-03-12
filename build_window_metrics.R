# --------------------------------------------------------------------
# build_window_metrics()
# --------------------------------------------------------------------
# Purpose
#   • Takes a *single* dataset (deployments + observations already
#     cropped to its 12-month anchor) and a shared window template
#     (`window_tbl`), then produces window-level detection metrics.
#   • Internal workflow:
#       1. Collapse raw observations to independent events
#          (≥ `independence_mins` gap).
#       2. Intersect every camera deployment with each window to
#          compute effort (whole trap-days).
#       3. Join events to windows; build camera × occasion detection
#          matrices at the chosen `occ_grain_days` using breaks
#          `start_date + 1  …  end_date − 1`.
#       4. Apply post-history thresholds  
#          (≥ `min_events` **and** ≥ `min_occasions_pos`
#           **and** ≥ `min_sites_pos`).
#       5. Fit three metrics for every retained species × window:  
#            – Occupancy p̂ (unmarked, intercept-only)  
#            – Time-to-event λ → P (right-censored exponential)  
#            – Naïve rate R → P (events / trap-days)  
#          Standard errors for P’s via the delta method.
#       6. Assemble tidy outputs.
#
# Required inputs  (all passed as arguments)
#   • window_grid    – tibble: window_id, start_date, end_date
#   • deploy_tbl     – camera_id, start, end   (UTC POSIXct)
#   • obs_tbl        – camera_id, species, timestamp (UTC POSIXct),
#                      already bird/pet-filtered & cropped to anchor
# --------------------------------------------------------------------


build_window_metrics <- function(
    window_grid, 
    deploy_tbl, 
    obs_tbl,
    a_start,            
    a_end,              
    independence_mins = 30, 
    occ_grain_days = 1,
    min_events = 20, 
    min_occasions_pos = 5, 
    min_sites_pos = 5,
    trap_array_km = NA_real_,          # NEW
    latitude_deg  = NA_real_,          # NEW
    drop_leap_day = TRUE)
{
  
  suppressPackageStartupMessages({
    library(tidyverse); library(lubridate); library(fuzzyjoin)
    library(unmarked);  library(zoo); library(geosphere)
  })
  
  ## ───── 1. Collapse raw detections to independent events ────────────
  events <- obs_tbl %>% 
    mutate(timestamp = with_tz(timestamp, tzone = "UTC")) %>% 
    { if (drop_leap_day) filter(., !(month(timestamp)==2 & mday(timestamp)==29)) else . } %>% 
    arrange(camera_id, species, timestamp) %>% 
    group_by(camera_id, species) %>% 
    mutate(gap = timestamp - lag(timestamp),
           new_event = is.na(gap) | gap > minutes(independence_mins)) %>% 
    filter(new_event) %>% 
    ungroup()
  
  ## ───── 2. Camera × window effort (whole days) ──────────────────────
  intersect_windows <- function(cam_start_virtual, cam_end_virtual) {
    window_grid %>%
      rowwise() %>%
      mutate(effort_days = {
        if (end_virtual > 365) {
          # Wrap-around case
          seg1 <- max(0, pmin(cam_end_virtual, 365) - pmax(cam_start_virtual, start_virtual) + 1L)
          seg2 <- max(0, pmin(cam_end_virtual, end_virtual - 365) - 0 + 1L)
          seg1 + seg2
        } else {
          # Normal case
          max(0, pmin(cam_end_virtual, end_virtual) - pmax(cam_start_virtual, start_virtual) + 1L)
        }
      }) %>%
      ungroup() %>%
      filter(effort_days > 0) %>%
      select(window_id, effort_days)
  }
  
  window_camera <- deploy_tbl %>%
    mutate(tmp = purrr::map2(start_virtual, end_virtual, intersect_windows)) %>%
    unnest(tmp) %>%
    select(camera_id, window_id, effort_days)
  
  
  ## ───── 3. Assign events to windows (non-equi join) ─────────────────
  ivl_tbl <- window_grid %>%
    transmute(
      window_id,
      start_virtual,
      end_virtual
    )
  
  window_events_raw <- events %>%
    fuzzyjoin::interval_left_join(
      ivl_tbl,
      by = c("doy_virtual" = "start_virtual",
             "doy_virtual" = "end_virtual")) %>%
    filter(!is.na(window_id))

  # keep the summary version only for quick tallies
  window_events_sum <- window_events_raw %>%       
    group_by(window_id, camera_id, species) %>%    
    summarise(n_events = n(), .groups = "drop")
  
  ## ───── 4. Build detection matrices & compute metrics ───────────────
  results_list <- vector("list", length = nrow(window_grid))
  names(results_list) <- window_grid$window_id
  
  for (i in seq_len(nrow(window_grid))) {
    win <- window_grid[i, ]
    
    ## 4a. occasion breaks -------------------------------------------------------
    if (win$end_date < win$start_date) {          # wrap-around window
      # segment 1: last part of anchor year
      seg1 <- seq(
        as.POSIXct(win$start_date, tz = "UTC"),
        as.POSIXct(a_end          , tz = "UTC") + days(1),
        by = paste0(occ_grain_days, " days")
      )
      
      # segment 2: first part of next anchor year
      seg2 <- seq(
        as.POSIXct(a_start        , tz = "UTC"),
        as.POSIXct(win$end_date   , tz = "UTC") + days(1),
        by = paste0(occ_grain_days, " days")
      )
      
      occ_breaks <- sort(unique(c(seg1, seg2)))   # POSIXct, increasing, no dups
      L_days     <- win$length_d                  # already stored in template
      
    } else {                                      # normal window
      occ_breaks <- seq(
        as.POSIXct(win$start_date, tz = "UTC"),
        as.POSIXct(win$end_date  , tz = "UTC") + days(1),
        by = paste0(occ_grain_days, " days")
      )
      L_days <- win$length_d
    }
    
    n_occ <- L_days               # used later for labelling
    
    
    ## 4b. subset effort & events for this window
    cam_eff <- window_camera %>% 
      filter(window_id == win$window_id)
    ev_sub  <- window_events_raw %>% filter(window_id == win$window_id)
    
    ## 4c. create list of species ↦ matrices
    species_vec <- unique(ev_sub$species)
    species_out <- vector("list", length(species_vec))
    names(species_out) <- species_vec
    
    for (sp in species_vec) {
      ev_sp <- ev_sub %>% filter(species == sp)
      
      # empty matrix of zeros
      mat <- matrix(0,
                    nrow = nrow(cam_eff),
                    ncol = L_days,
                    dimnames = list(cam_eff$camera_id,
                                    paste0("oc", seq_len(L_days))))
      
      # Fill detections
      if (nrow(ev_sp) > 0) {
        occ_ids <- cut(ev_sp$timestamp,
                       breaks = occ_breaks,
                       right = FALSE,
                       labels = FALSE)
        idx <- cbind(match(ev_sp$camera_id, cam_eff$camera_id), occ_ids)
        mat[idx] <- 1
      }
      
      ## 4d. compute summary counts for thresholds
      n_events_total   <- sum(mat)
      n_occ_pos        <- sum(colSums(mat) > 0)
      n_sites_pos      <- sum(rowSums(mat) > 0)
      
      if (n_events_total >= min_events &&
          n_occ_pos      >= min_occasions_pos &&
          n_sites_pos    >= min_sites_pos) {
        
        # ----- occupancy model
        umf  <- unmarkedFrameOccu(y = mat)
        
        all0 <- sum(mat) == 0
        all1 <- sum(mat) == length(mat)
        
        if (all0 || all1) {                        # degenerate matrix
          p_hat  <- if (all1) 1 else 0
          p_se   <- NA_real_
          psi_hat <- if (all1) 1 else 0            # occupancy mirrors detections
          psi_se  <- NA_real_
          conv    <- FALSE
        } else {
          fit <- tryCatch(occu(~1 ~ 1, data = umf), error = function(e) NULL)
          
          if (!is.null(fit) && fit@opt$convergence == 0) {
            det_est  <- predict(fit, "det")        # detection component
            state_est<- predict(fit, "state")      # occupancy component
            
            p_hat    <- det_est$Predicted[1]
            p_se     <- det_est$SE[1]
            psi_hat  <- state_est$Predicted[1]
            psi_se   <- state_est$SE[1]
            conv     <- TRUE
          } else {                                 # failed Hessian etc.
            fit2 <- tryCatch(occu(~1 ~ 1, data = umf, se = FALSE), error = function(e) NULL)
            p_hat   <- if (!is.null(fit2)) plogis(coef(fit2)["p"])           else NA_real_
            psi_hat <- if (!is.null(fit2)) plogis(coef(fit2)["(Intercept)"]) else NA_real_
            p_se  <- psi_se <- NA_real_
            conv  <- FALSE
          }
        }
        
        # ----- TTE (λ)
        first_time <- apply(mat, 1, function(x) {
          if (any(x == 1)) {
            occ_first <- which(x == 1)[1]
            (occ_first - 1) * occ_grain_days + 0.001   # tiny offset avoids zero exposure
          } else NA_real_
        })
        censored <- is.na(first_time)
        n_events_cam <- sum(!censored)
        exposure <- sum(first_time[!censored]) + sum(censored) * L_days
        
        lambda <- if (exposure > 0) n_events_cam / exposure else 0
        se_lambda <- if (exposure > 0 && n_events_cam > 0)
          sqrt(n_events_cam) / exposure else NA_real_
        p_tte  <- 1 - exp(-lambda * L_days)
        se_p_tte <- if (!is.na(se_lambda))
          exp(-lambda * L_days) * L_days * se_lambda else NA_real_
        
        # ----- Naïve rate
        trap_days <- sum(cam_eff$effort_days)
        rate <- n_events_total / trap_days
        se_rate <- if (n_events_total > 0) sqrt(n_events_total) / trap_days else NA_real_
        p_naive <- 1 - exp(-rate * L_days)
        se_p_naive <- if (!is.na(se_rate))
          exp(-rate * L_days) * L_days * se_rate else NA_real_
        
        species_out[[sp]] <- tibble(
          window_id        = win$window_id,
          window_len       = win$length_d,
          window_start     = win$start_date,
          species          = sp,
          n_events_total   = n_events_total,
          n_occasions_pos  = n_occ_pos,
          n_sites_pos      = n_sites_pos,
          trap_days_window = trap_days,
          n_sites          = nrow(cam_eff),
          # ---- study specific
          trap_array       = trap_array_km,
          latitude         = latitude_deg,
          # ---- occupancy
          psi_hat          = psi_hat,    
          psi_se           = psi_se,
          p_hat            = p_hat,
          p_se             = p_se,
          occ_converged    = conv,
          # ---- TTE
          lambda           = lambda,
          lambda_se        = se_lambda,
          p_tte            = p_tte,
          p_tte_se         = se_p_tte,
          # ---- naive
          rate             = rate,
          rate_se          = se_rate,
          p_naive          = p_naive,
          p_naive_se       = se_p_naive
        )
      } else {
        # species fails thresholds → NULL (handled later)
        species_out[[sp]] <- NULL
      }
    } # species loop
    
    results_list[[i]] <- bind_rows(species_out)
  } # window loop
  
  window_species <- bind_rows(results_list)
  if (nrow(window_species) == 0) {
    window_species <- tibble(window_id = character(),
                             species   = character())
  }
  
  # If nothing passed the thresholds, bail out early
  if (nrow(window_species) == 0)
    return(list(window_species  = tibble(),
                window_camera   = window_camera,
                dropped_species = window_events_sum))
  
  # ── 1. grab the “truth” row (= full-slice window_id == "FULL") ────────────────
  truth_tbl <- window_species %>% 
    filter(window_id == "FULL") %>% 
    select(
      species,
      psi_hat_full   = psi_hat,
      p_hat_full     = p_hat,
      p_tte_full     = p_tte,
      p_naive_full   = p_naive
    )
  
  # ── 2. attach truth + compute offsets, then drop the FULL row ─────────────────
  window_species <- window_species %>% 
    left_join(truth_tbl, by = "species") %>% 
    mutate(
      d_psi     = psi_hat - psi_hat_full,
      d_p       = p_hat   - p_hat_full,
      d_p_tte   = p_tte   - p_tte_full,
      d_p_naive = p_naive - p_naive_full
    ) %>% 
    filter(window_id != "FULL")   %>% # keep only the real windows
    select(-psi_hat_full, -p_hat_full, -p_tte_full, -p_naive_full)
  
  dropped_species <- window_events_sum %>%       # use the aggregated table
    anti_join(window_species, by = c("window_id","species"))
  
  list(
    window_species  = window_species,
    window_camera   = window_camera,
    dropped_species = dropped_species
  )
}



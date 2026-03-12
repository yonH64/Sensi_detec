detection_pdf <- function(archive_path,
                          output_pdf = NULL) {
  ## ── packages ────────────────────────────────────────────────
  suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
    library(gridExtra)
  })
  
  if (is.null(output_pdf)) {
    output_pdf <- file.path(getwd(),
                            paste0(tools::file_path_sans_ext(basename(archive_path)),
                                   "_mammal_detections_over_effort.pdf"))
  }
  
  ## ── load data ───────────────────────────────────────────────
  deploy <- read_csv_src(archive_path, "deployments\\.csv")
  obs    <- read_csv_src(archive_path, "observations\\.csv")
  
  # -------- column detection
  cols <- set_camtrap_cols(obs, deploy)
  list2env(cols, envir = environment())
  
  if (any(is.na(c(sci_col, date_col, depl_id_obs, depl_id_deploy, lat_col, lon_col, start_col, end_col))))
    stop("Missing required columns (species/date/deployment/lat/lon/start/end).")
  
  deploy_ok <- deploy %>% 
    rename(start = !!sym(start_col),
           end   = !!sym(end_col)) %>% 
    mutate(start = ymd_hms(start, tz = "UTC"),
           end   = ymd_hms(end,   tz = "UTC")) %>% 
    filter(end > start,
           end - start <= years(10))               # sanity-check duration
  
  ## ── full month sequence over the dataset span ───────────────
  month_min <- floor_date(min(deploy_ok$start, na.rm = TRUE), "month")
  month_max <- floor_date(max(deploy_ok$end,   na.rm = TRUE), "month")
  all_months <- tibble(month_date = seq.Date(month_min, month_max, by = "month"))
  
  ## ── monthly trap-nights (by year-month, no across-years sum) ─
  effort_month <- deploy_ok %>% 
    mutate(month_seq = map2(start, end,
                            ~ seq.Date(floor_date(.x, "month"),
                                       floor_date(.y, "month"),
                                       by = "month"))) %>% 
    unnest(month_seq) %>% 
    mutate(month_date  = month_seq,
           month_start = month_date,
           month_end   = (month_date + months(1)) - days(1),
           trap_days   = pmax(
             0L,
             as.integer(pmin(as.Date(end),   month_end) -
                          pmax(as.Date(start), month_start)) + 1L)) %>% 
    group_by(month_date) %>% 
    summarise(trap_nights = sum(trap_days), .groups = "drop") %>% 
    right_join(all_months, by = "month_date") %>% 
    mutate(trap_nights = replace_na(trap_nights, 0L))
  
  ## ── tidy observations ───────────────────────────────────────
  clean_empty <- function(df) {
    if (!"is_empty" %in% names(df)) return(rep(FALSE, nrow(df)))
    v <- df$is_empty
    if (is.list(v)) purrr::map_lgl(v, ~ identical(.x, TRUE)) else as.logical(v)
  }
  
  obs_det <- obs %>% 
    mutate(.empty = clean_empty(.)) %>% 
    filter(.empty == FALSE | is.na(.empty)) %>% 
    select(-.empty) %>% 
    filter(!is.na(.data[[sci_col]]),
           !str_starts(.data[[sci_col]],
                       regex("Undefined|Homo\\s+sapiens|Other", ignore_case = TRUE)))
  
  ## ── detections per species × year-month ─────────────────────
  dets_month <- obs_det %>% 
    mutate(month_date = as.Date(lubridate::floor_date(.data[[date_col]], "month"))) %>% 
    count(species = .data[[sci_col]], month_date, name = "n_detections") %>% 
    right_join(all_months, by = "month_date") %>% 
    complete(species, month_date = all_months$month_date, fill = list(n_detections = 0L))
  
  ## ── join effort (keep rate too, in case you need it later) ───
  rate <- dets_month %>% 
    left_join(effort_month, by = "month_date") %>% 
    mutate(month_date = as.Date(month_date),
           det_rate   = if_else(trap_nights > 0, n_detections / trap_nights, NA_real_))
  
  mamm_rate <- rate %>% 
    filter(!str_detect(species, taxa_filter)) %>% 
    # keep species rows that are actually present at least once
    group_by(species) %>% 
    filter(any(n_detections > 0, na.rm = TRUE)) %>% 
    ungroup()
  
  if (nrow(mamm_rate) == 0) {
    warning("No mammal detections after filtering; nothing to plot.")
    return(invisible(NULL))
  }
  
  ## ── build per‑species plots: bars = detections; line = effort ─
  plot_list <- mamm_rate %>% 
    split(.$species) %>% 
    imap(~{
      # scaling so the effort line fits on the bars axis
      ymax  <- max(.x$n_detections, na.rm = TRUE)
      emax  <- max(.x$trap_nights,  na.rm = TRUE)
      scale <- ifelse(emax > 0, ymax / emax, 1)
      
      ggplot(.x, aes(month_date)) +
        geom_col(aes(y = n_detections), width = 25) +
        geom_line(aes(y = trap_nights * scale),
                  linewidth = 0.6) +
        scale_y_continuous(
          name = "Detections",
          sec.axis = sec_axis(~ . / scale, name = "Trap-nights")
        ) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b %Y", expand = c(0,0)) +
        labs(x = NULL, title = .y) +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.minor = element_blank())
    })
  
  ## ── collect rates for overlaying (unchanged, but now year‑month) ─
  this_dataset_rates <- mamm_rate %>% 
    select(species, month = month_date, det_rate) %>% 
    mutate(dataset = tools::file_path_sans_ext(basename(archive_path)))
  
  if (exists("species_rates_bank", envir = .GlobalEnv)) {
    species_rates_bank <- get("species_rates_bank", envir = .GlobalEnv)
  } else {
    species_rates_bank <- list()
  }
  species_rates_bank[[this_dataset_rates$dataset[1]]] <- this_dataset_rates
  assign("species_rates_bank", species_rates_bank, envir = .GlobalEnv)
  
  ## ── export 4 per page ───────────────────────────────────────
  pdf(output_pdf, width = 8.3, height = 11.7)  # A4 portrait
  gridExtra::marrangeGrob(
    grobs = plot_list,
    nrow  = 3,
    ncol  = 1,
    top   = "Monthly detections per mammal species (bars) with effort (line)"
  ) %>% grid::grid.draw()
  dev.off()
  
  message("✅ PDF written to: ", output_pdf)
  invisible(output_pdf)
}



detection_pdf("C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/BFNP_201819")
detection_pdf("C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/BFNP_201920")
detection_pdf("C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/Tim1")
detection_pdf("C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/MICA")
detection_pdf("C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/MICA_N")
detection_pdf("C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/MICA_S")


overlay_plot <- function(species_name) {
  df <- bind_rows(species_rates_bank, .id = "dataset") %>%
    filter(species == species_name) %>%
    mutate(
      moy = factor(month.abb[lubridate::month(month)],
                   levels = month.abb, ordered = TRUE)
    ) %>%
    group_by(dataset, species, moy) %>%
    summarise(det_rate = mean(det_rate, na.rm = TRUE), .groups = "drop")
  
  ggplot(df, aes(x = moy, y = det_rate, colour = dataset, group = dataset)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1) +
    scale_x_discrete(drop = FALSE, limits = month.abb) +
    labs(title = species_name,
         y = "Detections / trap-night",
         x = "Month",
         colour = "Dataset") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# examples
overlay_plot("Capreolus capreolus")   # Roe deer
overlay_plot("Vulpes vulpes")         # Red fox
overlay_plot("Lynx lynx")             # Lynx
overlay_plot("Cervus elaphus")        # Red deer
overlay_plot("Lepus europaeus")       # Hare
overlay_plot("Sus scrofa")            # Wild boar
overlay_plot("Sciurus vulgaris")      # Squirrel
overlay_plot("Meles meles")           # Badger
overlay_plot("Martes martes")         # Stone marten


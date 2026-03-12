archive_path <- "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/BFNP_201819"
archive_path <- "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/BFNP_201920"


species_accum_pdf <- function(archive_path, output_pdf = NULL) {
  
  suppressPackageStartupMessages({
    library(tidyverse); library(lubridate); library(ggrepel)
  })
  
  
  ds <- tools::file_path_sans_ext(basename(archive_path))
  if (!exists("anchors", inherits = TRUE))
    stop("`anchors` data frame not found in the environment.")
  
  # ---- 1) load ------------------------------------------------------------
  deploy <- read_csv_src(archive_path, "deployments\\.csv")
  obs    <- read_csv_src(archive_path, "observations\\.csv")
  
  # -------- column detection
  cols <- set_camtrap_cols(obs, deploy)
  list2env(cols, envir = environment())
  
  if (any(is.na(c(sci_col, date_col, depl_id_obs, depl_id_deploy, lat_col, lon_col, start_col, end_col))))
    stop("Missing required columns (species/date/deployment/lat/lon/start/end).")
  
  # ---- 3) best slice from anchors ----------------------------------------
  best <- anchors %>%
    filter(dataset == ds) %>%
    arrange(desc(coverage_prop), desc(trap_days_12m), longest_gap_d) %>%
    slice_head(n = 1)
  
  if (nrow(best) == 0) stop("No rows for dataset ", ds, " in `anchors`.")
  anchor_start <- as.Date(best$anchor_start)
  anchor_end   <- as.Date(best$anchor_end)
  
  # ---- 4) deployments → monthly effort -----------------------------------
  deploy_m <- deploy %>%
    rename(start = !!sym(start_col), end = !!sym(end_col)) %>%
    mutate(start = ymd_hms(start, tz = "UTC"),
           end   = ymd_hms(end,   tz = "UTC")) %>%
    filter(end > start) %>%
    mutate(month_seq = map2(floor_date(start, "month"),
                            floor_date(end,   "month"),
                            ~ seq.Date(.x, .y, by = "month"))) %>%
    unnest(month_seq) %>%
    transmute(
      month_date  = month_seq,
      month_start = month_seq,
      month_end   = (month_seq + months(1)) - days(1),
      start_date  = as.Date(start),
      end_date    = as.Date(end),
      trap_days   = pmax(0L,
                         as.integer(pmin(end_date,   month_end) -
                                      pmax(start_date, month_start)) + 1L)
    ) %>%
    group_by(month_date) %>%
    summarise(trap_nights = sum(trap_days), .groups = "drop") %>%
    arrange(month_date)
  
  # ---- 5) observations → species list per month --------------------------
  obs_m <- obs %>%
    filter(!is.na(.data[[sci_col]]),
           !str_starts(.data[[sci_col]],
                       regex("Undefined|Homo\\s+sapiens|Other", TRUE))) %>%
    mutate(month_date = floor_date(.data[[date_col]], "month")) %>%
    group_by(month_date) %>%
    summarise(species_list = list(unique(.data[[sci_col]])), .groups = "drop")
  
  # ---- 6) merge & slice to the best 12-month window ----------------------
  month_df <- deploy_m %>%
    full_join(obs_m, by = "month_date") %>%
    arrange(month_date) %>%
    mutate(species_list = dplyr::coalesce(species_list, list(character())))
  
  if (nrow(month_df) < 12) stop("Dataset covers < 12 months.")
  
  if (!exists("month_bank", .GlobalEnv)) month_bank <<- list()
  
  month_bank[[ds]] <<- month_df %>%               # month_df already exists
    mutate(richness = lengths(species_list)) %>%  # species count that month
    select(month_date, richness)                  # keep only what we need
  
  start_m   <- lubridate::floor_date(anchor_start, "month")
  months_sel <- tibble(month_date = seq.Date(start_m, by = "month", length.out = 12))
  
  slice_df <- months_sel %>%
    left_join(month_df, by = "month_date") %>%
    mutate(trap_nights = replace_na(trap_nights, 0L)) %>%
    arrange(month_date)
  
  month_labels <- month.abb[ lubridate::month(slice_df$month_date) ]
  stopifnot(length(month_labels) == 12)
  print(month_labels)
  
  # ---- 7) simulate all 12 start months (overlay curves) ------------------
  # rotate indices 1..12 by start_offset 0..11
  build_curve <- function(start_offset) {
    idx <- ((start_offset + 0:11) %% 12) + 1
    seen <- character(0)
    richness <- map_int(idx, function(i) { seen <<- union(seen, slice_df$species_list[[i]]); length(seen) })
    tibble(
      start_offset = start_offset,
      start_label  = month_labels[((start_offset %% 12) + 1)],
      step         = 1:12,
      richness     = richness
    )
  }
  
  curves <- map_dfr(0:11, build_curve) %>% mutate(is_anchor = start_offset == 0L)
  effort_line <- slice_df %>% transmute(step = 1:12, trap_nights)
  scale_fac <- max(curves$richness, na.rm = TRUE) / max(effort_line$trap_nights, na.rm = TRUE)
  effort_line <- effort_line %>% mutate(trap_scaled = trap_nights * scale_fac)
  
  # ---- 8) plot (single panel, overlays) ----------------------------------
  curves <- curves %>%
    mutate(start_label = factor(start_label, levels = month_labels),
           is_anchor   = start_offset == 0L)
  
  end_pts <- curves %>%
    group_by(start_label) %>%
    summarise(y_end = richness[step == 12][1], .groups = "drop") %>%
    mutate(x_end = 12)
  
  pal <- setNames(scales::hue_pal()(length(month_labels)), month_labels)
  
  y_top <- max(curves$richness, na.rm = TRUE)
  
  p <- ggplot() +
    # all starts (colored)
    geom_line(data = curves,
              aes(step, richness, group = start_label, colour = start_label),
              linewidth = 0.8, alpha = 0.85) +
    # anchor highlighted
    geom_line(data = dplyr::filter(curves, is_anchor),
              aes(step, richness), colour = "black", linewidth = 1.3) +
    # effort (secondary axis; plotted on primary scale after scaling)
    geom_line(data = effort_line,
              aes(step, trap_scaled),
              colour = "firebrick", linetype = "dashed", linewidth = 0.9) +
    # anchor markers + dates
    geom_vline(xintercept = c(1, 12), linetype = "dotted", colour = "grey40") +
    annotate("text", x = 1,  y = y_top,  vjust = -0.6, label = format(anchor_start, "%Y-%m-%d"), size = 3.2) +
    annotate("text", x = 12, y = y_top,  vjust = -0.6, label = format(anchor_end,   "%Y-%m-%d"), size = 3.2) +
    # direct labels at the line ends
    geom_text_repel(data = end_pts,
                    aes(x_end, y_end, label = start_label, colour = start_label),
                    direction = "y", nudge_x = 0.3, hjust = 0,
                    segment.size = 0.2, size = 3, box.padding = 0.2,
                    min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE) +
    scale_colour_manual(values = pal, name = "Start month") +
    scale_x_continuous(breaks = 1:12, labels = month_labels, limits = c(1, 12.8)) +
    scale_y_continuous(
      name = "Cumulative species",
      sec.axis = sec_axis(~ . / scale_fac, name = "Monthly trap-nights")
    ) +
    labs(
      title = paste0("12-month accumulation vs. effort — ", ds),
      subtitle = paste0("Best slice: ", format(anchor_start, "%Y-%m-%d"),
                        " → ", format(anchor_end, "%Y-%m-%d"),
                        " (", month_labels[1], " start highlighted)")
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
  
  
  # ---- 9) save -----------------------------------------------------------
  if (is.null(output_pdf))
    output_pdf <- file.path(getwd(), paste0(ds, "_species_accumulation_best_slice.pdf"))
  
  ggsave(output_pdf, p, width = 8, height = 6)
  message("✅ PDF saved: ", output_pdf)
  
  invisible(list(curves = curves, effort = effort_line, plot = p))
}


species_accum_pdf("C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/BFNP_201819")
species_accum_pdf("C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/BFNP_201920")
species_accum_pdf("C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/MICA")
species_accum_pdf("C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/Tim1")

#####     bar plot of monthly (not cumulative) species richness    #####
all_months <- bind_rows(month_bank, .id = "dataset") %>%   # glue list → tibble
  group_by(dataset, month_date) %>%                        # de-dupe just in case
  summarise(richness = max(richness), .groups = "drop") %>%
  mutate(month_lab = factor(month.abb[month(month_date)],
                            levels = month.abb))           # Jan…Dec order

ggplot(all_months,
       aes(x = month_lab, y = richness, fill = dataset)) +
  geom_col(position = "dodge") +
  labs(x = "Month", y = "Species richness") +
  theme_minimal(base_size = 11)


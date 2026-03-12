library(tidyverse)
library(lubridate)
library(viridis)
library(glue)



plot_offset_heatmaps <- function(window_species_df,
                                 anchor_year = 2021,
                                 rng = .30,
                                 na_col = "grey85") {
  
  library(tidyverse); library(lubridate); library(glue)
  
  yr_start <- as.Date(sprintf("%d-01-01", anchor_year))
  yr_end   <- as.Date(sprintf("%d-12-31", anchor_year))
  
  df <- window_species_df %>% 
    mutate(start_doy  = parse_integer(str_extract(window_id, "(?<=d)\\d{3}")),
           length_d   = parse_integer(str_extract(window_id, "(?<=_L)\\d+")),
           start_date = yr_start + days(start_doy - 1L))
  
  for (sp in unique(df$species)) {
    
    dat_long <- df %>% 
      filter(species == sp) %>% 
      select(dataset, start_date, length_d,
             d_p, d_p_tte, d_p_naive) %>% 
      pivot_longer(d_p:d_p_naive,
                   names_to  = "metric",
                   values_to = "offset") %>% 
      filter(!is.na(offset))
    
    if (nrow(dat_long) == 0) {
      message("No offset values for ", sp, " – skipped.")
      next
    }
    
    facet_labs <- dat_long %>% 
      count(metric, dataset, name = "n_win") %>% 
      summarise(n_sets = n(),
                n_win  = sum(n_win), .by = metric) %>% 
      mutate(label = glue("{metric}  ({n_sets} datasets, {n_win} windows)")) %>% 
      with(setNames(label, metric))
    
    p <- ggplot(dat_long,
                aes(start_date, length_d, fill = offset)) +
      geom_tile() +
      #geom_contour(aes(z = offset), colour = "white",
      #             alpha = .4, binwidth = .02) +
      scale_x_date(breaks  = seq(yr_start, yr_end, by = "1 month"),
                   labels  = scales::label_date("%b"),
                   limits  = c(yr_start, yr_end),
                   expand  = c(0, 0),
                   name    = "Window start (Jan → Dec)") +
      scale_y_continuous(expand = c(0, 0),
                         name   = "Window length (days)") +
      scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                           midpoint = 0, limits = c(-rng, rng),
                           oob = scales::squish,
                           na.value = na_col,
                           name = "Offset\n(window – full)") +
      facet_wrap(~ metric, ncol = 1,
                 labeller = as_labeller(facet_labs)) +
      labs(title = glue("Window vs full-year offset — {sp}")) +
      theme_minimal(base_size = 14) +
      theme(panel.grid  = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p)                        # <-- explicit print!
    
    if (sp != dplyr::last(unique(df$species)))
      readline("Return for next species…")
  }
}

plot_heatmaps <- function(window_species_df,
                                  metric = c("d", "se"),
                                  anchor_year   = 2021,
                                  agg           = c("mean","median"),
                                  trans         = "identity",   # e.g. "sqrt" or "log1p"
                                  na_fill       = "grey10",
                                  na_outline    = "grey10",
                                  drop_fill     = "grey80",
                                  drop_outline  = "grey80") {
  suppressPackageStartupMessages({
    library(tidyverse); library(lubridate); library(glue)
  })
  agg <- match.arg(agg)
  metric <- match.arg(metric) 
  
  cols <- if (metric == "se") c("p_se","p_tte_se","p_naive_se")
  else                 c("d_p","d_p_tte","d_p_naive")
  legend_title <- if (metric == "se") paste(agg, "SE") else paste(agg, "offset")
  plot_title   <- if (metric == "se") "Standard error heatmaps" else "Offset heatmaps"
  caption_txt  <- if (metric == "se") "Dark grey × = dropped by threshold; light grey = no SE"
  else "Dark grey × = dropped by threshold; light grey = no offset"
  
  # ----- anchor year window
  yr_start <- as.Date(sprintf("%d-01-01", anchor_year))
  yr_end   <- as.Date(sprintf("%d-12-31", anchor_year))
 
  # ----- parse window_id -> start_date & length (days)
  base_df <- window_species_df %>%
    mutate(
      start_doy  = readr::parse_integer(stringr::str_extract(window_id, "(?<=d)\\d{3}")),
      length_d   = readr::parse_integer(stringr::str_extract(window_id, "(?<=_L)\\d+")),
      start_date = yr_start + days(start_doy - 1L)
    )
  
  # helper: prepare dropped windows to same grid (species, metric, start_date, length_d)
  prep_dropped <- function(df, metrics_present) {
    if (is.null(df) || nrow(df) == 0) return(tibble())
    out <- df %>%
      mutate(
        start_doy  = readr::parse_integer(stringr::str_extract(window_id, "(?<=d)\\d{3}")),
        length_d   = readr::parse_integer(stringr::str_extract(window_id, "(?<=_L)\\d+")),
        start_date = yr_start + days(start_doy - 1L)
      ) %>%
      select(any_of(c("species","dataset","metric","start_date","length_d"))) %>%
      distinct()
    if (!"metric" %in% names(out)) {
      out <- tidyr::crossing(out, tibble(metric = metrics_present))
    }
    out
  }
  
  # aggregator
  agg_fun <- switch(agg,
                    mean   = function(x) mean(x, na.rm = TRUE),
                    median = function(x) stats::median(x, na.rm = TRUE))
  
  # ----- loop over species
  sp_vec <- unique(base_df$species)
  for (sp in sp_vec) {
    
    dat_long <- base_df %>%
      filter(species == sp) %>%
      select(dataset, species, start_date, length_d, all_of(cols)) %>%
      pivot_longer(all_of(cols), names_to = "metric", values_to = "val")
    
    # one value per cell: aggregate across datasets (default = mean SE)
    grid <- dat_long %>%
      group_by(metric, start_date, length_d) %>%
      summarise(
        n_datasets = sum(!is.na(val)),
        val_agg    = agg_fun(val),   # mean/median across datasets
        .groups    = "drop"
      )
    
    grid_df <- grid  # avoid clash with graphics::grid()
    
    # dashed frontier: max non-wrapping length for each start date
    frontier <- tibble(
      start_date = seq(yr_start, yr_end, by = "1 day"),
      Lmax       = as.integer(yr_end - start_date + 1L)
    )
    
    # flag tiles that *must* wrap the calendar year (top-right region)
    grid_chk <- grid_df %>%
      mutate(Lmax     = as.integer(yr_end - start_date + 1L),
             cal_wrap = length_d > Lmax)
    
    corner_stats <- grid_chk %>%
      summarise(
        n_tiles      = dplyr::n(),
        n_wrap       = sum(cal_wrap, na.rm = TRUE),
        frac_wrap    = mean(cal_wrap, na.rm = TRUE),
        mean_val_all = mean(val_agg, na.rm = TRUE),
        sd_all       = sd(val_agg, na.rm = TRUE),
        mean_val_wrap= mean(val_agg[cal_wrap], na.rm = TRUE),
        sd_wrap      = sd(val_agg[cal_wrap], na.rm = TRUE),
        mean_val_nowrap = mean(val_agg[!cal_wrap], na.rm = TRUE),
        sd_nowrap = sd(val_agg[!cal_wrap], na.rm = TRUE),
        .by = metric
      )
    print(corner_stats)
    
    # identify NA cells (no SE available)
    dat_na <- grid %>% filter(is.na(val_agg)) %>% select(metric, start_date, length_d)
    
    dropped_sp <- prep_dropped(all_dropped_species %>% filter(species == sp),
                               metrics_present = unique(grid$metric))
    if (nrow(dropped_sp)) {
      dat_na <- anti_join(dat_na, dropped_sp %>% select(metric, start_date, length_d),
                          by = c("metric","start_date","length_d"))
    }
    
    # labels for facets
    facet_labs <- dat_long %>% 
      count(metric, dataset, name = "n_win") %>% 
      summarise(n_sets = n(),
                n_win  = sum(n_win), .by = metric) %>% 
      mutate(label = glue("{metric}  ({n_sets} datasets, {n_win} windows)")) %>% 
      with(setNames(label, metric))
    
    if (metric == "se") {
      limits_se <- c(0,0.02)
      } else {
      limits_d <- c(-1.5,1.5)
      }

    # ---- ggplot (handles metric = "se" or "d") ----
    p <- ggplot() +
      # dropped-by-threshold windows (dark grey + small “×”)
      { if (nrow(dropped_sp))
        geom_tile(data = dropped_sp, aes(start_date, length_d),
                  fill = drop_fill, colour = drop_outline, linewidth = 0.25) } +
      { if (nrow(dropped_sp))
        geom_point(data = dropped_sp, aes(start_date, length_d),
                   shape = 4, stroke = 0.35, size = 1.8, colour = drop_outline) } +
      # true NA cells (light grey)
      { if (nrow(dat_na))
        geom_tile(data = dat_na, aes(start_date, length_d),
                  fill = na_fill, colour = na_outline, linewidth = 0.2, alpha = 0.9) } +
      # main tiles
      geom_tile(data = dplyr::filter(grid, !is.na(val_agg)),
                aes(start_date, length_d, fill = val_agg)) +
      geom_line(data = frontier, aes(start_date, Lmax),
                inherit.aes = FALSE, linetype = "dotted",
                linewidth = 0.4, colour = "black") +
      scale_x_date(breaks = seq(yr_start, yr_end, by = "1 month"),
                   labels = scales::label_date("%b"),
                   limits = c(yr_start, yr_end),
                   expand = c(0, 0),
                   name   = "Window start") +
      scale_y_continuous(expand = c(0, 0), name = "Window length (days)") +
      facet_wrap(~ metric, ncol = 1,
                 labeller = if (length(facet_labs)) as_labeller(facet_labs) else "label_value") +
      labs(title   = glue("{plot_title} — {sp}"),
           caption = caption_txt) +
      theme_minimal(base_size = 14) +
      theme(panel.grid  = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right")
    
    # color scale depends on metric choice
    if (metric == "se") {
      p <- p + scale_fill_viridis_c(option = "magma", direction = -1,
                                    trans = trans, limits = limits_se,
                                    oob = scales::squish,
                                    name = legend_title)
    } else {
      p <- p + scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                                    midpoint = 0, limits = limits_d,
                                    oob = scales::squish,
                                    name = legend_title)
    }
    
    print(p)
    if (sp != dplyr::last(sp_vec)) readline("Return for next species…")
  }
}

plot_heatmaps_d_SE <- function(window_species_df,
                          anchor_year   = 2021,
                          agg           = c("mean","median"),
                          trans         = "identity",     # for fill; keep "identity" for signed offsets
                          limits        = c(-1, 1),       # color limits for offsets
                          na_fill       = "grey10",
                          na_outline    = "grey10",
                          drop_fill     = "grey80",
                          drop_outline  = "grey80",
                          alpha_clip_q  = 0.95) {         # clip alpha at this quantile
  suppressPackageStartupMessages({
    library(tidyverse); library(lubridate); library(glue)
  })
  agg <- match.arg(agg)
  
  # ----- anchor year window
  yr_start <- as.Date(sprintf("%d-01-01", anchor_year))
  yr_end   <- as.Date(sprintf("%d-12-31", anchor_year))
  
  # ----- parse window_id -> start_date & length (days)
  base_df <- window_species_df %>%
    mutate(
      start_doy  = readr::parse_integer(stringr::str_extract(window_id, "(?<=d)\\d{3}")),
      length_d   = readr::parse_integer(stringr::str_extract(window_id, "(?<=_L)\\d+")),
      start_date = yr_start + days(start_doy - 1L)
    )
  
  # helper: prepare dropped windows to same grid (species, metric, start_date, length_d)
  prep_dropped <- function(df, metrics_present) {
    if (is.null(df) || nrow(df) == 0) return(tibble())
    out <- df %>%
      mutate(
        start_doy  = readr::parse_integer(stringr::str_extract(window_id, "(?<=d)\\d{3}")),
        length_d   = readr::parse_integer(stringr::str_extract(window_id, "(?<=_L)\\d+")),
        start_date = yr_start + days(start_doy - 1L)
      ) %>%
      select(any_of(c("species","dataset","metric","start_date","length_d"))) %>%
      distinct()
    if (!"metric" %in% names(out)) {
      out <- tidyr::crossing(out, tibble(metric = c("d_p","d_p_tte","d_p_naive")))
    }
    out
  }
  
  # aggregator
  agg_fun <- switch(agg,
                    mean   = function(x) mean(x, na.rm = TRUE),
                    median = function(x) stats::median(x, na.rm = TRUE))
  
  # dropped windows df (optional, if present in env)
  dropped_all <- tryCatch(get("all_dropped_species", inherits = TRUE),
                          error = function(e) NULL)
  
  # ----- loop over species
  sp_vec <- unique(base_df$species)
  for (sp in sp_vec) {
    
    dat_long <- base_df %>%
      filter(species == sp) %>%
      select(dataset, species, start_date, length_d, all_of(c("d_p","d_p_tte","d_p_naive"))) %>%
      pivot_longer(all_of(c("d_p","d_p_tte","d_p_naive")),
                   names_to = "metric", values_to = "val")
    
    # mean offset + dispersion across datasets per cell
    grid <- dat_long %>%
      group_by(metric, start_date, length_d) %>%
      summarise(
        n_datasets = sum(!is.na(val)),
        mean_val   = agg_fun(val),
        sd_val     = ifelse(n_datasets > 1, sd(val, na.rm = TRUE), NA_real_),
        se_val     = ifelse(n_datasets > 1, sd_val / sqrt(n_datasets), NA_real_),
        .groups    = "drop"
      ) %>%
      mutate(disp_val = se_val)  # alpha encodes SE across datasets
    
    # NA cells (no mean estimate)
    dat_na <- grid %>% filter(is.na(mean_val)) %>% select(metric, start_date, length_d)
    
    # dropped windows for this species
    dropped_sp <- if (is.null(dropped_all)) tibble() else
      prep_dropped(dropped_all %>% filter(species == sp),
                   metrics_present = unique(grid$metric))
    
    # avoid overdraw: remove dropped from NA layer
    if (nrow(dropped_sp)) {
      dat_na <- anti_join(dat_na, dropped_sp %>% select(metric, start_date, length_d),
                          by = c("metric","start_date","length_d"))
    }
    
    # facet labels
    facet_labs <- grid %>%
      summarise(n_cells = sum(!is.na(mean_val)), .by = metric) %>%
      mutate(label = glue("{metric}  ({n_cells} cells)")) %>%
      { setNames(.$label, .$metric) }
    
    # alpha clipping to reduce dominance of extreme SEs
    alpha_clip <- stats::quantile(grid$disp_val, probs = alpha_clip_q, na.rm = TRUE)
    
    # ---- ggplot (offsets with variance in alpha) ----
    p <- ggplot() +
      # dropped-by-threshold windows (dark grey + small “×”)
      { if (nrow(dropped_sp))
        geom_tile(data = dropped_sp, aes(start_date, length_d),
                  fill = drop_fill, colour = drop_outline, linewidth = 0.25) } +
      { if (nrow(dropped_sp))
        geom_point(data = dropped_sp, aes(start_date, length_d),
                   shape = 4, stroke = 0.35, size = 1.8, colour = drop_outline) } +
      # true NA cells
      { if (nrow(dat_na))
        geom_tile(data = dat_na, aes(start_date, length_d),
                  fill = na_fill, colour = na_outline, linewidth = 0.2, alpha = 0.9) } +
      # mean offset (fill) + SE across datasets (alpha)
      geom_tile(data = dplyr::filter(grid, !is.na(mean_val)),
                aes(start_date, length_d, fill = mean_val, alpha = disp_val)) +
      scale_x_date(breaks = seq(yr_start, yr_end, by = "1 month"),
                   labels = scales::label_date("%b"),
                   limits = c(yr_start, yr_end),
                   expand = c(0, 0),
                   name   = "Window start (Jan → Dec)") +
      scale_y_continuous(expand = c(0, 0), name = "Window length (days)") +
      scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                           midpoint = 0, limits = limits, trans = trans,
                           oob = scales::squish,
                           name = paste(agg, "offset")) +
      scale_alpha_continuous(range = c(0.35, 1),
                             limits = c(0, alpha_clip),
                             oob = scales::squish,
                             name = "SE across datasets",
                             guide = guide_legend(order = 2)) +
      facet_wrap(~ metric, ncol = 1,
                 labeller = if (length(facet_labs)) as_labeller(facet_labs) else "label_value") +
      labs(title   = glue("Offset heatmaps — {sp}"),
           caption = "Color = mean(offset) • Alpha = SE across datasets (clipped at 95th pct).  Dark grey × = dropped; light grey = no data") +
      theme_minimal(base_size = 14) +
      theme(panel.grid  = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right")
    
    print(p)
    if (sp != dplyr::last(sp_vec)) readline("Return for next species…")
  }
}


plot_offset_heatmaps(all_window_species)

## run
plot_heatmaps(all_window_species, metric = "d", agg = "mean")
plot_heatmaps(all_window_species, metric = "se", agg = "mean")

plot_heatmaps_d_SE(all_window_species,
                   agg = "mean")

plot_effort_vs_est <- function(window_species_df,
                               point_alpha = .6,
                               min_trapdays = 1) {
  
  suppressPackageStartupMessages(library(tidyverse))
  
  # ── keep only windows that have an estimate and >0 effort ──────────────
  long_df <- window_species_df %>% 
    filter(trap_days_window >= min_trapdays) %>% 
    select(dataset, species, trap_days_window,
           p_hat, p_tte, p_naive) %>% 
    pivot_longer(p_hat:p_naive,
                 names_to  = "metric",
                 values_to = "prob") %>% 
    filter(!is.na(prob))
  
  # ── scatter with log-x & one facet per metric ───────────────────────────
  ggplot(long_df,
         aes(trap_days_window, prob,
             colour = metric)) +
    geom_point(alpha = point_alpha, size = 2, show.legend = FALSE) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                se = FALSE, colour = "black", linewidth = .7) +
    scale_x_log10(labels = scales::comma_format(accuracy = 1),
                  name   = "Trap-days in window (log scale)") +
    scale_y_continuous(limits = c(0, 1),
                       name   = "Estimated detection probability") +
    scale_colour_brewer(palette = "Set1") +
    facet_wrap(~ metric, ncol = 1,
               labeller = as_labeller(c(
                 p_hat   = "Occupancy p̂",
                 p_tte   = "TTE-derived p",
                 p_naive = "Naïve rate p"
               ))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
}
# call it
plot_effort_vs_est(all_window_species)

#### ── Histogram of detection metrics ─────────────────────────────────── ####
library(tidyverse)

# choose the variables ────
vars <- c(psi_hat = "ψ  (occupancy)",
          p_hat   = "p  (occupancy det.)",
          p_tte   = "P  (time-to-event)",
          p_naive = "P  (naïve rate)")

# gather, drop missing, restrict to [0,1] ────
hist_df <- all_window_species %>%
  select(any_of(names(vars))) %>%          # keep only the 4 columns
  pivot_longer(cols   = everything(),
               names_to  = "metric",
               values_to = "value") %>%
  mutate(metric = recode(metric, !!!vars)) %>%   # nicer facet labels
  filter(!is.na(value) & value >= 0 & value <= 1)

# plot ────────────────────────────────────────────────
ggplot(hist_df, aes(value)) +
  geom_histogram(colour = "white", fill = "#3182bd", bins = 40) +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") +
  scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1.11, .2)) +
  labs(x = "Estimated probability", y = "Count (windows × species)",
       title = "Distribution of detection-probability estimates") +
  theme_minimal(base_size = 13)



#### ── Top-N windows” calendar strip ──────────────────────────────────── ####
anchor_year <- 2021
year0 <- ymd(sprintf("%d-01-01", anchor_year), tz = "UTC")

strip_df <- all_window_species %>%                     # ← your table
  mutate(
    # decode window_id once ─────────────────────────────
    start_doy = parse_integer(str_extract(window_id, "(?<=d)\\d{3}")),
    length_d  = parse_integer(str_extract(window_id, "(?<=_L)\\d+")),
    window_start = as.Date(year0 + days(start_doy - 1L)),
    window_len   = length_d,
    month        = floor_date(window_start, "month")
  ) %>% 
  group_by(month) %>%                                 # rank within month
  mutate(rank = rank(p_hat, ties.method = "first")) %>% 
  slice_min(rank, n = 3, with_ties = TRUE) %>%        # keep Top-3
  ungroup() %>% 
  mutate(day0 = window_start,
         day1 = window_start + days(window_len - 1L))

## ── split windows that wrap 31 Dec → 1 Jan ─────────────────────────
strip_df_fixed <- strip_df %>% 
  mutate(end_doy = start_doy + length_d - 1L,
         wraps   = end_doy > 365) %>% 
  tidyr::uncount(if_else(wraps, 2L, 1L), .remove = FALSE) %>% 
  group_by(across(everything())) %>% 
  mutate(
    part      = row_number(),
    day0_fix  = if_else(part == 1 & wraps, start_doy, 1L),
    day1_fix  = if_else(part == 1 & wraps, 365L,
                        if_else(part == 2 & wraps, end_doy - 365L, end_doy)),
    date0     = as.Date(sprintf("%d-%03d", anchor_year, day0_fix), "%Y-%j"),
    date1     = as.Date(sprintf("%d-%03d", anchor_year, day1_fix), "%Y-%j")
  ) %>% 
  ungroup()

## ── plot ───────────────────────────────────────────────────────────
ggplot(strip_df_fixed) +
  geom_segment(aes(x = date0, xend = date1,
                   y = species, yend = species, colour = factor(rank)),
               size = 6, lineend = "butt") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b",
               expand = c(0, 0), name = "Calendar month") +
  scale_colour_viridis_d(option = "C", end = .7, direction = -1,
                         name = "Top-3 rank") +
  labs(y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))



#### ──  Effort-vs-Accuracy curve ──────────────────────────────────────── ####
curve_df <- all_window_species %>% 
  mutate(delta_p = abs(p_hat - max(p_hat, na.rm=TRUE))) %>%
  mutate(trap_days = length_d * n_sites)   # rough effort proxy

ggplot(curve_df,
       aes(x = trap_days, y = delta_p,
           colour = method)) +             # if you stored method
  geom_point(alpha = .4, size = 2) +
  geom_smooth(se = FALSE) +
  scale_x_log10(labels = scales::comma) +
  labs(x = "Trap-days (log-scale)",
       y = "Δ p versus best (%)",
       colour = "Metric",
       title = "Effort vs detection-accuracy") +
  facet_wrap(~ species, scales="free_x") +
  theme_minimal(base_size=12)


#### ── Dropped-count bar ──────────────────────────────────────────────── ####

dropped_df <- all_dropped_species %>% count(species)

ggplot(dropped_df,
       aes(reorder(species, n), n)) +
  geom_col(fill = "firebrick") +
  coord_flip() +
  labs(x = NULL, y = "# windows dropped",
       title = "Threshold failures by species") +
  theme_minimal(base_size = 12)


#### ── Heat-map of dropped species windows ────────────────────────────── ####

plot_dropped_heatmaps <- function(drop_df, anchor_year = 2021) {
  
  year0 <- lubridate::ymd(sprintf("%d-01-01", anchor_year))
  
  drop_df %>% 
    mutate(start_doy = readr::parse_integer(str_extract(window_id,"(?<=d)\\d+")),
           len_d     = readr::parse_integer(str_extract(window_id,"(?<=_L)\\d+")),
           start_dt  = as.Date(year0 + days(start_doy - 1))) %>% 
    count(species, start_dt, len_d, name = "n_drop") %>% 
    ggplot(aes(start_dt, len_d, fill = n_drop)) +
    geom_tile() +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_viridis_c(option = "B", name = "# datasets\n(dropped)") +
    facet_wrap(~ species) +
    labs(x = "Window start", y = "Window length (days)",
         title = "Where each species fails the thresholds") +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_dropped_heatmaps(all_dropped_species)

#### ── Effort vs estimate scatter plot ────────────────────────────────── ####
plot_effort_vs_est <- function(window_species_df,
                               point_alpha = .6,
                               min_trapdays = 1) {
  
  suppressPackageStartupMessages(library(tidyverse))
  
  # ── keep only windows that have an estimate and >0 effort ──────────────
  long_df <- window_species_df %>% 
    filter(trap_days_window >= min_trapdays) %>% 
    select(dataset, species, trap_days_window,
           p_hat, p_tte, p_naive) %>% 
    pivot_longer(p_hat:p_naive,
                 names_to  = "metric",
                 values_to = "prob") %>% 
    filter(!is.na(prob))
  
  # ── scatter with log-x & one facet per metric ───────────────────────────
  ggplot(long_df,
         aes(trap_days_window, prob,
             colour = metric)) +
    geom_point(alpha = point_alpha, size = 2, show.legend = FALSE) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                se = FALSE, colour = "black", linewidth = .7) +
    scale_x_log10(labels = scales::comma_format(accuracy = 1),
                  name   = "Trap-days in window (log scale)") +
    scale_y_continuous(limits = c(0, 1),
                       name   = "Estimated detection probability") +
    scale_colour_brewer(palette = "Set1") +
    facet_wrap(~ metric, ncol = 1,
               labeller = as_labeller(c(
                 p_hat   = "Occupancy p̂",
                 p_tte   = "TTE-derived p",
                 p_naive = "Naïve rate p"
               ))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
}

# call it
plot_effort_vs_est(all_window_species)



library(tidyverse)   # dplyr, tidyr, ggplot2
library(GGally)      # nice scatter-plot matrices (ggpairs)

#### ── Correlations between detection metrics ─────────────────────────── ####
###  1. pick the four columns & drop rows with any NA ────────────
met_df <- all_window_species %>% 
  select(psi_hat, p_hat, p_tte, p_naive) %>% 
  drop_na()

##  2. overall correlation matrix ────────────
cor_mat <- cor(met_df, use = "pairwise.complete.obs")
cor_mat
#>            psi_hat     p_hat     p_tte   p_naive
#> psi_hat  1.0000000 0.83 ...  ...

##  3. quick visual check ────────────────────
corrplot::corrplot(cor_mat, method = "circle")   # if you prefer base graphics

##  4. correlations per species (optional) ────────────────
cor_by_species <- all_window_species %>% 
  group_by(species) %>% 
  summarise(corr = list(cor(select(cur_data(),
                                   psi_hat, p_hat, p_tte, p_naive),
                            use = "pairwise.complete.obs")),
            .groups = "drop")

# show the matrix for one species
cor_by_species$corr[[1]]


#### ── Correlations between psi and effort / design covariates ────────── ####
library(tidyverse)

vars <- c("n_events_total", "n_occasions_pos", "n_sites_pos",
          "trap_days_window", "trap_array", "window_len")

all_window_species %>% 
  select(psi_hat, all_of(vars)) %>% 
  ggpairs(lower = list(continuous = wrap("smooth", method = "gam", se = FALSE)),
          progress = FALSE)

all_window_species %>% 
  ggplot()

archive_path <- "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/BFNP_201819"
archive_path <- "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/BFNP_201920"
archive_path <- "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/Tim1"
archive_path <- "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/MICA"



summarise_species_events <- function(
    archive_path,
    independence_mins = 30
) {
  suppressPackageStartupMessages({
    library(tidyverse); library(lubridate); library(purrr)
  })
    
  ## ── load deployments & observations ─────────────────────────
  deploy <- read_csv_src(archive_path, "deployments\\.csv")
  obs    <- read_csv_src(archive_path, "observations\\.csv")
  
  # # -------- column detection
  cols <- set_camtrap_cols(obs, deploy)
  list2env(cols, envir = environment())
  
  if (any(is.na(c(sci_col, date_col, depl_id_obs, depl_id_deploy, lat_col, lon_col, start_col, end_col))))
    stop("Missing required columns (species/date/deployment/lat/lon/start/end).")
  
  ## ── mammals only & collapse to independent events ───────────
  obs_det <- obs %>% 
    filter(!is.na(.data[[sci_col]]),
           !str_detect(.data[[sci_col]], taxa_filter)) %>% 
    mutate(timestamp = ymd_hms(.data[[date_col]], tz = "UTC")) %>% 
    arrange(.data[[depl_id_obs]], .data[[sci_col]], timestamp) %>% 
    group_by(.data[[depl_id_obs]], .data[[sci_col]]) %>% 
    mutate(gap = timestamp - lag(timestamp),
           new_event = is.na(gap) | gap > minutes(independence_mins)) %>% 
    filter(new_event) %>% 
    ungroup()
  
  ## ── rarity metrics ──────────────────────────────────────────
  total_cameras <- n_distinct(deploy[[depl_id_deploy]])

  names(deploy)                       # see actual column names
  depl_id_obs                             # what the helper guessed from obs
  setdiff(names(deploy), names(obs))  # reveals fields that exist only in deploy
  
  
  summary_tbl <- obs_det %>% 
    group_by(species = .data[[sci_col]]) %>% 
    summarise(
      n_events   = n(),
      n_cameras  = n_distinct(.data[[depl_id_obs]]),
      prop_cameras = n_cameras / total_cameras,
      .groups = "drop"
    ) %>% 
    arrange(desc(n_events)) %>% 
    mutate(rank = row_number(),
           dataset = tools::file_path_sans_ext(basename(archive_path)))
  }


## 1. List your dataset paths
dataset_paths <- c(
  "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/BFNP_201819",
  "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/BFNP_201920",
  "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/Tim1",
  "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/MICA")
## 2. Build a single master table
species_rarity <- map_df(dataset_paths, summarise_species_events)

## Rank–abundance curve for one dataset -------------------------
library(ggplot2)
library(ggrepel)

plot_rank_abundance <- function(df, dataset_id) {
  dat <- df %>% dplyr::filter(dataset == dataset_id)
  
  ggplot(dat, aes(rank, n_events)) +
    geom_line(colour = "steelblue", linewidth = 1) +
    geom_point(size = 2, colour = "steelblue") +
    scale_x_continuous(
      breaks = seq_len(max(dat$rank)),
      minor_breaks = NULL
    ) +
    scale_y_log10(
      breaks  = 10^(0:3),
      labels  = scales::comma
    ) +
    geom_hline(
      yintercept = 20,
      linetype   = "dashed",
      colour     = "grey60"
    ) +
    ## label **every** point  -------------------------------
  geom_text_repel(
    aes(label = species),
    size              = 3,
    min.segment.length = 0,
    box.padding        = 0.3,
    max.overlaps       = Inf          # don’t drop any labels
  ) +
    labs(
      title = paste(dataset_id, "– Rank–abundance"),
      x     = "Species rank (common → rare)",
      y     = "Independent events (log10)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank()
    )
}


plot_rank_abundance(species_rarity, "BFNP_201819")
plot_rank_abundance(species_rarity, "BFNP_201920")
plot_rank_abundance(species_rarity, "Tim1")
plot_rank_abundance(species_rarity, "MICA")

## Heat-map across datasets ------------------------------------
plot_heatmap_prevalence <- function(df) {
  # global order by total events
  species_order <- df %>% 
    group_by(species) %>% 
    summarise(total_events = sum(n_events), .groups = "drop") %>% 
    arrange(total_events) %>% 
    pull(species)
  
  df %>% 
    mutate(species = factor(species, levels = species_order)) %>% 
    ggplot(aes(dataset, species, fill = prop_cameras)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Prop. cameras", na.value = "grey90") +
    labs(x = "Dataset", y = NULL,
         title = "Camera prevalence heat-map per species") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_heatmap_prevalence(species_rarity)


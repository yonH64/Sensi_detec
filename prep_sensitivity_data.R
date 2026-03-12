# prep_sensitivity_data.R
# Prepares the sensitivity surface analysis data from the window metric outputs.
# Input:  all_window_species.rds, all_window_richness.rds (from Full1.R)
#         species_meta.csv, dataset_metadata.csv
# Output: sensitivity_species_data.rds, sensitivity_richness_data.rds

library(tidyverse)

# ── Load source data ──────────────────────────────────────────────────────────
all_window_species  <- readRDS("all_window_species.rds")
all_window_richness <- readRDS("all_window_richness.rds")
species_meta  <- read.csv("species_meta.csv")
dataset_meta  <- read.csv("dataset_metadata.csv")

# BE-Leuven is missing from dataset_metadata.csv (extracted manually)
be_leuven_env <- tibble(
  bio4_temp_seasonality  = 559.2,
  bio12_annual_precip    = 774,
  bio15_precip_seasonality = 11.2,
  ndvi_amplitude         = 0.449
)

# ── Helper: join environmental covariates with slice-name matching ────────────
join_env_covariates <- function(df, meta, be_fix) {
  env_cols <- c("bio4_temp_seasonality", "bio12_annual_precip",
                "bio15_precip_seasonality", "ndvi_amplitude")
  # centroid_lat follows the same join logic as the env covariates
  join_cols <- c(env_cols, "centroid_lat")
  
  df <- df |>
    mutate(dataset_base = str_remove(dataset, "_slice\\d+$"))
  
  # First-pass join on exact dataset name
  df <- df |>
    left_join(
      meta |> select(dataset_name, all_of(join_cols)) |>
        rename(dataset = dataset_name),
      by = "dataset"
    )
  
  # Second-pass: fill missing from base-name match (sliced datasets)
  lookup <- meta |>
    select(dataset_name, all_of(join_cols)) |>
    rename(dataset_base = dataset_name)
  
  for (col in join_cols) {
    missing_idx <- which(is.na(df[[col]]))
    if (length(missing_idx) > 0) {
      fill_vals <- lookup[[col]][match(df$dataset_base[missing_idx], lookup$dataset_base)]
      df[[col]][missing_idx] <- fill_vals
    }
  }
  
  # Third-pass: BE-Leuven manual fix (env covariates only)
  for (col in env_cols) {
    be_rows <- grepl("^BE-Leuven", df$dataset) & is.na(df[[col]])
    if (any(be_rows)) df[[col]][be_rows] <- be_fix[[col]]
  }
  
  df
}

# ── Species-level data ────────────────────────────────────────────────────────
sens_species <- all_window_species |>
  # Use sliding windows only (exclude named protocol windows)
  filter(!window_id %in% c("SNAP_EU_CORE", "SNAP_EU_BUFFER", "EOW_EARLY", "EOW_LATE")) |>
  mutate(
    # Window design parameters
    day_start  = as.integer(str_extract(window_id, "(?<=d)\\d+")),
    day_center = (day_start + window_len / 2) %% 365,
    
    # Absolute deviations
    abs_d_lambda       = abs(d_lambda),
    abs_d_rate         = abs(d_rate),
    abs_d_matched_rate = abs(d_matched_rate),
    
    # Effort covariates (log-transformed)
    l_trapdays = log1p(trap_days_window),
    l_nsites   = log1p(n_sites),
    
    # Grouping
    ds_sp = interaction(dataset, species, drop = TRUE)
  ) |>
  # Join species traits
  left_join(
    species_meta |> select(species, guild_major, guild_minor_habitat, guild_minor_diet),
    by = "species"
  ) |>
  # Join environmental covariates
  join_env_covariates(dataset_meta, be_leuven_env) |>
  mutate(
    # Standardise continuous covariates
    s_bio4       = as.numeric(scale(bio4_temp_seasonality)),
    s_latitude   = as.numeric(scale(latitude)),
    s_trap_array = as.numeric(scale(log1p(trap_array))),
    
    # Factor coding
    guild_major         = factor(guild_major),
    guild_minor_habitat = factor(guild_minor_habitat),
    guild_minor_diet    = factor(guild_minor_diet),
    species_f  = factor(species),
    dataset_f  = factor(dataset),
    ds_sp_f    = factor(ds_sp)
  )

stopifnot(sum(is.na(sens_species$bio4_temp_seasonality)) == 0)
stopifnot(sum(is.na(sens_species$guild_major)) == 0)

# ── Richness-level data ──────────────────────────────────────────────────────
sens_richness <- all_window_richness |>
  filter(!window_id %in% c("SNAP_EU_CORE", "SNAP_EU_BUFFER", "EOW_EARLY", "EOW_LATE")) |>
  mutate(
    day_start  = as.integer(str_extract(window_id, "(?<=d)\\d+")),
    day_center = (day_start + window_len / 2) %% 365,
    
    abs_d_sr       = abs(d_sr),
    abs_d_sr_raref = abs(d_sr_raref),
    
    l_trapdays = log1p(trap_days_window),
    l_nsites   = log1p(n_sites),
    
    dataset_f = factor(dataset)
  ) |>
  join_env_covariates(dataset_meta, be_leuven_env) |>
  mutate(
    s_bio4     = as.numeric(scale(bio4_temp_seasonality)),
    s_latitude = as.numeric(scale(centroid_lat))
  )

stopifnot(sum(is.na(sens_richness$bio4_temp_seasonality)) == 0)

# ── Save ─────────────────────────────────────────────────────────────────────
saveRDS(sens_species,  "sensitivity_species_data.rds")
saveRDS(sens_richness, "sensitivity_richness_data.rds")

cat("Species-level:", nrow(sens_species), "rows,",
    n_distinct(sens_species$species), "species,",
    n_distinct(sens_species$dataset), "datasets\n")
cat("Richness-level:", nrow(sens_richness), "rows,",
    n_distinct(sens_richness$dataset), "datasets\n")

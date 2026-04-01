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
dataset_overview_metrics  <- read.csv("dataset_metadata.csv")

# ── Helper: join centroid_lat with slice-name matching ────────────────────────
join_centroid_lat <- function(df, meta) {
  df <- df |>
    mutate(dataset_base = str_remove(dataset, "_slice\\d+$"))

  # First-pass join on exact dataset name
  df <- df |>
    left_join(
      meta |> select(dataset_name, centroid_lat) |>
        rename(dataset = dataset_name),
      by = "dataset"
    )

  # Second-pass: fill missing from base-name match (sliced datasets)
  lookup <- meta |>
    select(dataset_name, centroid_lat) |>
    rename(dataset_base = dataset_name)

  missing_idx <- which(is.na(df$centroid_lat))
  if (length(missing_idx) > 0) {
    fill_vals <- lookup$centroid_lat[match(df$dataset_base[missing_idx], lookup$dataset_base)]
    df$centroid_lat[missing_idx] <- fill_vals
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
    # 1 obs out of ~222k has exact zero; offset for Gamma regression
    abs_d_rate         = pmax(abs(d_rate), 1e-10),
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
  mutate(
    # Standardise continuous covariates
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

stopifnot(sum(is.na(sens_species$guild_major)) == 0)

# ── Shared scaling constants (species-level = canonical reference) ─────────
# Both species and richness data use the same centering/scaling for latitude
# so that s_latitude = 0 maps to the same physical latitude in both models.
lat_center <- mean(sens_species$latitude)
lat_sd     <- sd(sens_species$latitude)
trap_array_center <- mean(log1p(sens_species$trap_array))
trap_array_sd     <- sd(log1p(sens_species$trap_array))

# Re-apply species scaling with explicit constants (for reproducibility)
sens_species <- sens_species |>
  mutate(
    s_latitude   = (latitude - lat_center) / lat_sd,
    s_trap_array = (log1p(trap_array) - trap_array_center) / trap_array_sd
  )

# Save canonical constants
write.csv(
  data.frame(
    covariate = c("latitude", "log1p(trap_array)",
                   "log1p(trap_days_window)", "log1p(n_sites)"),
    transformation = c("(x - mean) / sd", "(log1p(x) - mean) / sd",
                        "log1p(x)", "log1p(x)"),
    center = c(round(lat_center, 4), round(trap_array_center, 4), NA, NA),
    scale  = c(round(lat_sd, 4), round(trap_array_sd, 4), NA, NA),
    model_name = c("s_latitude", "s_trap_array", "l_trapdays", "l_nsites")
  ),
  "covariate_scaling_constants.csv", row.names = FALSE
)

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
  join_centroid_lat(dataset_overview_metrics) |>
  mutate(
    # Use SAME scaling constants as species data
    s_latitude = (centroid_lat - lat_center) / lat_sd
  )

# ── Save ─────────────────────────────────────────────────────────────────────
saveRDS(sens_species,  "sensitivity_species_data.rds")
saveRDS(sens_richness, "sensitivity_richness_data.rds")

cat("Species-level:", nrow(sens_species), "rows,",
    n_distinct(sens_species$species), "species,",
    n_distinct(sens_species$dataset), "datasets\n")
cat("Richness-level:", nrow(sens_richness), "rows,",
    n_distinct(sens_richness$dataset), "datasets\n")

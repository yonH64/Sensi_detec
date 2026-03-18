# threshold_sensitivity_sites.R
# ============================================================================
# Sensitivity analysis: does the min_sites_pos species-inclusion threshold
# (5 vs 10 vs 15) affect the main results?
#
# Tests three scenarios:
#   Baseline (min_sites_pos = 5):  existing data
#   Strict   (min_sites_pos = 10): subset of existing data
#   V.Strict (min_sites_pos = 15): subset of existing data
#
# For each scenario, fits the M6 species-specific GAM surface for |d_lambda|
# and compares: deviance explained, parametric coefficients, duration curves,
# seasonal profiles, and full 2D surface correlations.
#
# A lenient threshold (min_sites_pos = 3) was evaluated diagnostically but
# skipped: it would add only ~3,100 rows (2.5%) and zero new species after
# the 50-row minimum filter, so a pipeline re-run is not justified.
#
# Prerequisites:
#   - all_window_species.rds, all_dropped_species.rds on disk
#   - sensitivity_species_data.rds on disk (baseline)
#   - sensitivity_gam_models.rds on disk (baseline M6)
#   - species_meta.csv, dataset_metadata.csv on disk
#
# Output:
#   threshold_sensitivity_sites_summary.csv
#   figures/threshold_sensitivity_sites_comparison.pdf
# ============================================================================

library(tidyverse)
library(mgcv)
library(patchwork)

N_THREADS <- 4
cc_knots  <- list(day_start = c(0, 365))

cat("=== min_sites_pos Threshold Sensitivity Analysis ===\n\n")

# ── Shared data prep utilities ───────────────────────────────────────────────

species_meta  <- read.csv("species_meta.csv")
dataset_meta  <- read.csv("dataset_metadata.csv")

be_leuven_env <- tibble(
  bio4_temp_seasonality   = 559.2,
  bio12_annual_precip     = 774,
  bio15_precip_seasonality = 11.2,
  ndvi_amplitude          = 0.449
)

join_env_covariates <- function(df, meta, be_fix) {
  env_cols  <- c("bio4_temp_seasonality", "bio12_annual_precip",
                 "bio15_precip_seasonality", "ndvi_amplitude")
  join_cols <- c(env_cols, "centroid_lat")
  
  df <- df |> mutate(dataset_base = str_remove(dataset, "_slice\\d+$"))
  
  df <- df |>
    left_join(
      meta |> select(dataset_name, all_of(join_cols)) |>
        rename(dataset = dataset_name),
      by = "dataset"
    )
  
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
  
  for (col in env_cols) {
    be_rows <- grepl("^BE-Leuven", df$dataset) & is.na(df[[col]])
    if (any(be_rows)) df[[col]][be_rows] <- be_fix[[col]]
  }
  
  df
}

prep_species_data <- function(ws_data, min_rows_per_species = 50) {
  out <- ws_data |>
    filter(!window_id %in% c("SNAP_EU_CORE", "SNAP_EU_BUFFER",
                              "EOW_EARLY", "EOW_LATE")) |>
    mutate(
      day_start          = as.integer(str_extract(window_id, "(?<=d)\\d+")),
      day_center         = (day_start + window_len / 2) %% 365,
      abs_d_lambda       = abs(d_lambda),
      abs_d_rate         = abs(d_rate),
      abs_d_matched_rate = abs(d_matched_rate),
      l_trapdays         = log1p(trap_days_window),
      l_nsites           = log1p(n_sites),
      ds_sp              = interaction(dataset, species, drop = TRUE)
    ) |>
    left_join(
      species_meta |> select(species, guild_major, guild_minor_habitat, guild_minor_diet),
      by = "species"
    ) |>
    join_env_covariates(dataset_meta, be_leuven_env)
  
  # Drop species with missing guild or too few rows
  spp_counts <- out |> count(species)
  keep_spp   <- spp_counts |> filter(n >= min_rows_per_species) |> pull(species)
  
  out <- out |>
    filter(!is.na(guild_major), species %in% keep_spp) |>
    mutate(
      s_bio4       = as.numeric(scale(bio4_temp_seasonality)),
      s_latitude   = as.numeric(scale(latitude)),
      s_trap_array = as.numeric(scale(log1p(trap_array))),
      guild_major         = factor(guild_major),
      guild_minor_habitat = factor(guild_minor_habitat),
      guild_minor_diet    = factor(guild_minor_diet),
      species_f  = factor(species),
      dataset_f  = factor(dataset),
      ds_sp_f    = factor(ds_sp)
    )
  
  out
}

fit_m6_lambda <- function(data, label = "") {
  cat(sprintf("  Fitting M6 abs_d_lambda %-30s ... ", label))
  t0 <- proc.time()
  fit <- bam(
    abs_d_lambda ~
      te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6),
         by = species_f) +
      species_f +
      s_bio4 + l_trapdays + l_nsites + s_latitude + s_trap_array +
      s(ds_sp_f, bs = "re"),
    data     = data,
    family   = Gamma(link = "log"),
    method   = "fREML",
    discrete = TRUE,
    nthreads = N_THREADS,
    knots    = cc_knots
  )
  elapsed <- round((proc.time() - t0)[3], 1)
  cat(sprintf("done (%ss, dev.expl = %.1f%%)\n", elapsed,
              100 * summary(fit)$dev.expl))
  fit
}


# ── 1. Diagnostic: Characterise exclusion margins ────────────────────────────

cat("── Step 1: Diagnostic ─────────────────────────────────────────\n\n")

all_window_species  <- readRDS("all_window_species.rds")
all_dropped_species <- readRDS("all_dropped_species.rds")

cat(sprintf("Retained: %d rows, %d species\n",
            nrow(all_window_species), n_distinct(all_window_species$species)))
cat(sprintf("Dropped:  %d rows, %d species\n\n",
            nrow(all_dropped_species), n_distinct(all_dropped_species$species)))

# n_sites_pos distribution in retained data
cat("Retained data — n_sites_pos breakdown:\n")
all_window_species |>
  mutate(bin = cut(n_sites_pos,
                   breaks = c(0, 3, 5, 7, 10, 15, 20, Inf),
                   labels = c("1-3", "4-5", "6-7", "8-10", "11-15", "16-20", "20+"),
                   right = TRUE)) |>
  count(bin) |> mutate(pct = round(100 * n / sum(n), 1)) |> print()

# n_sites_pos distribution in dropped data
cat("\nDropped data — n_sites_pos breakdown:\n")
all_dropped_species |>
  mutate(bin = cut(n_sites_pos,
                   breaks = c(0, 2, 3, 5, 10, Inf),
                   labels = c("1-2", "3", "4-5", "6-10", "10+"),
                   right = TRUE)) |>
  count(bin) |> mutate(pct = round(100 * n / sum(n), 1)) |> print()

# Potential gains at min_sites_pos = 3
# (dropped only because of min_sites_pos, with n_sites_pos >= 3)
potential_at_3 <- all_dropped_species |>
  filter(n_sites_pos >= 3, n_sites_pos < 5,
         !failed_min_events, !failed_min_occasions_pos)
potential_spp_3 <- unique(potential_at_3$species)
spp_at_5 <- unique(all_window_species$species)
spp_only_new_3 <- setdiff(potential_spp_3, spp_at_5)

cat(sprintf("\nLenient (min_sites_pos=3): %d additional species-window combos, %d species already in data, %d new species: %s\n",
            nrow(potential_at_3),
            length(intersect(potential_spp_3, spp_at_5)),
            length(spp_only_new_3),
            if (length(spp_only_new_3) > 0) paste(spp_only_new_3, collapse = ", ") else "(none)"))

# Losses at stricter thresholds
cat(sprintf("Strict (min_sites_pos=10): lose %d rows (%.1f%%), %d species retained\n",
            sum(all_window_species$n_sites_pos < 10),
            100 * mean(all_window_species$n_sites_pos < 10),
            n_distinct(filter(all_window_species, n_sites_pos >= 10)$species)))

cat(sprintf("V.Strict (min_sites_pos=15): lose %d rows (%.1f%%), %d species retained\n\n",
            sum(all_window_species$n_sites_pos < 15),
            100 * mean(all_window_species$n_sites_pos < 15),
            n_distinct(filter(all_window_species, n_sites_pos >= 15)$species)))


# ── 2. Strict thresholds (post-hoc filter) ──────────────────────────────────

cat("── Step 2: Strict thresholds ──────────────────────────────────\n\n")

# min_sites_pos = 10
sens_strict10 <- all_window_species |>
  filter(n_sites_pos >= 10) |>
  prep_species_data()

cat(sprintf("  min_sites_pos=10: %d rows, %d species, %d datasets\n",
            nrow(sens_strict10), n_distinct(sens_strict10$species),
            n_distinct(sens_strict10$dataset)))

mod_strict10 <- fit_m6_lambda(sens_strict10, "(strict, min_sites_pos=10)")

# min_sites_pos = 15
sens_strict15 <- all_window_species |>
  filter(n_sites_pos >= 15) |>
  prep_species_data()

n_spp_15 <- n_distinct(sens_strict15$species)
cat(sprintf("  min_sites_pos=15: %d rows, %d species, %d datasets\n",
            nrow(sens_strict15), n_spp_15,
            n_distinct(sens_strict15$dataset)))

mod_strict15 <- NULL
if (n_spp_15 >= 10) {
  mod_strict15 <- fit_m6_lambda(sens_strict15, "(v.strict, min_sites_pos=15)")
} else {
  cat("  Skipping M6 fit: too few species (<10) at min_sites_pos=15.\n")
}


# ── 3. Load baseline ────────────────────────────────────────────────────────

cat("\n── Step 3: Loading baseline (min_sites_pos = 5) ──────────────────\n\n")

sens_baseline <- readRDS("sensitivity_species_data.rds")
all_models    <- readRDS("sensitivity_gam_models.rds")
mod_baseline  <- all_models$detection$abs_d_lambda

cat(sprintf("  %d rows, %d species (dev.expl = %.1f%%)\n\n",
            nrow(sens_baseline), n_distinct(sens_baseline$species),
            100 * summary(mod_baseline)$dev.expl))


# ── 4. Compare results ──────────────────────────────────────────────────────

cat("── Step 4: Comparison ─────────────────────────────────────────\n\n")

# 4a. Model-level summary
labels   <- c("Baseline (5)", "Strict (10)")
models   <- list(mod_baseline, mod_strict10)
datasets <- list(sens_baseline, sens_strict10)

if (!is.null(mod_strict15)) {
  labels   <- c(labels, "V.Strict (15)")
  models   <- c(models, list(mod_strict15))
  datasets <- c(datasets, list(sens_strict15))
}

model_comp <- tibble(
  threshold      = labels,
  min_sites_pos  = as.integer(str_extract(labels, "\\d+")),
  n_species      = map_int(datasets, ~ n_distinct(.x$species)),
  n_rows         = map_int(datasets, nrow),
  dev_expl       = round(map_dbl(models, ~ summary(.x)$dev.expl), 3),
  r_sq           = round(map_dbl(models, ~ summary(.x)$r.sq), 3)
)

cat("Model-level comparison:\n")
print(model_comp)

# 4b. Parametric coefficients
extract_parametric <- function(mod, label) {
  s <- summary(mod)
  ptable <- as.data.frame(s$p.table)
  ptable$term <- rownames(ptable)
  ptable$threshold <- label
  as_tibble(ptable)
}

params <- map2_dfr(models, labels, extract_parametric) |>
  filter(term %in% c("s_bio4", "l_trapdays", "l_nsites", "s_latitude", "s_trap_array")) |>
  select(threshold, term, Estimate, `Std. Error`, `Pr(>|t|)`) |>
  arrange(term, threshold)

cat("\nParametric coefficients:\n")
params |>
  mutate(across(c(Estimate, `Std. Error`), ~ round(.x, 4)),
         `Pr(>|t|)` = formatC(`Pr(>|t|)`, format = "g", digits = 3)) |>
  print(n = 30)

# 4c. Predictions on shared species
shared_species <- Reduce(intersect, map(datasets, ~ unique(.x$species)))

cat(sprintf("\nShared species across all thresholds: %d\n", length(shared_species)))

make_pred_grid <- function(sens_data, shared_spp) {
  med_covs <- list(
    s_bio4       = 0,
    l_trapdays   = median(sens_data$l_trapdays),
    l_nsites     = median(sens_data$l_nsites),
    s_latitude   = 0,
    s_trap_array = 0,
    ds_sp_f      = sens_data$ds_sp_f[1]
  )
  expand.grid(
    day_start  = seq(1, 358, by = 7),
    window_len = seq(15, 120, by = 1),
    species_f  = factor(shared_spp, levels = levels(sens_data$species_f)),
    stringsAsFactors = FALSE
  ) |>
    as_tibble() |>
    bind_cols(as_tibble(med_covs))
}

grids <- map(datasets, ~ make_pred_grid(.x, shared_species))
preds <- map2(models, grids, ~ predict(.x, newdata = .y, type = "response",
                                        exclude = "s(ds_sp_f)"))
for (i in seq_along(grids)) grids[[i]]$pred <- preds[[i]]

# Duration curves
duration_comp <- map2_dfr(grids, labels, ~ {
  .x |>
    summarise(mean_pred = mean(pred), .by = window_len) |>
    mutate(threshold = .y)
})

cat("\nDuration curve key values (mean |d_lambda|):\n")
duration_comp |>
  filter(window_len %in% c(15, 30, 60, 90, 120)) |>
  pivot_wider(names_from = threshold, values_from = mean_pred) |>
  mutate(across(-window_len, ~ round(.x, 5))) |>
  print()

# Seasonal profiles (60-day windows)
seasonal_comp <- map2_dfr(grids, labels, ~ {
  .x |>
    filter(window_len == 60) |>
    summarise(mean_pred = mean(pred), .by = day_start) |>
    mutate(threshold = .y)
})

# Surface correlations
wide_surface <- map2_dfr(grids, labels, ~ {
  .x |>
    select(day_start, window_len, species_f, pred) |>
    mutate(threshold = .y)
}) |>
  pivot_wider(names_from = threshold, values_from = pred)

wide_season <- seasonal_comp |>
  pivot_wider(names_from = threshold, values_from = mean_pred)

cat("\nSurface correlations vs baseline:\n")
for (lbl in setdiff(labels, "Baseline (5)")) {
  if (lbl %in% names(wide_surface)) {
    cat(sprintf("  Full 2D surface — %s: %.4f\n", lbl,
                cor(wide_surface[[lbl]], wide_surface[["Baseline (5)"]])))
  }
  if (lbl %in% names(wide_season)) {
    cat(sprintf("  Seasonal (60d)  — %s: %.4f\n", lbl,
                cor(wide_season[[lbl]], wide_season[["Baseline (5)"]])))
  }
}

# Add correlations to summary table
model_comp$surface_cor_vs_baseline <- map_dbl(labels, ~ {
  if (.x == "Baseline (5)") return(1.0)
  if (.x %in% names(wide_surface)) {
    round(cor(wide_surface[[.x]], wide_surface[["Baseline (5)"]]), 4)
  } else NA_real_
})
model_comp$seasonal_60d_cor_vs_baseline <- map_dbl(labels, ~ {
  if (.x == "Baseline (5)") return(1.0)
  if (.x %in% names(wide_season)) {
    round(cor(wide_season[[.x]], wide_season[["Baseline (5)"]]), 4)
  } else NA_real_
})


# ── 5. Figures ──────────────────────────────────────────────────────────────

threshold_colors <- c("Baseline (5)" = "#1B1B1B",
                      "Strict (10)"  = "#B2182B",
                      "V.Strict (15)" = "#E08214")
threshold_colors <- threshold_colors[labels]

p_duration <- ggplot(duration_comp, aes(window_len, mean_pred, color = threshold)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = threshold_colors) +
  labs(x = "Window duration (days)", y = "Mean |d_lambda|",
       color = "min_sites_pos", title = "Duration curves") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p_seasonal <- ggplot(seasonal_comp, aes(day_start, mean_pred, color = threshold)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = threshold_colors) +
  scale_x_continuous(breaks = c(1, 91, 182, 274),
                     labels = c("Jan", "Apr", "Jul", "Oct")) +
  labs(x = "Window start", y = "Mean |d_lambda|",
       color = "min_sites_pos", title = "Seasonal profiles (60-day windows)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

fig_threshold <- p_duration + p_seasonal +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dir.create("figures", showWarnings = FALSE)
ggsave("figures/threshold_sensitivity_sites_comparison.pdf",
       fig_threshold, width = 9, height = 4.5)
cat("\nSaved figures/threshold_sensitivity_sites_comparison.pdf\n")


# ── 6. Save summary ────────────────────────────────────────────────────────

write.csv(model_comp, "threshold_sensitivity_sites_summary.csv", row.names = FALSE)
cat("Saved threshold_sensitivity_sites_summary.csv\n\n")

cat("── Final summary ──────────────────────────────────────────────\n")
print(model_comp)
cat("\nDone.\n")

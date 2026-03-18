# threshold_sensitivity.R
# ============================================================================
# Sensitivity analysis: does the min_events species-inclusion threshold
# (10 vs 20 vs 30) affect the main results?
#
# Tests three scenarios:
#   Lenient  (min_events = 10): re-runs the full pipeline to include
#                                previously-dropped species × windows
#   Baseline (min_events = 20): existing data
#   Strict   (min_events = 30): subset of existing data
#
# For each scenario, fits the M6 species-specific GAM surface for |d_lambda|
# and compares: deviance explained, parametric coefficients, duration curves,
# seasonal profiles, and full 2D surface correlations.
#
# Prerequisites:
#   - helpers.R and Full1.R sourced (for dataset_wrapper1, build_window_metrics_fast1)
#   - anchors, ds_paths, spp_keep available in the session
#   - furrr/future plan set up for parallel pipeline re-run
#   - all_window_species.rds, all_dropped_species.rds on disk
#   - sensitivity_species_data.rds on disk (baseline)
#   - sensitivity_gam_models.rds on disk (baseline M6)
#   - species_meta.csv, dataset_metadata.csv on disk
#
# Output:
#   threshold_sensitivity_summary.csv
#   figures/threshold_sensitivity_comparison.pdf
# ============================================================================

library(tidyverse)
library(mgcv)
library(furrr)
library(patchwork)

N_THREADS <- 4
cc_knots  <- list(day_start = c(0, 365))

cat("=== Threshold Sensitivity Analysis ===\n\n")

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

# Rows at the margins
cat("Retained data — n_events_total breakdown:\n")
all_window_species |>
  mutate(bin = cut(n_events_total, breaks = c(0, 10, 20, 30, 50, Inf),
                   labels = c("1-10", "11-20", "21-30", "31-50", "50+"),
                   right = TRUE)) |>
  count(bin) |> mutate(pct = round(100 * n / sum(n), 1)) |> print()

cat("\nDropped data — n_events_total breakdown:\n")
all_dropped_species |>
  mutate(bin = cut(n_events_total, breaks = c(0, 5, 10, 20, 30, Inf),
                   labels = c("1-5", "6-10", "11-20", "21-30", "30+"),
                   right = TRUE)) |>
  count(bin) |> mutate(pct = round(100 * n / sum(n), 1)) |> print()

# Potential gains/losses
spp_at_20 <- unique(all_window_species$species)
spp_at_30 <- unique(filter(all_window_species, n_events_total >= 30)$species)
potential_new <- all_dropped_species |>
  filter(n_events_total >= 10, !failed_min_occasions_pos, !failed_min_sites_pos) |>
  pull(species) |> unique()
spp_only_new <- setdiff(potential_new, spp_at_20)

cat(sprintf("\nStrict (30): lose %d rows, same %d species (none dropped entirely)\n",
            sum(all_window_species$n_events_total < 30),
            length(spp_at_30)))
cat(sprintf("Lenient (10): potential %d new species: %s\n\n",
            length(spp_only_new), paste(spp_only_new, collapse = ", ")))


# ── 2. Strict threshold (min_events = 30) ───────────────────────────────────

cat("── Step 2: Strict threshold (min_events = 30) ────────────────\n\n")

sens_strict <- all_window_species |>
  filter(n_events_total >= 30) |>
  prep_species_data()

cat(sprintf("  %d rows, %d species, %d datasets\n",
            nrow(sens_strict), n_distinct(sens_strict$species),
            n_distinct(sens_strict$dataset)))

mod_strict <- fit_m6_lambda(sens_strict, "(strict, min_events=30)")


# ── 3. Lenient threshold (min_events = 10) ──────────────────────────────────

cat("\n── Step 3: Lenient threshold (min_events = 10) ───────────────\n\n")
cat("  Re-running pipeline across", length(ds_paths), "datasets...\n")

plan(multisession, workers = 4)

t0 <- proc.time()
wrapped_me10 <- furrr::future_map(
  ds_paths,
  dataset_wrapper1,
  min_events = 10,
  .options = furrr::furrr_options(seed = TRUE)
)
t_pipeline <- round((proc.time() - t0)[3] / 60, 1)
cat(sprintf("  Pipeline done in %.1f minutes.\n", t_pipeline))

aws_me10 <- purrr::map_df(wrapped_me10, "window_species")

sens_lenient <- aws_me10 |> prep_species_data(min_rows_per_species = 50)

cat(sprintf("  %d rows, %d species, %d datasets\n",
            nrow(sens_lenient), n_distinct(sens_lenient$species),
            n_distinct(sens_lenient$dataset)))

mod_lenient <- fit_m6_lambda(sens_lenient, "(lenient, min_events=10)")


# ── 4. Load baseline ────────────────────────────────────────────────────────

cat("\n── Loading baseline (min_events = 20) ──────────────────────────\n\n")

sens_baseline <- readRDS("sensitivity_species_data.rds")
all_models    <- readRDS("sensitivity_gam_models.rds")
mod_baseline  <- all_models$detection$abs_d_lambda

cat(sprintf("  %d rows, %d species (dev.expl = %.1f%%)\n\n",
            nrow(sens_baseline), n_distinct(sens_baseline$species),
            100 * summary(mod_baseline)$dev.expl))


# ── 5. Compare results ─────────────────────────────────────────────────────

cat("── Step 4: Comparison ─────────────────────────────────────────\n\n")

# 5a. Model-level summary
model_comp <- tibble(
  threshold   = c("Lenient (10)", "Baseline (20)", "Strict (30)"),
  min_events  = c(10, 20, 30),
  n_species   = c(n_distinct(sens_lenient$species),
                  n_distinct(sens_baseline$species),
                  n_distinct(sens_strict$species)),
  n_rows      = c(nrow(sens_lenient), nrow(sens_baseline), nrow(sens_strict)),
  dev_expl    = round(c(summary(mod_lenient)$dev.expl,
                        summary(mod_baseline)$dev.expl,
                        summary(mod_strict)$dev.expl), 3),
  r_sq        = round(c(summary(mod_lenient)$r.sq,
                        summary(mod_baseline)$r.sq,
                        summary(mod_strict)$r.sq), 3)
)

cat("Model-level comparison:\n")
print(model_comp)

# 5b. Parametric coefficients
extract_parametric <- function(mod, label) {
  s <- summary(mod)
  ptable <- as.data.frame(s$p.table)
  ptable$term <- rownames(ptable)
  ptable$threshold <- label
  as_tibble(ptable)
}

params <- bind_rows(
  extract_parametric(mod_lenient, "Lenient (10)"),
  extract_parametric(mod_baseline, "Baseline (20)"),
  extract_parametric(mod_strict, "Strict (30)")
) |>
  filter(term %in% c("s_bio4", "l_trapdays", "l_nsites", "s_latitude", "s_trap_array")) |>
  select(threshold, term, Estimate, `Std. Error`, `Pr(>|t|)`) |>
  arrange(term, threshold)

cat("\nParametric coefficients:\n")
params |>
  mutate(across(c(Estimate, `Std. Error`), ~ round(.x, 4)),
         `Pr(>|t|)` = formatC(`Pr(>|t|)`, format = "g", digits = 3)) |>
  print(n = 20)

# 5c. Predictions on shared species
shared_species <- Reduce(intersect, list(
  unique(sens_lenient$species),
  unique(sens_baseline$species),
  unique(sens_strict$species)
))

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

grid_lenient  <- make_pred_grid(sens_lenient, shared_species)
grid_baseline <- make_pred_grid(sens_baseline, shared_species)
grid_strict   <- make_pred_grid(sens_strict, shared_species)

grid_lenient$pred  <- predict(mod_lenient,  newdata = grid_lenient,  type = "response", exclude = "s(ds_sp_f)")
grid_baseline$pred <- predict(mod_baseline, newdata = grid_baseline, type = "response", exclude = "s(ds_sp_f)")
grid_strict$pred   <- predict(mod_strict,   newdata = grid_strict,   type = "response", exclude = "s(ds_sp_f)")

# Duration curves
duration_comp <- bind_rows(
  grid_lenient  |> summarise(mean_pred = mean(pred), .by = window_len) |> mutate(threshold = "Lenient (10)"),
  grid_baseline |> summarise(mean_pred = mean(pred), .by = window_len) |> mutate(threshold = "Baseline (20)"),
  grid_strict   |> summarise(mean_pred = mean(pred), .by = window_len) |> mutate(threshold = "Strict (30)")
)

cat("\nDuration curve key values (mean |d_lambda|):\n")
duration_comp |>
  filter(window_len %in% c(15, 30, 60, 90, 120)) |>
  pivot_wider(names_from = threshold, values_from = mean_pred) |>
  mutate(across(-window_len, ~ round(.x, 5))) |>
  print()

# Surface correlations
wide_surface <- bind_rows(
  grid_lenient  |> select(day_start, window_len, species_f, pred) |> mutate(threshold = "lenient"),
  grid_baseline |> select(day_start, window_len, species_f, pred) |> mutate(threshold = "baseline"),
  grid_strict   |> select(day_start, window_len, species_f, pred) |> mutate(threshold = "strict")
) |>
  pivot_wider(names_from = threshold, values_from = pred)

seasonal_comp <- bind_rows(
  grid_lenient  |> filter(window_len == 60) |> summarise(mean_pred = mean(pred), .by = day_start) |> mutate(threshold = "Lenient (10)"),
  grid_baseline |> filter(window_len == 60) |> summarise(mean_pred = mean(pred), .by = day_start) |> mutate(threshold = "Baseline (20)"),
  grid_strict   |> filter(window_len == 60) |> summarise(mean_pred = mean(pred), .by = day_start) |> mutate(threshold = "Strict (30)")
)

wide_season <- seasonal_comp |>
  pivot_wider(names_from = threshold, values_from = mean_pred)

cat("\nSurface correlations vs baseline:\n")
cat("  Full 2D surface — Lenient:", round(cor(wide_surface$lenient, wide_surface$baseline), 4), "\n")
cat("  Full 2D surface — Strict: ", round(cor(wide_surface$strict,  wide_surface$baseline), 4), "\n")
cat("  Seasonal (60d)  — Lenient:", round(cor(wide_season$`Lenient (10)`, wide_season$`Baseline (20)`), 4), "\n")
cat("  Seasonal (60d)  — Strict: ", round(cor(wide_season$`Strict (30)`,  wide_season$`Baseline (20)`), 4), "\n")

# Add correlations to summary table
model_comp$surface_cor_vs_baseline <- c(
  round(cor(wide_surface$lenient, wide_surface$baseline), 4),
  1.0,
  round(cor(wide_surface$strict, wide_surface$baseline), 4)
)
model_comp$seasonal_60d_cor_vs_baseline <- c(
  round(cor(wide_season$`Lenient (10)`, wide_season$`Baseline (20)`), 4),
  1.0,
  round(cor(wide_season$`Strict (30)`, wide_season$`Baseline (20)`), 4)
)


# ── 6. Figures ──────────────────────────────────────────────────────────────

p_duration <- ggplot(duration_comp, aes(window_len, mean_pred, color = threshold)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("Lenient (10)" = "#2166AC",
                                "Baseline (20)" = "#1B1B1B",
                                "Strict (30)" = "#B2182B")) +
  labs(x = "Window duration (days)", y = "Mean |d_lambda|",
       color = "min_events", title = "Duration curves") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p_seasonal <- ggplot(seasonal_comp, aes(day_start, mean_pred, color = threshold)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("Lenient (10)" = "#2166AC",
                                "Baseline (20)" = "#1B1B1B",
                                "Strict (30)" = "#B2182B")) +
  scale_x_continuous(breaks = c(1, 91, 182, 274),
                     labels = c("Jan", "Apr", "Jul", "Oct")) +
  labs(x = "Window start", y = "Mean |d_lambda|",
       color = "min_events", title = "Seasonal profiles (60-day windows)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

fig_threshold <- p_duration + p_seasonal +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dir.create("figures", showWarnings = FALSE)
ggsave("figures/threshold_sensitivity_comparison.pdf",
       fig_threshold, width = 9, height = 4.5)
cat("\nSaved figures/threshold_sensitivity_comparison.pdf\n")


# ── 7. Save summary ────────────────────────────────────────────────────────

write.csv(model_comp, "threshold_sensitivity_summary.csv", row.names = FALSE)
cat("Saved threshold_sensitivity_summary.csv\n\n")

cat("── Final summary ──────────────────────────────────────────────\n")
print(model_comp)
cat("\nDone.\n")

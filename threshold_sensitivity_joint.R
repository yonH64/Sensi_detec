# threshold_sensitivity_joint.R
# ============================================================================
# Joint threshold sensitivity analysis: do the species-inclusion thresholds
# (min_events, min_sites_pos, min_occasions_pos) affect the main results?
#
# Tests three joint configurations:
#   Lenient  (min_events=10, min_sites_pos=3, min_occasions_pos=3)
#   Baseline (min_events=20, min_sites_pos=5, min_occasions_pos=5)
#   Strict   (min_events=30, min_sites_pos=10, min_occasions_pos=5)
#
# Strategy: run the full pipeline once at the lenient config to capture all
# possible species × window combinations. Derive baseline and strict via
# post-hoc filtering on n_events_total, n_sites_pos, n_occasions_pos columns
# (which are retained in all_window_species).
#
# Prerequisites:
#   - helpers.R and Full1.R sourced
#   - anchors, ds_paths, spp_keep, window_grid available in session
#   - sensitivity_species_data.rds on disk (baseline reference)
#   - sensitivity_gam_models.rds on disk (baseline M6)
#   - species_meta.csv, dataset_metadata.csv on disk
#
# Output:
#   threshold_sensitivity_joint_summary.csv
#   figures/threshold_sensitivity_joint.pdf
# ============================================================================

library(tidyverse)
library(mgcv)
library(furrr)
library(patchwork)

N_THREADS <- 4
cc_knots  <- list(day_start = c(0, 365))

cat("=== Joint Threshold Sensitivity Analysis ===\n\n")

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
    left_join(meta |> select(dataset_name, all_of(join_cols)) |>
                rename(dataset = dataset_name), by = "dataset")
  lookup <- meta |> select(dataset_name, all_of(join_cols)) |>
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
      abs_d_rate         = pmax(abs(d_rate), 1e-10),
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


# ── 1. Re-run pipeline at lenient thresholds ─────────────────────────────────

cat("── Step 1: Lenient pipeline re-run (10/3/3) ─────────────────────\n\n")
cat("  Re-running pipeline across", length(ds_paths), "datasets...\n")

plan(multisession, workers = 4)

t0 <- proc.time()
wrapped_lenient <- furrr::future_map(
  ds_paths,
  dataset_wrapper1,
  min_events        = 10,
  min_sites_pos     = 3,
  min_occasions_pos = 3,
  .options = furrr::furrr_options(seed = TRUE)
)
t_pipeline <- round((proc.time() - t0)[3] / 60, 1)
cat(sprintf("  Pipeline done in %.1f minutes.\n", t_pipeline))

aws_lenient <- purrr::map_df(wrapped_lenient, "window_species")

cat(sprintf("  Lenient raw: %d rows, %d species\n",
            nrow(aws_lenient), n_distinct(aws_lenient$species)))


# ── 2. Derive three configs via post-hoc filtering ───────────────────────────

cat("\n── Step 2: Post-hoc filtering ──────────────────────────────────\n\n")

# Lenient (10/3/3): already have
sens_lenient <- aws_lenient |> prep_species_data()
cat(sprintf("  Lenient (10/3/3): %d rows, %d species, %d datasets\n",
            nrow(sens_lenient), n_distinct(sens_lenient$species),
            n_distinct(sens_lenient$dataset)))

# Baseline (20/5/5): filter from lenient
aws_baseline <- aws_lenient |>
  filter(n_events_total >= 20, n_sites_pos >= 5, n_occasions_pos >= 5)
sens_baseline <- aws_baseline |> prep_species_data()
cat(sprintf("  Baseline (20/5/5): %d rows, %d species, %d datasets\n",
            nrow(sens_baseline), n_distinct(sens_baseline$species),
            n_distinct(sens_baseline$dataset)))

# Strict (30/10/5): filter from lenient
aws_strict <- aws_lenient |>
  filter(n_events_total >= 30, n_sites_pos >= 10, n_occasions_pos >= 5)
sens_strict <- aws_strict |> prep_species_data()
cat(sprintf("  Strict (30/10/5): %d rows, %d species, %d datasets\n",
            nrow(sens_strict), n_distinct(sens_strict$species),
            n_distinct(sens_strict$dataset)))


# ── 3. Load reference baseline (from current pipeline) ──────────────────────

cat("\n── Step 3: Reference baseline comparison ──────────────────────\n\n")

sens_reference <- readRDS("sensitivity_species_data.rds")
all_models     <- readRDS("sensitivity_gam_models.rds")
mod_reference  <- all_models$detection$abs_d_lambda

cat(sprintf("  Reference (current pipeline): %d rows, %d species (dev.expl=%.1f%%)\n\n",
            nrow(sens_reference), n_distinct(sens_reference$species),
            100 * summary(mod_reference)$dev.expl))

# Check that baseline from filter matches reference from pipeline
cat(sprintf("  Baseline vs Reference row count: %d vs %d (delta=%d)\n",
            nrow(sens_baseline), nrow(sens_reference),
            nrow(sens_baseline) - nrow(sens_reference)))
cat(sprintf("  Baseline vs Reference species: %d vs %d\n",
            n_distinct(sens_baseline$species), n_distinct(sens_reference$species)))


# ── 4. Fit M6 for each config ───────────────────────────────────────────────

cat("\n── Step 4: Model fitting ──────────────────────────────────────\n\n")

mod_lenient  <- fit_m6_lambda(sens_lenient,  "(lenient, 10/3/3)")
mod_baseline <- fit_m6_lambda(sens_baseline, "(baseline, 20/5/5)")
mod_strict   <- fit_m6_lambda(sens_strict,   "(strict, 30/10/5)")


# ── 5. Compare results ──────────────────────────────────────────────────────

cat("\n── Step 5: Comparison ─────────────────────────────────────────\n\n")

# 5a. Model-level summary
model_comp <- tibble(
  config      = c("Lenient (10/3/3)", "Baseline (20/5/5)", "Strict (30/10/5)"),
  min_events  = c(10, 20, 30),
  min_sites   = c(3, 5, 10),
  min_occ     = c(3, 5, 5),
  n_species   = c(n_distinct(sens_lenient$species),
                  n_distinct(sens_baseline$species),
                  n_distinct(sens_strict$species)),
  n_rows      = c(nrow(sens_lenient), nrow(sens_baseline), nrow(sens_strict)),
  dev_expl    = round(c(summary(mod_lenient)$dev.expl,
                        summary(mod_baseline)$dev.expl,
                        summary(mod_strict)$dev.expl), 3),
  r_sq        = round(c(summary(mod_lenient)$r.sq,
                        summary(mod_baseline)$r.sq,
                        summary(mod_strict)$r.sq), 3),
  AIC         = c(round(AIC(mod_lenient)), round(AIC(mod_baseline)),
                  round(AIC(mod_strict)))
)

cat("Model-level comparison:\n")
print(model_comp)

# 5b. Parametric coefficients
extract_parametric <- function(mod, label) {
  s <- summary(mod)
  ptable <- as.data.frame(s$p.table)
  ptable$term <- rownames(ptable)
  ptable$config <- label
  as_tibble(ptable)
}

params <- bind_rows(
  extract_parametric(mod_lenient, "Lenient (10/3/3)"),
  extract_parametric(mod_baseline, "Baseline (20/5/5)"),
  extract_parametric(mod_strict, "Strict (30/10/5)")
) |>
  filter(term %in% c("s_bio4", "l_trapdays", "l_nsites", "s_latitude", "s_trap_array")) |>
  select(config, term, Estimate, `Std. Error`, `Pr(>|t|)`) |>
  arrange(term, config)

cat("\nParametric coefficients:\n")
params |>
  mutate(across(c(Estimate, `Std. Error`), ~ round(.x, 4)),
         `Pr(>|t|)` = formatC(`Pr(>|t|)`, format = "g", digits = 3)) |>
  print(n = 20)

# 5c. Predictions on shared species → surface correlations
shared_species <- Reduce(intersect, list(
  unique(sens_lenient$species),
  unique(sens_baseline$species),
  unique(sens_strict$species)
))

cat(sprintf("\n  Shared species across all configs: %d\n", length(shared_species)))

make_pred_grid <- function(sens_data, shared_spp) {
  expand.grid(
    day_start  = seq(1, 358, by = 7),
    window_len = seq(15, 183, by = 1),
    species_f  = factor(shared_spp, levels = levels(sens_data$species_f)),
    stringsAsFactors = FALSE
  ) |>
    as_tibble() |>
    mutate(
      s_bio4       = 0,
      l_trapdays   = median(sens_data$l_trapdays),
      l_nsites     = median(sens_data$l_nsites),
      s_latitude   = 0,
      s_trap_array = 0,
      ds_sp_f      = sens_data$ds_sp_f[1]
    )
}

grid_lenient  <- make_pred_grid(sens_lenient, shared_species)
grid_baseline <- make_pred_grid(sens_baseline, shared_species)
grid_strict   <- make_pred_grid(sens_strict, shared_species)

grid_lenient$pred  <- predict(mod_lenient,  newdata = grid_lenient,  type = "response", exclude = "s(ds_sp_f)")
grid_baseline$pred <- predict(mod_baseline, newdata = grid_baseline, type = "response", exclude = "s(ds_sp_f)")
grid_strict$pred   <- predict(mod_strict,   newdata = grid_strict,   type = "response", exclude = "s(ds_sp_f)")

# Duration curves
duration_comp <- bind_rows(
  grid_lenient  |> summarise(mean_pred = mean(pred), .by = window_len) |> mutate(config = "Lenient (10/3/3)"),
  grid_baseline |> summarise(mean_pred = mean(pred), .by = window_len) |> mutate(config = "Baseline (20/5/5)"),
  grid_strict   |> summarise(mean_pred = mean(pred), .by = window_len) |> mutate(config = "Strict (30/10/5)")
)

cat("\nDuration curve key values (mean |d_lambda|):\n")
duration_comp |>
  filter(window_len %in% c(15, 30, 60, 90, 120, 183)) |>
  pivot_wider(names_from = config, values_from = mean_pred) |>
  mutate(across(-window_len, ~ round(.x, 5))) |>
  print()

# Surface correlations
wide_surface <- bind_rows(
  grid_lenient  |> select(day_start, window_len, species_f, pred) |> mutate(config = "lenient"),
  grid_baseline |> select(day_start, window_len, species_f, pred) |> mutate(config = "baseline"),
  grid_strict   |> select(day_start, window_len, species_f, pred) |> mutate(config = "strict")
) |>
  pivot_wider(names_from = config, values_from = pred)

seasonal_comp <- bind_rows(
  grid_lenient  |> filter(window_len == 60) |> summarise(mean_pred = mean(pred), .by = day_start) |> mutate(config = "Lenient (10/3/3)"),
  grid_baseline |> filter(window_len == 60) |> summarise(mean_pred = mean(pred), .by = day_start) |> mutate(config = "Baseline (20/5/5)"),
  grid_strict   |> filter(window_len == 60) |> summarise(mean_pred = mean(pred), .by = day_start) |> mutate(config = "Strict (30/10/5)")
)

wide_season <- seasonal_comp |>
  pivot_wider(names_from = config, values_from = mean_pred)

cat("\nSurface correlations vs baseline:\n")
cat("  Full 2D surface — Lenient vs Baseline:", round(cor(wide_surface$lenient, wide_surface$baseline), 4), "\n")
cat("  Full 2D surface — Strict vs Baseline: ", round(cor(wide_surface$strict,  wide_surface$baseline), 4), "\n")
cat("  Seasonal (60d)  — Lenient vs Baseline:", round(cor(wide_season$`Lenient (10/3/3)`, wide_season$`Baseline (20/5/5)`), 4), "\n")
cat("  Seasonal (60d)  — Strict vs Baseline: ", round(cor(wide_season$`Strict (30/10/5)`,  wide_season$`Baseline (20/5/5)`), 4), "\n")

model_comp$surface_cor_vs_baseline <- c(
  round(cor(wide_surface$lenient, wide_surface$baseline), 4),
  1.0,
  round(cor(wide_surface$strict, wide_surface$baseline), 4)
)
model_comp$seasonal_60d_cor_vs_baseline <- c(
  round(cor(wide_season$`Lenient (10/3/3)`, wide_season$`Baseline (20/5/5)`), 4),
  1.0,
  round(cor(wide_season$`Strict (30/10/5)`, wide_season$`Baseline (20/5/5)`), 4)
)


# ── 6. Figures ──────────────────────────────────────────────────────────────

p_duration <- ggplot(duration_comp, aes(window_len, mean_pred, color = config)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("Lenient (10/3/3)" = "#2166AC",
                                "Baseline (20/5/5)" = "#1B1B1B",
                                "Strict (30/10/5)" = "#B2182B")) +
  labs(x = "Window duration (days)", y = "Mean |d_lambda|",
       color = "Threshold config", title = "Duration curves") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p_seasonal <- ggplot(seasonal_comp, aes(day_start, mean_pred, color = config)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("Lenient (10/3/3)" = "#2166AC",
                                "Baseline (20/5/5)" = "#1B1B1B",
                                "Strict (30/10/5)" = "#B2182B")) +
  scale_x_continuous(breaks = c(1, 91, 182, 274),
                     labels = c("Jan", "Apr", "Jul", "Oct")) +
  labs(x = "Window start", y = "Mean |d_lambda|",
       color = "Threshold config", title = "Seasonal profiles (60-day windows)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

fig <- p_duration + p_seasonal +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dir.create("figures", showWarnings = FALSE)
ggsave("figures/threshold_sensitivity_joint.pdf", fig, width = 9, height = 4.5)
cat("\nSaved figures/threshold_sensitivity_joint.pdf\n")


# ── 7. Save summary ────────────────────────────────────────────────────────

write.csv(model_comp, "threshold_sensitivity_joint_summary.csv", row.names = FALSE)
cat("Saved threshold_sensitivity_joint_summary.csv\n\n")

cat("── Final summary ──────────────────────────────────────────────\n")
print(model_comp)
cat("\nDone.\n")

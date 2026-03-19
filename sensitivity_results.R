# sensitivity_results.R
# ============================================================================
# Extract derived quantities from the sensitivity surface models to answer
# research questions Q1–Q9. Outputs CSV summary tables and key objects.
#
# Input:  sensitivity_models_env.RData  (from models_sensitivity_surface.R)
#         interannual_cv_lambda_full.csv (for Q9 noise floor)
#
# Output: Q1_duration_effect.csv
#         Q2_seasonal_profiles.csv
#         Q3_surface_predictions.csv
#         Q4_species_guild_surfaces.csv
#         Q5_bio4_effect.csv
#         Q6_snapshot_europe.csv
#         Q7_optimal_timing.csv
#         Q8_richness_surface.csv
#         Q9_noise_floor.csv
#         model_comparison_table.csv
# ============================================================================

library(mgcv)
library(dplyr)
library(tidyr)
library(tibble)

cat("=== Sensitivity Results Extraction ===\n\n")

# ── Load ─────────────────────────────────────────────────────────────────────

load("sensitivity_models_env.RData")

mod_lambda      <- all_models$detection$abs_d_lambda
mod_rate        <- all_models$detection$abs_d_rate
mod_matched     <- all_models$detection$abs_d_matched_rate
mod_rate_signed <- all_models$detection$d_rate_signed
mod_richness    <- all_models$community$d_sr_raref
mod_prop        <- all_models$community$prop_sr_full
mod_rho         <- all_models$community$rho_lambda

species_levels <- levels(sens_species$species_f)
guild_info     <- sens_species |>
  distinct(species_f, guild_major, guild_minor_habitat, guild_minor_diet)


# ── Shared prediction scaffolds ─────────────────────────────────────────────

# Median covariate values for "typical" predictions
median_covs <- list(
  s_bio4       = 0,
  l_trapdays   = median(sens_species$l_trapdays),
  l_nsites     = median(sens_species$l_nsites),
  s_latitude   = 0,
  s_trap_array = 0,
  ds_sp_f      = sens_species$ds_sp_f[1]
)

# Species-level grid (fine resolution)
sp_grid <- expand.grid(
  day_start  = seq(1, 358, by = 7),
  window_len = seq(15, 183, by = 1),
  species_f  = factor(species_levels, levels = species_levels),
  stringsAsFactors = FALSE
) |>
  as_tibble() |>
  bind_cols(as_tibble(median_covs))

# Dataset-level grid for richness models
rich_grid_base <- expand.grid(
  day_start  = seq(1, 358, by = 7),
  window_len = seq(15, 183, by = 1)
) |>
  as_tibble() |>
  mutate(
    s_bio4     = 0,
    l_trapdays = median(sens_richness$l_trapdays),
    l_nsites   = median(sens_richness$l_nsites),
    s_latitude = 0,
    dataset_f  = sens_richness$dataset_f[1]
  )


# ── Helper: predict from a detection model per species ───────────────────────

predict_detection <- function(model, grid = sp_grid) {
  grid$pred <- predict(model, newdata = grid, type = "response",
                       exclude = "s(ds_sp_f)")
  grid |> left_join(guild_info, by = "species_f")
}

predict_richness <- function(model, grid = rich_grid_base) {
  grid$pred <- predict(model, newdata = grid, type = "response",
                       exclude = "s(dataset_f)")
  grid
}


# ============================================================================
# MODEL COMPARISON TABLE
# ============================================================================

cat("── Model comparison table ──────────────────────────────────────\n")

comp_models <- all_models$comparison
comp_tbl <- tibble(
  model    = names(comp_models),
  AIC      = sapply(comp_models, AIC),
  dev_expl = sapply(comp_models, \(m) summary(m)$dev.expl),
  edf      = sapply(comp_models, \(m) round(sum(m$edf), 1)),
  r_sq     = sapply(comp_models, \(m) summary(m)$r.sq)
) |>
  arrange(AIC) |>
  mutate(delta_AIC = round(AIC - min(AIC)),
         AIC = round(AIC),
         dev_expl = round(dev_expl, 3),
         r_sq = round(r_sq, 3))

print(comp_tbl, n = 10)
write.csv(comp_tbl, "model_comparison_table.csv", row.names = FALSE)
cat("\n")


# ============================================================================
# Q1: DURATION EFFECT
# ============================================================================
# How does deviation decrease with window duration?

cat("── Q1: Duration effect ────────────────────────────────────────\n")

preds_lambda <- predict_detection(mod_lambda)

q1 <- preds_lambda |>
  summarise(
    mean_pred = mean(pred),
    median_pred = median(pred),
    q25 = quantile(pred, 0.25),
    q75 = quantile(pred, 0.75),
    .by = c(window_len, guild_major)
  ) |>
  arrange(guild_major, window_len)

# Also compute overall (across all species)
q1_overall <- preds_lambda |>
  summarise(
    mean_pred = mean(pred),
    median_pred = median(pred),
    q25 = quantile(pred, 0.25),
    q75 = quantile(pred, 0.75),
    .by = window_len
  ) |>
  mutate(guild_major = "All") |>
  arrange(window_len)

q1_out <- bind_rows(q1, q1_overall)
write.csv(q1_out, "Q1_duration_effect.csv", row.names = FALSE)
cat("  Saved Q1_duration_effect.csv\n")

# Key insight: how much deviation drops from 15 → 60 → 120 days
q1_overall |>
  filter(window_len %in% c(15, 30, 60, 90, 120, 150, 183)) |>
  mutate(across(c(mean_pred, median_pred), ~round(., 4))) |>
  print()
cat("\n")


# ============================================================================
# Q2: SEASONAL PROFILES
# ============================================================================
# How does deviation vary across the year, at fixed durations?

cat("── Q2: Seasonal profiles ──────────────────────────────────────\n")

q2 <- preds_lambda |>
  filter(window_len %in% c(15, 30, 60, 90, 120, 150, 183)) |>
  summarise(
    mean_pred = mean(pred),
    .by = c(day_start, window_len, guild_major)
  )

q2_overall <- preds_lambda |>
  filter(window_len %in% c(15, 30, 60, 90, 120, 150, 183)) |>
  summarise(mean_pred = mean(pred), .by = c(day_start, window_len)) |>
  mutate(guild_major = "All")

q2_out <- bind_rows(q2, q2_overall)
write.csv(q2_out, "Q2_seasonal_profiles.csv", row.names = FALSE)
cat("  Saved Q2_seasonal_profiles.csv\n")

# Worst season vs best season at 60 days
q2_overall |>
  filter(window_len == 60) |>
  summarise(
    best_start  = day_start[which.min(mean_pred)],
    best_pred   = round(min(mean_pred), 4),
    worst_start = day_start[which.max(mean_pred)],
    worst_pred  = round(max(mean_pred), 4),
    ratio       = round(max(mean_pred) / min(mean_pred), 1)
  ) |> print()
cat("\n")


# ============================================================================
# Q3: FULL SURFACE PREDICTIONS
# ============================================================================
# The 2D surface (timing × duration) for each metric.

cat("── Q3: Full surface predictions ───────────────────────────────\n")

# Coarser grid for file size
surf_grid <- expand.grid(
  day_start  = seq(1, 358, by = 7),
  window_len = seq(15, 183, by = 7),
  species_f  = factor(species_levels, levels = species_levels),
  stringsAsFactors = FALSE
) |>
  as_tibble() |>
  bind_cols(as_tibble(median_covs))

q3 <- list()
for (metric_name in c("abs_d_lambda", "abs_d_rate", "abs_d_matched_rate")) {
  model <- all_models$detection[[metric_name]]
  surf_grid$pred <- predict(model, newdata = surf_grid, type = "response",
                            exclude = "s(ds_sp_f)")
  q3[[metric_name]] <- surf_grid |>
    left_join(guild_info, by = "species_f") |>
    summarise(mean_pred = mean(pred), .by = c(day_start, window_len, guild_major)) |>
    mutate(metric = metric_name)
}

# Signed rate (Gaussian model — predictions can be negative)
surf_grid$pred <- predict(mod_rate_signed, newdata = surf_grid, type = "response",
                          exclude = "s(ds_sp_f)")
q3$d_rate_signed <- surf_grid |>
  left_join(guild_info, by = "species_f") |>
  summarise(mean_pred = mean(pred), .by = c(day_start, window_len, guild_major)) |>
  mutate(metric = "d_rate_signed")

q3_out <- bind_rows(q3)
write.csv(q3_out, "Q3_surface_predictions.csv", row.names = FALSE)
cat("  Saved Q3_surface_predictions.csv (", nrow(q3_out), "rows)\n\n")


# ============================================================================
# Q4: SPECIES/GUILD SURFACE VARIATION
# ============================================================================
# How different are species-specific surfaces?

cat("── Q4: Species/guild surface variation ────────────────────────\n")

preds_sp <- predict_detection(mod_lambda, grid = surf_grid)

# Species-level summary: deviation at key durations × seasons
q4 <- preds_sp |>
  mutate(
    season = case_when(
      day_start >= 1   & day_start <= 78  ~ "Winter",
      day_start >= 85  & day_start <= 162 ~ "Spring",
      day_start >= 169 & day_start <= 246 ~ "Summer",
      day_start >= 253 & day_start <= 337 ~ "Autumn",
      TRUE ~ "Winter"
    )
  ) |>
  summarise(
    mean_pred = mean(pred),
    .by = c(species_f, guild_major, guild_minor_habitat, guild_minor_diet,
            window_len, season)
  )

write.csv(q4, "Q4_species_guild_surfaces.csv", row.names = FALSE)
cat("  Saved Q4_species_guild_surfaces.csv\n")

# Which species show the most timing sensitivity at 60 days?
q4 |>
  filter(window_len == 64) |>
  summarise(
    season_range = max(mean_pred) - min(mean_pred),
    worst_season = season[which.max(mean_pred)],
    .by = c(species_f, guild_major)
  ) |>
  arrange(desc(season_range)) |>
  mutate(season_range = round(season_range, 4)) |>
  head(10) |>
  print()
cat("\n")


# ============================================================================
# Q5: TEMPERATURE SEASONALITY EFFECT
# ============================================================================
# How does BIO4 modulate deviation?

cat("── Q5: Temperature seasonality (BIO4) effect ─────────────────\n")

# Extract parametric coefficient for s_bio4 from each detection model
q5 <- tibble(
  metric = c("abs_d_lambda", "abs_d_rate", "abs_d_matched_rate", "d_rate_signed"),
  model  = list(mod_lambda, mod_rate, mod_matched, mod_rate_signed)
) |>
  rowwise() |>
  mutate(
    coef   = coef(model)["s_bio4"],
    se     = summary(model)$p.table["s_bio4", "Std. Error"],
    t_val  = summary(model)$p.table["s_bio4", "t value"],
    p_val  = summary(model)$p.table["s_bio4", "Pr(>|t|)"]
  ) |>
  ungroup() |>
  select(-model) |>
  mutate(across(c(coef, se, t_val), ~round(., 4)),
         p_val = signif(p_val, 3))

print(q5)
write.csv(q5, "Q5_bio4_effect.csv", row.names = FALSE)
cat("  Saved Q5_bio4_effect.csv\n\n")


# ============================================================================
# Q6: PROTOCOL EVALUATION
# ============================================================================
# Where do the named protocols sit on the sensitivity surface?
# Protocols are evaluated as derived predictions at their exact coordinates.

cat("── Q6: Protocol positions on the sensitivity surface ──────────\n")

# Protocol definitions (from helpers.R::protocol_windows())
#   SNAP_EU_CORE:   Sep 1  (DOY 244), 61 days
#   SNAP_EU_BUFFER: Aug 18 (DOY 230), 89 days
#   EOW_EARLY:      Aug 2  (DOY 214), 60 days  (pre-rut half)
#   EOW_LATE:       Oct 1  (DOY 274), 60 days  (post-rut half)
protocol_coords <- tibble(
  protocol   = c("CORE", "BUFFER", "EOW_EARLY", "EOW_LATE"),
  day_start  = c(244, 230, 214, 274),
  window_len = c(61, 89, 60, 60)
)

protocol_grid <- tidyr::crossing(
  protocol_coords,
  tibble(species_f = factor(species_levels, levels = species_levels))
) |>
  bind_cols(as_tibble(median_covs) |> slice(rep(1, n())))

# Predict for all detection metrics
q6_list <- list()
for (metric_name in c("abs_d_lambda", "abs_d_rate", "abs_d_matched_rate")) {
  model <- all_models$detection[[metric_name]]
  protocol_grid$pred <- predict(model, newdata = protocol_grid, type = "response",
                                exclude = "s(ds_sp_f)")
  q6_list[[metric_name]] <- protocol_grid |>
    left_join(guild_info, by = "species_f") |>
    select(protocol, day_start, window_len, species_f, guild_major, pred) |>
    mutate(metric = metric_name)
}
q6_all <- bind_rows(q6_list)

# Summary per guild × protocol × metric
q6_summary <- q6_all |>
  summarise(mean_pred = mean(pred), .by = c(protocol, guild_major, metric)) |>
  pivot_wider(names_from = protocol, values_from = mean_pred) |>
  mutate(across(c(CORE, BUFFER, EOW_EARLY, EOW_LATE), ~round(., 5)))

print(q6_summary, n = 25)

# Overall summary across guilds (species-averaged)
q6_overall <- q6_all |>
  summarise(mean_pred = mean(pred), .by = c(protocol, metric)) |>
  pivot_wider(names_from = protocol, values_from = mean_pred) |>
  mutate(across(c(CORE, BUFFER, EOW_EARLY, EOW_LATE), ~round(., 5)))

cat("\n  Overall (species-averaged):\n")
print(q6_overall)

# Optimal timing at each protocol's duration
protocol_durations <- unique(protocol_coords$window_len)
cat("\n  Optimal timing vs protocol positions:\n")
for (metric_name in c("abs_d_lambda", "abs_d_rate", "abs_d_matched_rate")) {
  model <- all_models$detection[[metric_name]]
  tmp_grid <- expand.grid(
    day_start = seq(1, 358, by = 7),
    window_len = protocol_durations,
    species_f = factor(species_levels, levels = species_levels)
  ) |>
    as_tibble() |>
    bind_cols(as_tibble(median_covs))
  tmp_grid$pred <- predict(model, newdata = tmp_grid, type = "response",
                           exclude = "s(ds_sp_f)")
  opt <- tmp_grid |>
    summarise(pred = mean(pred), .by = c(day_start, window_len)) |>
    summarise(
      best_day = day_start[which.min(pred)],
      best_pred = round(min(pred), 5),
      .by = window_len
    ) |>
    arrange(window_len) |>
    mutate(metric = metric_name)
  for (i in seq_len(nrow(opt))) {
    cat(sprintf("  %s: %dd optimal start=day %d (pred=%.5f)\n",
                metric_name, opt$window_len[i], opt$best_day[i], opt$best_pred[i]))
  }
}

write.csv(q6_summary, "Q6_protocol_evaluation.csv", row.names = FALSE)
cat("  Saved Q6_protocol_evaluation.csv\n\n")


# ============================================================================
# Q7: OPTIMAL WINDOW TIMING
# ============================================================================
# For each species/guild, what timing minimises deviation at key durations?

cat("── Q7: Optimal window timing ──────────────────────────────────\n")

q7 <- preds_sp |>
  filter(window_len %in% c(15, 29, 64, 92, 120, 155, 183)) |>
  summarise(
    best_day_start  = day_start[which.min(pred)],
    min_pred        = min(pred),
    worst_day_start = day_start[which.max(pred)],
    max_pred        = max(pred),
    timing_ratio    = max(pred) / min(pred),
    .by = c(species_f, guild_major, guild_minor_diet, window_len)
  ) |>
  arrange(window_len, desc(timing_ratio))

write.csv(q7, "Q7_optimal_timing.csv", row.names = FALSE)
cat("  Saved Q7_optimal_timing.csv\n")

# Guild-level optimum at ~60 days
q7 |>
  filter(window_len == 64) |>
  summarise(
    best_day = round(median(best_day_start)),
    worst_day = round(median(worst_day_start)),
    median_ratio = round(median(timing_ratio), 1),
    .by = guild_major
  ) |>
  print()
cat("\n")


# ============================================================================
# Q8: RICHNESS & SIMILARITY SURFACE
# ============================================================================
# How do community-level metrics vary across the surface?

cat("── Q8: Richness & similarity surface ──────────────────────────\n")

# Coarser grid for richness
rich_grid_coarse <- expand.grid(
  day_start  = seq(1, 358, by = 7),
  window_len = seq(15, 183, by = 7)
) |>
  as_tibble() |>
  mutate(
    s_bio4 = 0, l_trapdays = median(sens_richness$l_trapdays),
    l_nsites = median(sens_richness$l_nsites), s_latitude = 0,
    dataset_f = sens_richness$dataset_f[1]
  )

# Prep beta-model data frames with adjusted columns
rich_grid_beta <- rich_grid_coarse
rich_grid_rho  <- rich_grid_coarse

q8 <- list()

# d_sr_raref (Gaussian)
rich_grid_coarse$pred <- predict(mod_richness, newdata = rich_grid_coarse,
                                 type = "response", exclude = "s(dataset_f)")
q8$d_sr_raref <- rich_grid_coarse |>
  select(day_start, window_len, pred) |>
  mutate(metric = "d_sr_raref")

# prop_sr_full (Beta — need the adjusted response name in data but prediction is fine)
rich_grid_beta$pred <- predict(mod_prop, newdata = rich_grid_beta,
                               type = "response", exclude = "s(dataset_f)")
q8$prop_sr_full <- rich_grid_beta |>
  select(day_start, window_len, pred) |>
  mutate(metric = "prop_sr_full")

# rho_lambda (Beta on transformed scale — back-transform: rho = pred * 2 - 1)
# Model was fitted on windows with >= 5 shared species; use that subset's medians
sens_rich_rho_filt <- sens_richness |>
  filter(!is.na(rho_lambda), n_shared_spp >= 5)
rich_grid_rho <- rich_grid_rho |>
  mutate(l_trapdays = median(sens_rich_rho_filt$l_trapdays),
         l_nsites   = median(sens_rich_rho_filt$l_nsites),
         dataset_f  = sens_rich_rho_filt$dataset_f[1])
rich_grid_rho$pred <- predict(mod_rho, newdata = rich_grid_rho,
                              type = "response", exclude = "s(dataset_f)")
q8$rho_lambda <- rich_grid_rho |>
  select(day_start, window_len, pred) |>
  mutate(pred = pred * 2 - 1,  # back-transform to [-1, 1]
         metric = "rho_lambda")

q8_out <- bind_rows(q8)
write.csv(q8_out, "Q8_richness_surface.csv", row.names = FALSE)
cat("  Saved Q8_richness_surface.csv\n")

# Best timing for each community metric at 60 days
q8_out |>
  filter(window_len == 64) |>
  summarise(
    best_start = day_start[which.max(abs(pred) * ifelse(metric == "d_sr_raref", -1, 1))],
    best_val   = round(pred[which.max(abs(pred) * ifelse(metric == "d_sr_raref", -1, 1))], 3),
    .by = metric
  ) |>
  print()


# ============================================================================
# Q9: BENCHMARK NOISE FLOOR & SIGNAL-TO-NOISE RATIO
# ============================================================================
# How large are window deviations relative to inter-annual benchmark variability?
# Uses multi-year sites (BE-Leuven, NO-evenstadlia, SI-serknica, SP-donana)
# to define a species-specific noise floor.
#
# SNR = |d_lambda| / SD_benchmark
# When SNR < 1, the deviation is smaller than year-to-year benchmark fluctuation.

cat("── Q9: Benchmark noise floor ──────────────────────────────────\n")

# ── 9a. Load inter-annual CV data ──
cv_data <- read.csv("interannual_cv_lambda_full.csv")
cat("  Inter-annual CV data:", nrow(cv_data), "species × site combinations\n")

# ── 9b. Recover lambda_full and compute SNR per observation ──
# Restrict to multi-year datasets where we have inter-annual SD
multi_year_bases <- unique(cv_data$base_dataset)

noise_floor <- all_window_species |>
  filter(
    window_id != "FULL",
    !window_id %in% c("SNAP_EU_CORE", "SNAP_EU_BUFFER", "EOW_EARLY", "EOW_LATE")
  ) |>
  mutate(base_dataset = gsub("_slice\\d+$", "", dataset)) |>
  filter(base_dataset %in% multi_year_bases) |>
  left_join(
    cv_data |> select(base_dataset, species, sd_lambda, cv),
    by = c("base_dataset", "species")
  ) |>
  filter(!is.na(sd_lambda)) |>
  mutate(
    abs_d_lambda = abs(d_lambda),
    snr = abs_d_lambda / sd_lambda
  )

cat("  Observations with SNR data:", nrow(noise_floor),
    "(", n_distinct(noise_floor$species), "species )\n")

# ── 9c. Summarise SNR by window duration ──
snr_by_duration <- noise_floor |>
  summarise(
    median_snr    = median(snr),
    q25_snr       = quantile(snr, 0.25),
    q75_snr       = quantile(snr, 0.75),
    pct_below_1   = mean(snr < 1) * 100,
    pct_below_2   = mean(snr < 2) * 100,
    median_abs_d  = median(abs_d_lambda),
    median_sd     = median(sd_lambda),
    n_obs         = n(),
    .by = window_len
  ) |>
  mutate(level = "Overall") |>
  arrange(window_len)

# ── 9d. Per-species SNR summary at long windows (≥85d) ──
snr_by_species <- noise_floor |>
  filter(window_len >= 85) |>
  summarise(
    median_snr  = median(snr),
    pct_below_1 = mean(snr < 1) * 100,
    cv_pct      = first(cv),
    n_obs       = n(),
    .by = c(base_dataset, species)
  ) |>
  arrange(desc(pct_below_1))

# ── 9e. Combine and save ──
q9_out <- bind_rows(
  snr_by_duration |>
    mutate(base_dataset = NA_character_, species = NA_character_,
           section = "duration_summary"),
  snr_by_species |>
    mutate(window_len = NA_integer_, q25_snr = NA_real_, q75_snr = NA_real_,
           pct_below_2 = NA_real_, median_abs_d = NA_real_, median_sd = NA_real_,
           level = NA_character_, section = "species_summary")
)

write.csv(q9_out, "Q9_noise_floor.csv", row.names = FALSE)
cat("  Saved Q9_noise_floor.csv\n")

# Print key results
cat("\n  --- SNR by window duration (overall) ---\n")
snr_by_duration |>
  filter(window_len %in% c(15, 29, 57, 85, 120, 155, 183)) |>
  select(window_len, median_snr, pct_below_1, n_obs) |>
  mutate(median_snr = round(median_snr, 1),
         pct_below_1 = round(pct_below_1, 2)) |>
  print()

cat("\n  --- Species approaching noise floor (≥85d windows, pct SNR<1 > 0) ---\n")
snr_by_species |>
  filter(pct_below_1 > 0) |>
  mutate(median_snr = round(median_snr, 1),
         pct_below_1 = round(pct_below_1, 1),
         cv_pct = round(cv_pct, 1)) |>
  print(n = 20)
cat("\n")


cat("\n=== Done. All Q1-Q9 results saved. ===\n")

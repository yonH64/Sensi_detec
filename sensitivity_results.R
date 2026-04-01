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
#         Q5_latitude_effect.csv
#         Q6_protocol_evaluation.csv + Q6_protocol_species.csv
#         Q7_optimal_timing.csv
#         Q8_richness_surface.csv
#         Q9_noise_floor.csv
#         Q10_metric_agreement.csv + Q10_metric_agreement_species.csv
#         model_comparison_table.csv
# ============================================================================

library(mgcv)
library(dplyr)
library(tidyr)
library(tibble)

cat("=== Sensitivity Results Extraction ===\n\n")

# Everything runs inside local() so intermediate objects don't persist.
# Returns a compact summary list assigned to `sensitivity_summary`.
sensitivity_summary <- local({

# ── Load ─────────────────────────────────────────────────────────────────────

load("sensitivity_models_env.RData", envir = environment())

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

# Median covariate values for "typical" predictions.
# ds_sp_f: mgcv requires a valid factor level even for excluded RE terms;
# the arbitrary level is nullified by exclude = "s(ds_sp_f)" in predict().
median_covs <- list(
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


# ── Posterior simulation helpers ───────────────────────────────────────────────
# Three-layer design:
#   compute_pred_draws()  — draws from coefficient posterior → response-scale matrix
#   aggregate_draws()     — groups pred_draws and returns CI summary
#   posterior_ci_single() — per-row CIs (no grouping; for richness / species-level)
#
# For Q1-Q3 (lambda), draws are computed ONCE and reused across aggregations.
# For other metrics and Q6/Q8, compute per-model.

N_DRAWS <- 1000L
POSTERIOR_SEED <- 4721L

compute_pred_draws <- function(model, grid, n_draws = N_DRAWS,
                                seed = POSTERIOR_SEED,
                                exclude = "s(ds_sp_f)") {
  set.seed(seed)
  Xp <- predict(model, newdata = grid, type = "lpmatrix", exclude = exclude)
  beta_draws <- MASS::mvrnorm(n_draws, coef(model), vcov(model))
  eta_draws <- Xp %*% t(beta_draws)
  pred_draws <- model$family$linkinv(eta_draws)
  rm(Xp, beta_draws, eta_draws); gc(verbose = FALSE)
  pred_draws   # n_grid × n_draws matrix
}

aggregate_draws <- function(pred_draws, grid, group_vars) {
  grid_aug <- grid |> mutate(.row = row_number())
  groups <- grid_aug |> distinct(across(all_of(group_vars))) |>
    mutate(.gid = row_number())
  grid_aug <- grid_aug |> left_join(groups, by = group_vars)

  n_groups <- nrow(groups)
  n_draws <- ncol(pred_draws)
  avg_draws <- matrix(NA_real_, n_groups, n_draws)
  for (g in seq_len(n_groups)) {
    rows <- grid_aug$.row[grid_aug$.gid == g]
    avg_draws[g, ] <- colMeans(pred_draws[rows, , drop = FALSE])
  }

  groups |>
    mutate(
      ci_lo     = apply(avg_draws, 1, quantile, probs = 0.025),
      ci_hi     = apply(avg_draws, 1, quantile, probs = 0.975),
      post_sd   = apply(avg_draws, 1, sd),
      post_mean = rowMeans(avg_draws),
      rse       = post_sd / abs(post_mean)
    ) |>
    select(-.gid)
}

posterior_ci_single <- function(model, grid, n_draws = N_DRAWS,
                                 seed = POSTERIOR_SEED,
                                 exclude = "s(dataset_f)") {
  pred_draws <- compute_pred_draws(model, grid, n_draws, seed, exclude)
  tibble(
    ci_lo = apply(pred_draws, 1, quantile, probs = 0.025),
    ci_hi = apply(pred_draws, 1, quantile, probs = 0.975),
    se    = apply(pred_draws, 1, sd)
  )
}


# ── Coarse species grid (shared by Q1-Q4, Q10) ───────────────────────────────
# 52 day_start × 25 window_len × 29 species = ~37,700 rows.
# Fine enough for surfaces; small enough for posterior simulation (~720 MB peak).
surf_grid <- expand.grid(
  day_start  = seq(1, 358, by = 7),
  window_len = seq(15, 183, by = 7),
  species_f  = factor(species_levels, levels = species_levels),
  stringsAsFactors = FALSE
) |>
  as_tibble() |>
  bind_cols(as_tibble(median_covs))

# Version with guild info for aggregation
surf_grid_g <- surf_grid |> left_join(guild_info, by = "species_f")

# Standard key durations (all on the coarse grid seq(15, 183, by = 7))
key_durations <- c(15, 29, 57, 85, 120, 155, 183)


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

# Point predictions on coarse grid (reused in Q2-Q4)
preds_lambda <- predict_detection(mod_lambda, grid = surf_grid)

# Posterior draws for lambda on surf_grid (reused in Q1-Q3)
cat("  Computing posterior draws for lambda...\n")
draws_lambda <- compute_pred_draws(mod_lambda, surf_grid)
cat("  Done (", ncol(draws_lambda), "draws ×", nrow(draws_lambda), "grid points)\n")

# Per-guild summaries
q1 <- preds_lambda |>
  summarise(
    mean_pred = mean(pred),
    median_pred = median(pred),
    q25 = quantile(pred, 0.25),
    q75 = quantile(pred, 0.75),
    .by = c(window_len, guild_major)
  ) |>
  arrange(guild_major, window_len)

q1_ci <- aggregate_draws(draws_lambda, surf_grid_g,
                          c("window_len", "guild_major")) |>
  select(window_len, guild_major, ci_lo, ci_hi, rse)

q1 <- q1 |> left_join(q1_ci, by = c("window_len", "guild_major"))

# Overall (across all species + day_start positions)
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

q1_ci_all <- aggregate_draws(draws_lambda, surf_grid_g, "window_len") |>
  select(window_len, ci_lo, ci_hi, rse)

q1_overall <- q1_overall |> left_join(q1_ci_all, by = "window_len")

q1_out <- bind_rows(q1, q1_overall)
write.csv(q1_out, "Q1_duration_effect.csv", row.names = FALSE)
cat("  Saved Q1_duration_effect.csv\n")

# Key insight: how much deviation drops from 15 → 60 → 120 days
q1_overall |>
  filter(window_len %in% c(15, 29, 57, 92, 120, 148, 183)) |>
  mutate(across(c(mean_pred, median_pred, ci_lo, ci_hi), ~round(., 4)),
         rse = round(rse, 3)) |>
  print()
cat("\n")


# ============================================================================
# Q2: SEASONAL PROFILES
# ============================================================================
# How does deviation vary across the year, at fixed durations?

cat("── Q2: Seasonal profiles ──────────────────────────────────────\n")

# Point estimates + IQR (species quantile bands)
q2 <- preds_lambda |>
  filter(window_len %in% key_durations) |>
  summarise(
    mean_pred = mean(pred),
    q25 = quantile(pred, 0.25),
    q75 = quantile(pred, 0.75),
    .by = c(day_start, window_len, guild_major)
  )

q2_overall <- preds_lambda |>
  filter(window_len %in% key_durations) |>
  summarise(
    mean_pred = mean(pred),
    q25 = quantile(pred, 0.25),
    q75 = quantile(pred, 0.75),
    .by = c(day_start, window_len)
  ) |>
  mutate(guild_major = "All")

# Posterior CIs on species-averaged predictions (reuse draws_lambda)
q2_ci <- aggregate_draws(draws_lambda, surf_grid_g,
                          c("day_start", "window_len", "guild_major")) |>
  filter(window_len %in% key_durations) |>
  select(day_start, window_len, guild_major, ci_lo, ci_hi)

q2_ci_all <- aggregate_draws(draws_lambda, surf_grid_g,
                              c("day_start", "window_len")) |>
  filter(window_len %in% key_durations) |>
  select(day_start, window_len, ci_lo, ci_hi)

q2 <- q2 |> left_join(q2_ci, by = c("day_start", "window_len", "guild_major"))
q2_overall <- q2_overall |> left_join(q2_ci_all, by = c("day_start", "window_len"))

q2_out <- bind_rows(q2, q2_overall)
write.csv(q2_out, "Q2_seasonal_profiles.csv", row.names = FALSE)
cat("  Saved Q2_seasonal_profiles.csv\n")

# Worst season vs best season at ~60 days
q2_overall |>
  filter(window_len == 57) |>
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
# Uses surf_grid (defined in scaffold section above).

cat("── Q3: Full surface predictions ───────────────────────────────\n")

# Helper: compute point estimates (mean, IQR) + posterior CI for one metric
q3_one_metric <- function(model, metric_name, draws = NULL) {
  surf_grid$pred <- predict(model, newdata = surf_grid, type = "response",
                            exclude = "s(ds_sp_f)")
  pts <- surf_grid |>
    left_join(guild_info, by = "species_f") |>
    summarise(
      mean_pred = mean(pred),
      q25 = quantile(pred, 0.25),
      q75 = quantile(pred, 0.75),
      .by = c(day_start, window_len, guild_major)
    ) |>
    mutate(metric = metric_name)

  # Posterior CIs
  if (is.null(draws)) {
    cat("    Computing posterior draws for", metric_name, "...\n")
    draws <- compute_pred_draws(model, surf_grid)
  }
  ci <- aggregate_draws(draws, surf_grid_g,
                         c("day_start", "window_len", "guild_major")) |>
    select(day_start, window_len, guild_major, ci_lo, ci_hi)

  pts |> left_join(ci, by = c("day_start", "window_len", "guild_major"))
}

q3 <- list()

# Lambda — reuse existing draws
q3$abs_d_lambda <- q3_one_metric(mod_lambda, "abs_d_lambda", draws = draws_lambda)
rm(draws_lambda); gc(verbose = FALSE)

# Other absolute deviation metrics (compute draws fresh, discard after)
q3$abs_d_rate <- q3_one_metric(mod_rate, "abs_d_rate")
q3$abs_d_matched_rate <- q3_one_metric(mod_matched, "abs_d_matched_rate")

# Signed rate (Gaussian — predictions can be negative)
q3$d_rate_signed <- q3_one_metric(mod_rate_signed, "d_rate_signed")

q3_out <- bind_rows(q3)
write.csv(q3_out, "Q3_surface_predictions.csv", row.names = FALSE)
cat("  Saved Q3_surface_predictions.csv (", nrow(q3_out), "rows)\n\n")


# ============================================================================
# Q4: SPECIES/GUILD SURFACE VARIATION
# ============================================================================
# How different are species-specific surfaces?

cat("── Q4: Species/guild surface variation ────────────────────────\n")

preds_sp <- preds_lambda  # reuse Q1 predictions (same model + grid)

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
# Q5: LATITUDE EFFECT
# ============================================================================
# How does latitude modulate deviation?
# (BIO4 was dropped due to collinearity with latitude (r = 0.77) and LOO
#  instability; latitude absorbs the shared environmental signal.)

cat("── Q5: Latitude effect ────────────────────────────────────────\n")

# Extract parametric coefficient for s_latitude from each detection model
q5 <- tibble(
  metric = c("abs_d_lambda", "abs_d_rate", "abs_d_matched_rate", "d_rate_signed"),
  model  = list(mod_lambda, mod_rate, mod_matched, mod_rate_signed)
) |>
  rowwise() |>
  mutate(
    coef   = coef(model)["s_latitude"],
    se     = summary(model)$p.table["s_latitude", "Std. Error"],
    t_val  = summary(model)$p.table["s_latitude", "t value"],
    p_val  = summary(model)$p.table["s_latitude", "Pr(>|t|)"]
  ) |>
  ungroup() |>
  select(-model) |>
  mutate(across(c(coef, se, t_val), ~round(., 4)),
         p_val = signif(p_val, 3))

print(q5)
write.csv(q5, "Q5_latitude_effect.csv", row.names = FALSE)
cat("  Saved Q5_latitude_effect.csv\n\n")


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

protocol_grid_g <- protocol_grid |>
  left_join(guild_info, by = "species_f")

# ── Predict + posterior CIs per metric ──
q6_list <- list()
q6_ci_guild_list <- list()
q6_ci_overall_list <- list()

for (metric_name in c("abs_d_lambda", "abs_d_rate", "abs_d_matched_rate")) {
  model <- all_models$detection[[metric_name]]

  # Point predictions + species-level SE via delta method
  p <- predict(model, newdata = protocol_grid, type = "link",
               se.fit = TRUE, exclude = "s(ds_sp_f)")
  linkinv <- model$family$linkinv
  pred_resp <- linkinv(p$fit)
  # Delta method SE on response scale: d/deta[linkinv(eta)] * se_eta
  se_resp <- abs(pred_resp) * p$se.fit   # for log link, d/deta[exp(eta)] = exp(eta)

  q6_list[[metric_name]] <- protocol_grid_g |>
    mutate(
      pred  = pred_resp,
      se    = se_resp,
      ci_lo = linkinv(p$fit - 1.96 * p$se.fit),
      ci_hi = linkinv(p$fit + 1.96 * p$se.fit),
      metric = metric_name
    ) |>
    select(protocol, day_start, window_len, species_f, guild_major,
           pred, se, ci_lo, ci_hi, metric)

  # Posterior CIs on species-averaged predictions (small grid — fast)
  cat("    Posterior draws for", metric_name, "...\n")
  draws_prot <- compute_pred_draws(model, protocol_grid)

  q6_ci_guild_list[[metric_name]] <-
    aggregate_draws(draws_prot, protocol_grid_g,
                     c("protocol", "guild_major")) |>
    mutate(metric = metric_name)

  q6_ci_overall_list[[metric_name]] <-
    aggregate_draws(draws_prot, protocol_grid_g, "protocol") |>
    mutate(metric = metric_name)

  rm(draws_prot); gc(verbose = FALSE)
}

q6_all <- bind_rows(q6_list)

# ── Guild-level summary (long format with CIs) ──
q6_guild_pts <- q6_all |>
  summarise(
    mean_pred = mean(pred),
    q25 = quantile(pred, 0.25),
    q75 = quantile(pred, 0.75),
    .by = c(protocol, guild_major, metric)
  )

q6_ci_guild <- bind_rows(q6_ci_guild_list) |>
  select(protocol, guild_major, metric, ci_lo, ci_hi, rse)

q6_summary <- q6_guild_pts |>
  left_join(q6_ci_guild, by = c("protocol", "guild_major", "metric"))

# ── Overall summary (species-averaged, long format with CIs) ──
q6_overall_pts <- q6_all |>
  summarise(
    mean_pred = mean(pred),
    q25 = quantile(pred, 0.25),
    q75 = quantile(pred, 0.75),
    .by = c(protocol, metric)
  ) |>
  mutate(guild_major = "All")

q6_ci_overall <- bind_rows(q6_ci_overall_list) |>
  select(protocol, metric, ci_lo, ci_hi, rse)

q6_overall <- q6_overall_pts |>
  left_join(q6_ci_overall, by = c("protocol", "metric"))

# Console printout in wide format (for readability)
q6_summary |>
  select(guild_major, metric, protocol, mean_pred) |>
  pivot_wider(names_from = protocol, values_from = mean_pred) |>
  mutate(across(c(CORE, BUFFER, EOW_EARLY, EOW_LATE), ~round(., 5))) |>
  print(n = 25)

cat("\n  Overall (species-averaged):\n")
q6_overall |>
  select(metric, protocol, mean_pred) |>
  pivot_wider(names_from = protocol, values_from = mean_pred) |>
  mutate(across(c(CORE, BUFFER, EOW_EARLY, EOW_LATE), ~round(., 5))) |>
  print()

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

# Save long-format summary + species-level with SEs
q6_summary_out <- bind_rows(q6_summary, q6_overall)
write.csv(q6_summary_out, "Q6_protocol_evaluation.csv", row.names = FALSE)
write.csv(q6_all, "Q6_protocol_species.csv", row.names = FALSE)
cat("  Saved Q6_protocol_evaluation.csv (long format) + Q6_protocol_species.csv\n\n")


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
    l_trapdays = median(sens_richness$l_trapdays),
    l_nsites = median(sens_richness$l_nsites), s_latitude = 0,
    dataset_f = sens_richness$dataset_f[1]
  )

q8 <- list()

# ── d_sr_raref (Gaussian) ──
cat("    Posterior CIs for d_sr_raref...\n")
rich_grid_coarse$pred <- predict(mod_richness, newdata = rich_grid_coarse,
                                 type = "response", exclude = "s(dataset_f)")
ci_raref <- posterior_ci_single(mod_richness, rich_grid_coarse)
q8$d_sr_raref <- rich_grid_coarse |>
  select(day_start, window_len, pred) |>
  bind_cols(ci_raref) |>
  mutate(metric = "d_sr_raref")

# ── prop_sr_full (Beta) ──
cat("    Posterior CIs for prop_sr_full...\n")
rich_grid_beta <- rich_grid_coarse |> select(-pred)
rich_grid_beta$pred <- predict(mod_prop, newdata = rich_grid_beta,
                               type = "response", exclude = "s(dataset_f)")
ci_prop <- posterior_ci_single(mod_prop, rich_grid_beta)
q8$prop_sr_full <- rich_grid_beta |>
  select(day_start, window_len, pred) |>
  bind_cols(ci_prop) |>
  mutate(metric = "prop_sr_full")

# ── rho_lambda (Beta on transformed scale) ──
# Model was fitted on windows with >= 5 shared species; use that subset's medians
cat("    Posterior CIs for rho_lambda...\n")
sens_rich_rho_filt <- sens_richness |>
  filter(!is.na(rho_lambda), n_shared_spp >= 5)
rich_grid_rho <- rich_grid_coarse |>
  select(day_start, window_len) |>
  mutate(
    l_trapdays = median(sens_rich_rho_filt$l_trapdays),
    l_nsites   = median(sens_rich_rho_filt$l_nsites),
    s_latitude = 0,
    dataset_f  = sens_rich_rho_filt$dataset_f[1]
  )
rich_grid_rho$pred <- predict(mod_rho, newdata = rich_grid_rho,
                              type = "response", exclude = "s(dataset_f)")

# Posterior draws for rho — back-transform from (0,1) beta to (-1,1) correlation
ci_rho_raw <- posterior_ci_single(mod_rho, rich_grid_rho)
# Back-transform: rho = pred * 2 - 1
q8$rho_lambda <- rich_grid_rho |>
  select(day_start, window_len, pred) |>
  mutate(
    pred  = pred * 2 - 1,
    ci_lo = ci_rho_raw$ci_lo * 2 - 1,
    ci_hi = ci_rho_raw$ci_hi * 2 - 1,
    se    = ci_rho_raw$se * 2,  # SE scales linearly with the affine transform
    metric = "rho_lambda"
  )

q8_out <- bind_rows(q8)
write.csv(q8_out, "Q8_richness_surface.csv", row.names = FALSE)
cat("  Saved Q8_richness_surface.csv\n")

# Best timing for each community metric at 64 days
# d_sr_raref: best = closest to zero; prop_sr_full / rho_lambda: best = highest
q8_out |>
  filter(window_len == 64) |>
  summarise(
    best_start = case_when(
      first(metric) == "d_sr_raref" ~ day_start[which.min(abs(pred))],
      TRUE ~ day_start[which.max(pred)]
    ),
    best_val = case_when(
      first(metric) == "d_sr_raref" ~ round(pred[which.min(abs(pred))], 3),
      TRUE ~ round(pred[which.max(pred)], 3)
    ),
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

cat("  Saved Q9_noise_floor.csv\n\n")


# ============================================================================
# Q10: METRIC AGREEMENT
# ============================================================================
# Do the three absolute deviation metrics (lambda, rate, matched_rate) produce
# consistent sensitivity surfaces?

cat("── Q10: Metric agreement ──────────────────────────────────────\n")

# Predict all three metrics on the same coarse species grid
q10_grid <- surf_grid  # reuse from Q3

q10_preds <- list()
for (metric_name in c("abs_d_lambda", "abs_d_rate", "abs_d_matched_rate")) {
  model <- all_models$detection[[metric_name]]
  q10_grid$pred <- predict(model, newdata = q10_grid, type = "response",
                           exclude = "s(ds_sp_f)")
  q10_preds[[metric_name]] <- q10_grid |>
    transmute(day_start, window_len, species_f, pred) |>
    rename(!!metric_name := pred)
}

q10_merged <- q10_preds$abs_d_lambda |>
  left_join(q10_preds$abs_d_rate, by = c("day_start", "window_len", "species_f")) |>
  left_join(q10_preds$abs_d_matched_rate, by = c("day_start", "window_len", "species_f"))

# Pairwise correlations (species-level predicted surfaces)
pairs <- list(
  c("abs_d_lambda", "abs_d_rate"),
  c("abs_d_lambda", "abs_d_matched_rate"),
  c("abs_d_rate", "abs_d_matched_rate")
)

q10_overall <- tibble(
  metric_a = sapply(pairs, `[`, 1),
  metric_b = sapply(pairs, `[`, 2),
  pearson_r = sapply(pairs, \(p)
    round(cor(q10_merged[[p[1]]], q10_merged[[p[2]]]), 3)),
  spearman_rho = sapply(pairs, \(p)
    round(cor(q10_merged[[p[1]]], q10_merged[[p[2]]], method = "spearman"), 3))
)

# Per-species surface correlations
q10_by_species <- q10_merged |>
  summarise(
    r_lambda_rate    = cor(abs_d_lambda, abs_d_rate),
    r_lambda_matched = cor(abs_d_lambda, abs_d_matched_rate),
    r_rate_matched   = cor(abs_d_rate, abs_d_matched_rate),
    .by = species_f
  ) |>
  left_join(guild_info, by = "species_f")

q10_out <- list(
  overall = q10_overall,
  by_species = q10_by_species |>
    mutate(across(starts_with("r_"), ~round(., 3)))
)

print(q10_overall)
cat("\n  Per-species surface correlation summary:\n")
q10_by_species |>
  summarise(across(starts_with("r_"), \(x)
    paste0("median=", round(median(x), 3), " [", round(min(x), 3), ", ",
           round(max(x), 3), "]"))) |>
  glimpse()

write.csv(q10_overall, "Q10_metric_agreement.csv", row.names = FALSE)
write.csv(q10_out$by_species, "Q10_metric_agreement_species.csv", row.names = FALSE)
cat("  Saved Q10_metric_agreement.csv + Q10_metric_agreement_species.csv\n\n")


# ============================================================================
# CONSOLIDATED SUMMARY
# ============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║              KEY RESULTS SUMMARY                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# ── Model performance ──
cat("── Model performance ──────────────────────────────────────────\n")
det_models <- list(
  abs_d_lambda       = mod_lambda,
  abs_d_rate         = mod_rate,
  abs_d_matched_rate = mod_matched,
  d_rate_signed      = mod_rate_signed
)
com_models <- list(
  d_sr_raref   = mod_richness,
  prop_sr_full = mod_prop,
  rho_lambda   = mod_rho
)

model_perf <- tibble(
  model = c(names(det_models), names(com_models)),
  family = c(sapply(det_models, \(m) m$family$family),
             sapply(com_models, \(m) m$family$family)),
  N = c(sapply(det_models, \(m) m$df.null + 1),
        sapply(com_models, \(m) m$df.null + 1)),
  dev_expl = c(sapply(det_models, \(m) summary(m)$dev.expl),
               sapply(com_models, \(m) summary(m)$dev.expl))
) |>
  mutate(dev_expl_pct = paste0(round(dev_expl * 100, 1), "%"))

print(model_perf |> select(model, family, N, dev_expl_pct))

# ── All parametric coefficients (detection models) ──
cat("\n── Parametric coefficients (detection models) ─────────────────\n")
param_names <- c("s_latitude", "s_trap_array", "l_trapdays", "l_nsites")

param_tbl <- purrr::map_dfr(names(det_models), \(metric_name) {
  m <- det_models[[metric_name]]
  pt <- summary(m)$p.table
  available <- intersect(param_names, rownames(pt))
  tibble(
    metric    = metric_name,
    term      = available,
    estimate  = round(pt[available, "Estimate"], 4),
    se        = round(pt[available, "Std. Error"], 4),
    t_value   = round(pt[available, "t value"], 2),
    p_value   = signif(pt[available, "Pr(>|t|)"], 3)
  )
})

print(param_tbl, n = 20)

# ── Duration effect (headline numbers) ──
cat("\n── Duration effect (mean |d_lambda|, all species, with 95% CI) ─\n")
q1_overall |>
  filter(window_len %in% key_durations) |>
  mutate(across(c(mean_pred, median_pred, ci_lo, ci_hi), ~round(., 4)),
         rse = round(rse, 3)) |>
  select(window_len, mean_pred, ci_lo, ci_hi, rse) |>
  print()

# ── Seasonal timing (~57-day windows) ──
cat("\n── Seasonal timing at ~57d (lambda, all species) ──────────────\n")
q2_overall |>
  filter(window_len == 57) |>
  summarise(
    best_start  = day_start[which.min(mean_pred)],
    best_pred   = round(min(mean_pred), 5),
    worst_start = day_start[which.max(mean_pred)],
    worst_pred  = round(max(mean_pred), 5),
    ratio       = round(max(mean_pred) / min(mean_pred), 1)
  ) |> print()

# ── Protocol evaluation (overall) ──
cat("\n── Protocol evaluation (species-averaged |d_lambda|, with 95% CI) ─\n")
q6_overall |>
  filter(metric == "abs_d_lambda") |>
  mutate(across(c(mean_pred, ci_lo, ci_hi), ~round(., 5)),
         rse = round(rse, 3)) |>
  select(protocol, mean_pred, ci_lo, ci_hi, rse) |>
  print()

# ── Richness recovery at 64d ──
cat("\n── Richness recovery at ~64 days ───────────────────────────────\n")
q8_out |>
  filter(window_len == 64) |>
  summarise(
    best_start = case_when(
      first(metric) == "d_sr_raref" ~ day_start[which.min(abs(pred))],
      TRUE ~ day_start[which.max(pred)]
    ),
    best_val = case_when(
      first(metric) == "d_sr_raref" ~ round(pred[which.min(abs(pred))], 3),
      TRUE ~ round(pred[which.max(pred)], 3)
    ),
    .by = metric
  ) |> print()

# ── Metric agreement ──
cat("\n── Metric agreement (surface correlations) ─────────────────────\n")
print(q10_overall)

# ── Noise floor ──
cat("\n── Noise floor (SNR by window duration) ────────────────────────\n")
snr_by_duration |>
  filter(window_len %in% c(15, 57, 85, 120, 183)) |>
  select(window_len, median_snr, pct_below_1, n_obs) |>
  mutate(median_snr = round(median_snr, 1),
         pct_below_1 = round(pct_below_1, 2)) |>
  print()

cat("\n  Species approaching noise floor (≥85d, pct SNR<1 > 0):\n")
snr_by_species |>
  filter(pct_below_1 > 0) |>
  mutate(median_snr = round(median_snr, 1),
         pct_below_1 = round(pct_below_1, 1),
         cv_pct = round(cv_pct, 1)) |>
  print(n = 20)

cat("\n=== Done. All Q1-Q10 results saved. ===\n")

# Return compact summary (the only object that persists)
list(
  model_performance = model_perf |> select(model, family, N, dev_expl_pct),
  parametric_coefficients = param_tbl,
  latitude_effect = q5,
  metric_agreement = q10_overall
)

}) # end local()

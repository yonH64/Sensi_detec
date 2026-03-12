############################################################
#   PROTOCOL EVALUATION: MSE-BASED ANALYSIS
#   SPECIES-LEVEL DETECTION METRICS (TTE / LRR / SCOV)
#
#   OBJECTIVE: Compare protocols using bias-variance tradeoff
#   via Mean Squared Error (MSE) rather than absolute deviation
#
#   MSE = Bias² + Variance
#   where:
#     - Bias = (estimate - truth)
#     - Variance = SE²
#
#   METRICS ANALYZED:
#   1. MSE of TTE-derived detection probability  (mse_p_tte)
#   2. MSE of event rate (log scale)             (mse_log_rate)
#   3. MSE of spatial detection coverage         (mse_spatial_cov)
#
#   NOTE: richness MSE models (mse_sr_raref) are in
#         models_mse_richness.R
#
#   UPDATED: 2026-02-20
#   - Added l_n_species covariate (portfolio effect)
#   - Removed family-based models (insufficient data: 9 families)
#   - Added posterior predictive checks for validation
#   - Added richness × protocol interaction models
#
############################################################

# ──────────────────────────────────────────────────────────
# OVERVIEW OF MODELS
# ──────────────────────────────────────────────────────────
#
# RESPONSE VARIABLES (MSE from FULL):
#   - mse_p_tte        : (p_tte - p_tte_full)² + p_tte_se²         [0,∞)
#   - mse_log_rate     : (log_rate - log_rate_full)² + log_rate_se² [0,∞)
#   - mse_spatial_cov  : (spatial_cov - spatial_cov_full)² + spatial_cov_se² [0,∞)
#
# FIXED EFFECTS (in all models):
#   - protocol_bin            : BUFFER=1, CORE=0 (focal predictor)
#   - log1p(trap_days_window) : Sampling effort (trap-days)
#   - log1p(n_sites)          : Number of camera sites
#   - scale(latitude)         : Geographic latitude
#   - scale(log1p(trap_array)): Spatial extent of camera array
#   - log1p(n_species)        : Species richness (portfolio effect)
#
# RANDOM EFFECTS:
#   - (1 | ds_sp)                    : Dataset × Species interaction
#   - (1 | species)                  : Species-level variation (base models)
#   - (1 + protocol_bin || species)  : Species-specific protocol effects
#   - guild_major                    : Functional guild (fixed interaction)
#
# STATISTICAL FAMILY:
#   - Gamma(link="log")  : MSE responses (positive continuous, right-skewed)
#
# INTERPRETATION:
#   - Smaller MSE = better overall approximation (less bias + less variance)
#   - Negative protocol_bin coefficient = BUFFER has lower MSE (better)
#   - Positive covariate coefficient = larger values → larger MSE (worse)
#   - Negative l_n_species coefficient = portfolio effect reduces MSE
#
############################################################

# ---- Packages ----
library(dplyr)
library(brms)
library(readr)
library(tidyr)
library(tibble)
library(posterior)
library(purrr)


############################################################
# 1) DATA PREPARATION
############################################################

prep_mse_summaries <- function(all_window_species,
                                protocol_keep = c("SNAP_EU_CORE","SNAP_EU_BUFFER")) {
  
  req <- c("dataset","species","window_id",
           "mse_p_tte","mse_log_rate","mse_spatial_cov",
           "trap_days_window","n_sites","latitude","trap_array",
           "trap_days_full","n_sites_full")
  stopifnot(all(req %in% names(all_window_species)))
  
  df <- all_window_species %>%
    mutate(
      protocol       = factor(window_id, levels = protocol_keep),
      protocol_bin   = as.integer(protocol == "SNAP_EU_BUFFER"),
      ds_sp          = interaction(dataset, species, drop = TRUE),
      l_trapdays     = log1p(trap_days_window),
      l_nsites       = log1p(n_sites),
      s_latitude     = as.numeric(scale(latitude)),
      s_trap_array   = as.numeric(scale(log1p(trap_array))),
      w_truth        = log1p(trap_days_full)
    ) %>%
    filter(!is.na(protocol),
           is.finite(mse_p_tte),
           is.finite(mse_log_rate),
           is.finite(mse_spatial_cov))
  
  # Aggregate to one row per dataset × species × protocol
  df_sum <- df %>%
    group_by(dataset, species, ds_sp, protocol, protocol_bin) %>%
    summarise(
      mse_p_tte       = median(mse_p_tte,       na.rm = TRUE),
      mse_log_rate    = median(mse_log_rate,    na.rm = TRUE),
      mse_spatial_cov = median(mse_spatial_cov, na.rm = TRUE),
      l_trapdays      = median(l_trapdays,      na.rm = TRUE),
      l_nsites        = median(l_nsites,        na.rm = TRUE),
      s_latitude      = median(s_latitude,      na.rm = TRUE),
      s_trap_array    = median(s_trap_array,    na.rm = TRUE),
      w_truth         = median(w_truth,         na.rm = TRUE),
      .groups = "drop"
    )
  
  # Species richness (portfolio effect covariate)
  richness <- df %>%
    group_by(dataset, protocol) %>%
    summarise(n_species = n_distinct(species), .groups = "drop")
  
  df_sum %>%
    left_join(richness, by = c("dataset", "protocol")) %>%
    mutate(l_n_species = log1p(n_species))
}


############################################################
# 2) MODEL REGISTRY — DETECTION MSE ONLY
############################################################

model_specs <- list(
  
  # ════════════════════════════════════════════════════════
  # TTE DETECTION PROBABILITY MSE MODELS
  # ════════════════════════════════════════════════════════
  
  mse_tte_base = list(
    resp   = "mse_p_tte",
    family = Gamma(link = "log"),
    formula = mse_p_tte ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 | species)
  ),
  
  mse_tte_spxprot = list(
    resp   = "mse_p_tte",
    family = Gamma(link = "log"),
    formula = mse_p_tte ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 + protocol_bin || species)
  ),
  
  mse_tte_guild = list(
    resp   = "mse_p_tte",
    family = Gamma(link = "log"),
    formula = mse_p_tte ~ 1 + protocol_bin * guild_major +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 | species)
  ),
  
  mse_tte_richness_xprot = list(
    resp   = "mse_p_tte",
    family = Gamma(link = "log"),
    formula = mse_p_tte ~ 1 + protocol_bin * l_n_species +
      l_trapdays + l_nsites + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  
  # ════════════════════════════════════════════════════════
  # EVENT RATE MSE MODELS (LOG SCALE)
  # ════════════════════════════════════════════════════════
  
  mse_lrr_base = list(
    resp   = "mse_log_rate",
    family = Gamma(link = "log"),
    formula = mse_log_rate ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 | species)
  ),
  
  mse_lrr_spxprot = list(
    resp   = "mse_log_rate",
    family = Gamma(link = "log"),
    formula = mse_log_rate ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 + protocol_bin || species)
  ),
  
  mse_lrr_guild = list(
    resp   = "mse_log_rate",
    family = Gamma(link = "log"),
    formula = mse_log_rate ~ 1 + protocol_bin * guild_major +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 | species)
  ),
  
  mse_lrr_richness_xprot = list(
    resp   = "mse_log_rate",
    family = Gamma(link = "log"),
    formula = mse_log_rate ~ 1 + protocol_bin * l_n_species +
      l_trapdays + l_nsites + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  
  # ════════════════════════════════════════════════════════
  # SPATIAL COVERAGE MSE MODELS
  # ════════════════════════════════════════════════════════
  
  mse_scov_base = list(
    resp   = "mse_spatial_cov",
    family = Gamma(link = "log"),
    formula = mse_spatial_cov ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 | species)
  ),
  
  mse_scov_spxprot = list(
    resp   = "mse_spatial_cov",
    family = Gamma(link = "log"),
    formula = mse_spatial_cov ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 + protocol_bin || species)
  ),
  
  mse_scov_guild = list(
    resp   = "mse_spatial_cov",
    family = Gamma(link = "log"),
    formula = mse_spatial_cov ~ 1 + protocol_bin * guild_major +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 | species)
  ),
  
  mse_scov_richness_xprot = list(
    resp   = "mse_spatial_cov",
    family = Gamma(link = "log"),
    formula = mse_spatial_cov ~ 1 + protocol_bin * l_n_species +
      l_trapdays + l_nsites + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  )
)


############################################################
# 3) GENERIC FITTER
############################################################

fit_from_spec <- function(df_mse, spec_name,
                          iter          = 4000,
                          warmup        = 2000,
                          chains        = 4,
                          cores         = 4,
                          adapt_delta   = 0.98,
                          max_treedepth = 12,
                          seed          = 1) {
  
  sp <- model_specs[[spec_name]]
  stopifnot(!is.null(sp))
  
  brm(
    formula   = sp$formula,
    data      = df_mse,
    family    = sp$family,
    iter      = iter,
    warmup    = warmup,
    chains    = chains,
    cores     = cores,
    seed      = seed,
    control   = list(adapt_delta   = adapt_delta,
                     max_treedepth = max_treedepth),
    save_pars = save_pars(all = TRUE)
  )
}


############################################################
# 4) HELPER FUNCTIONS
############################################################

# Compute effect summaries for fixed effects
compute_effect_summary <- function(fit, model_name) {
  
  draws <- as_draws_df(fit)
  pars  <- names(draws)[grepl("^b_", names(draws))]
  
  out <- lapply(pars, function(p) {
    x <- draws[[p]]
    tibble(
      model     = model_name,
      parameter = p,
      median    = median(x),
      q2.5      = quantile(x, 0.025),
      q97.5     = quantile(x, 0.975),
      p_gt0     = mean(x > 0),
      p_lt0     = mean(x < 0)
    )
  })
  
  bind_rows(out)
}

# Population-level protocol comparison (BUFFER vs CORE) on response scale
p_buffer_better_mse <- function(fit, df) {
  
  nd <- tibble(
    protocol_bin = c(0L, 1L),
    l_trapdays   = median(df$l_trapdays),
    l_nsites     = median(df$l_nsites),
    s_latitude   = median(df$s_latitude),
    s_trap_array = median(df$s_trap_array),
    l_n_species  = median(df$l_n_species)
  )
  
  mu <- posterior_epred(fit, newdata = nd, re_formula = NA)
  if (length(dim(mu)) == 3) mu <- mu[, , 1]
  
  mu_core   <- mu[, 1]
  mu_buffer <- mu[, 2]
  
  # Lower MSE is better
  mean(mu_buffer < mu_core)
}

# LOO-based model selection helper
select_best_model <- function(loo_list, delta_threshold = 4) {
  comp <- loo_compare(loo_list)
  
  comp_df <- as.data.frame(comp) %>%
    tibble::rownames_to_column("model") %>%
    tibble::as_tibble()
  
  best_model <- comp_df$model[1]
  
  if (nrow(comp_df) > 1) {
    delta_loo <- abs(comp_df$elpd_diff[2])
    if (delta_loo < delta_threshold) {
      message("  → Top 2 models within ΔLOO threshold (", round(delta_loo, 1),
              " < ", delta_threshold, ")")
      message("  → Consider parsimony: choose simpler model")
    } else {
      message("  → Clear winner: ", best_model,
              " (ΔLOO = ", round(delta_loo, 1), ")")
    }
  }
  
  list(
    best_model     = best_model,
    comparison     = comp_df,
    recommendation = if (exists("delta_loo") && delta_loo < delta_threshold)
      "Consider simpler model (tie)" else "Use best LOO model"
  )
}


############################################################
# 5) WORKFLOW
############################################################

# ---- Prepare data ----
df_mse <- prep_mse_summaries(all_window_species)

# ---- Join trait data ----
traits <- read_csv("species_meta.csv", show_col_types = FALSE) %>%
  select(species, guild_major, family)

df_mse <- df_mse %>%
  left_join(traits, by = "species") %>%
  mutate(
    guild_major = factor(guild_major),
    family      = factor(family)
  )

# ---- Verify sample size ----
message("MSE detection models sample size: N = ", nrow(df_mse))
message("  Datasets: ", n_distinct(df_mse$dataset))
message("  Species: ",  n_distinct(df_mse$species))
message("  Protocols: ", paste(unique(df_mse$protocol), collapse = ", "))


# ════════════════════════════════════════════════════════
# FIT TTE MSE MODELS
# ════════════════════════════════════════════════════════

fit_mse_tte_base           <- fit_from_spec(df_mse, "mse_tte_base")
pp_check(fit_mse_tte_base, ndraws = 100)

fit_mse_tte_spx            <- fit_from_spec(df_mse, "mse_tte_spxprot")
pp_check(fit_mse_tte_spx, ndraws = 100)

fit_mse_tte_guild          <- fit_from_spec(df_mse, "mse_tte_guild")
pp_check(fit_mse_tte_guild, ndraws = 100)

fit_mse_tte_richness_xprot <- fit_from_spec(df_mse, "mse_tte_richness_xprot")
pp_check(fit_mse_tte_richness_xprot, ndraws = 100)

# Compare models
loo_mse_tte_base           <- loo(fit_mse_tte_base,           moment_match = TRUE)
loo_mse_tte_spx            <- loo(fit_mse_tte_spx,            moment_match = TRUE)
loo_mse_tte_guild          <- loo(fit_mse_tte_guild,          moment_match = TRUE)
loo_mse_tte_richness_xprot <- loo(fit_mse_tte_richness_xprot, moment_match = TRUE)

message("\n═══ TTE MSE MODEL SELECTION ═══")
tte_mse_selection <- select_best_model(
  list(
    tte_base           = loo_mse_tte_base,
    tte_spx            = loo_mse_tte_spx,
    tte_guild          = loo_mse_tte_guild,
    tte_richness_xprot = loo_mse_tte_richness_xprot
  )
)

# Summarize base model
summary(fit_mse_tte_base)
draws_mse_tte <- as_draws_df(fit_mse_tte_base)

# Protocol effect (on log scale)
mean(draws_mse_tte$b_protocol_bin < 0)  # P(BUFFER reduces MSE)

# Response-scale comparison
p_buffer_better_mse(fit_mse_tte_base, df_mse)

# Covariate effects
cov_names <- c("l_trapdays","l_nsites","s_latitude","s_trap_array","l_n_species")
summ_cov_mse_tte <- map_df(cov_names, ~{
  col <- paste0("b_", .x)
  d   <- draws_mse_tte[[col]]
  tibble(
    term  = .x,
    mean  = mean(d),
    q025  = quantile(d, 0.025),
    q975  = quantile(d, 0.975),
    p_neg = mean(d < 0),
    p_pos = mean(d > 0)
  )
})
print(summ_cov_mse_tte)


# ════════════════════════════════════════════════════════
# FIT LOG-RATE MSE MODELS
# ════════════════════════════════════════════════════════

fit_mse_lrr_base           <- fit_from_spec(df_mse, "mse_lrr_base")
pp_check(fit_mse_lrr_base, ndraws = 100)

fit_mse_lrr_spx            <- fit_from_spec(df_mse, "mse_lrr_spxprot")
pp_check(fit_mse_lrr_spx, ndraws = 100)

fit_mse_lrr_guild          <- fit_from_spec(df_mse, "mse_lrr_guild")
pp_check(fit_mse_lrr_guild, ndraws = 100)

fit_mse_lrr_richness_xprot <- fit_from_spec(df_mse, "mse_lrr_richness_xprot")
pp_check(fit_mse_lrr_richness_xprot, ndraws = 100)

# Compare models
loo_mse_lrr_base           <- loo(fit_mse_lrr_base,           moment_match = TRUE)
loo_mse_lrr_spx            <- loo(fit_mse_lrr_spx,            moment_match = TRUE)
loo_mse_lrr_guild          <- loo(fit_mse_lrr_guild,          moment_match = TRUE)
loo_mse_lrr_richness_xprot <- loo(fit_mse_lrr_richness_xprot, moment_match = TRUE)

message("\n═══ LOG-RATE MSE MODEL SELECTION ═══")
lrr_mse_selection <- select_best_model(
  list(
    lrr_base           = loo_mse_lrr_base,
    lrr_spx            = loo_mse_lrr_spx,
    lrr_guild          = loo_mse_lrr_guild,
    lrr_richness_xprot = loo_mse_lrr_richness_xprot
  )
)

summary(fit_mse_lrr_base)
draws_mse_lrr <- as_draws_df(fit_mse_lrr_base)

mean(draws_mse_lrr$b_protocol_bin < 0)
p_buffer_better_mse(fit_mse_lrr_base, df_mse)

summ_cov_mse_lrr <- map_df(cov_names, ~{
  col <- paste0("b_", .x)
  d   <- draws_mse_lrr[[col]]
  tibble(
    term  = .x,
    mean  = mean(d),
    q025  = quantile(d, 0.025),
    q975  = quantile(d, 0.975),
    p_neg = mean(d < 0),
    p_pos = mean(d > 0)
  )
})
print(summ_cov_mse_lrr)


# ════════════════════════════════════════════════════════
# FIT SPATIAL COVERAGE MSE MODELS
# ════════════════════════════════════════════════════════

fit_mse_scov_base           <- fit_from_spec(df_mse, "mse_scov_base")
pp_check(fit_mse_scov_base, ndraws = 100)

fit_mse_scov_spx            <- fit_from_spec(df_mse, "mse_scov_spxprot")
pp_check(fit_mse_scov_spx, ndraws = 100)

fit_mse_scov_guild          <- fit_from_spec(df_mse, "mse_scov_guild")
pp_check(fit_mse_scov_guild, ndraws = 100)

fit_mse_scov_richness_xprot <- fit_from_spec(df_mse, "mse_scov_richness_xprot")
pp_check(fit_mse_scov_richness_xprot, ndraws = 100)

# Compare models
loo_mse_scov_base           <- loo(fit_mse_scov_base,           moment_match = TRUE)
loo_mse_scov_spx            <- loo(fit_mse_scov_spx,            moment_match = TRUE)
loo_mse_scov_guild          <- loo(fit_mse_scov_guild,          moment_match = TRUE)
loo_mse_scov_richness_xprot <- loo(fit_mse_scov_richness_xprot, moment_match = TRUE)

message("\n═══ SPATIAL COVERAGE MSE MODEL SELECTION ═══")
scov_mse_selection <- select_best_model(
  list(
    scov_base           = loo_mse_scov_base,
    scov_spx            = loo_mse_scov_spx,
    scov_guild          = loo_mse_scov_guild,
    scov_richness_xprot = loo_mse_scov_richness_xprot
  )
)

summary(fit_mse_scov_base)
draws_mse_scov <- as_draws_df(fit_mse_scov_base)

mean(draws_mse_scov$b_protocol_bin < 0)
p_buffer_better_mse(fit_mse_scov_base, df_mse)

summ_cov_mse_scov <- map_df(cov_names, ~{
  col <- paste0("b_", .x)
  d   <- draws_mse_scov[[col]]
  tibble(
    term  = .x,
    mean  = mean(d),
    q025  = quantile(d, 0.025),
    q975  = quantile(d, 0.975),
    p_neg = mean(d < 0),
    p_pos = mean(d > 0)
  )
})
print(summ_cov_mse_scov)


############################################################
# 6) CROSS-METRIC COMPARISON
############################################################

# Combine all effect summaries across the three base models
all_mse_effects <- bind_rows(
  compute_effect_summary(fit_mse_tte_base,  "mse_tte_base"),
  compute_effect_summary(fit_mse_lrr_base,  "mse_lrr_base"),
  compute_effect_summary(fit_mse_scov_base, "mse_scov_base")
) %>%
  mutate(
    parameter = gsub("b_", "", parameter),
    interpretation = case_when(
      p_lt0 > 0.95 ~ "Strong negative effect (reduces MSE)",
      p_gt0 > 0.95 ~ "Strong positive effect (increases MSE)",
      TRUE         ~ "Uncertain"
    )
  )

print(all_mse_effects)

# Protocol effectiveness across metrics
protocol_mse_effects <- all_mse_effects %>%
  filter(parameter == "protocol_bin") %>%
  select(model, median, q2.5, q97.5, p_lt0)
print(protocol_mse_effects)

# Richness (portfolio) effects across metrics
richness_mse_effects <- all_mse_effects %>%
  filter(parameter == "l_n_species") %>%
  select(model, median, q2.5, q97.5, p_lt0)
print(richness_mse_effects)


############################################################
# 7) COMPARISON WITH ABSOLUTE DEVIATION MODELS
############################################################
#
# INTERPRETATION GUIDE:
#
# If MSE models and absolute deviation models agree:
#   → Strong evidence for protocol preference
#
# If MSE shows BUFFER better but |Δ| shows tie:
#   → BUFFER reduces variance despite similar bias
#   → Prefer BUFFER for stability
#
# If |Δ| shows BUFFER better but MSE shows tie:
#   → BUFFER reduces bias but increases variance
#   → Trade-off depends on context
#
# If they disagree on direction:
#   → Investigate bias-variance decomposition
#   → Check for outliers or heterogeneity
#
# Portfolio effect (l_n_species):
#   → Should show negative coefficient in both MSE and |Δ| models
#   → If effect stronger in MSE: richness reduces variance more than bias
#   → If effect stronger in |Δ|: richness reduces bias more than variance
#
############################################################

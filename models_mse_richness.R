############################################################
#   PROTOCOL EVALUATION: MSE-BASED ANALYSIS
#   SPECIES RICHNESS METRIC (mse_sr_raref)
#
#   OBJECTIVE: Compare protocols using bias-variance tradeoff
#   via Mean Squared Error (MSE) for rarefied species richness
#
#   MSE = Bias² + Variance
#   where:
#     - Bias     = (sr_raref_window - sr_raref_full)²
#     - Variance = sr_raref_var  (analytical rarefaction variance,
#                                 Colwell 2012)
#
#   METRIC:
#   - mse_sr_raref : d_sr_raref² + sr_raref_var  [0,∞)
#     Jointly captures bias relative to full-year richness AND
#     sampling imprecision from the rarefaction estimator.
#
#   NOTE: absolute deviation richness models (abs_d_sr_raref,
#         prop_sr_full) are in
#         models_abs_richness_PARALLEL.R
#
############################################################

# ──────────────────────────────────────────────────────────
# OVERVIEW OF MODELS
# ──────────────────────────────────────────────────────────
#
# RESPONSE VARIABLE:
#   - mse_sr_raref : bias² + rarefaction variance  [0,∞)
#
# FIXED EFFECTS:
#   - protocol_bin            : BUFFER=1, CORE=0 (focal predictor)
#   - log1p(trap_days_window) : Sampling effort (trap-days)
#   - log1p(n_sites)          : Number of camera sites
#   - scale(latitude)         : Geographic latitude
#   - scale(log1p(trap_array)): Spatial extent of camera array
#   - log1p(n_species)        : Species richness (portfolio effect)
#     OR s_mean_p_full        : Community mean detectability (alternative)
#
# RANDOM EFFECTS:
#   - (1 | dataset)                 : Dataset-level intercept
#   - (1 + protocol_bin || dataset) : Dataset-specific protocol slope
#
# STATISTICAL FAMILY:
#   - Gaussian  : MSE is non-negative by construction but can be
#                 approximately Gaussian after the bias² + var sum;
#                 check pp_check for fit quality
#
############################################################

# ---- Packages ----
library(dplyr)
library(brms)
library(readr)
library(tibble)
library(posterior)
library(purrr)


############################################################
# 1) DATA PREPARATION
############################################################

# Reuses the same helper as models_abs_richness_PARALLEL.R to build
# dataset-level covariates from the species-level data.
prep_protocol_summaries_spp <- function(all_window_species,
                                        protocol_keep = c("SNAP_EU_CORE","SNAP_EU_BUFFER")) {
  
  req <- c("dataset","species","window_id",
           "p_tte","log_rate","spatial_cov",
           "trap_days_window","n_sites","latitude","trap_array",
           "trap_days_full","n_sites_full")
  stopifnot(all(req %in% names(all_window_species)))
  
  df <- all_window_species %>%
    mutate(
      protocol     = factor(window_id, levels = protocol_keep),
      protocol_bin = as.integer(protocol == "SNAP_EU_BUFFER"),
      ds_sp        = interaction(dataset, species, drop = TRUE),
      l_trapdays   = log1p(trap_days_window),
      l_nsites     = log1p(n_sites),
      s_latitude   = as.numeric(scale(latitude)),
      s_trap_array = as.numeric(scale(log1p(trap_array))),
      abs_d_p_tte       = abs(d_p_tte),
      abs_d_log_rate    = abs(d_log_rate),
      abs_d_spatial_cov = abs(d_spatial_cov)
    ) %>%
    filter(!is.na(protocol),
           is.finite(abs_d_p_tte),
           is.finite(abs_d_log_rate),
           is.finite(abs_d_spatial_cov))
  
  richness <- df %>%
    group_by(dataset, protocol) %>%
    summarise(n_species = n_distinct(species), .groups = "drop")
  
  p_tte_reference <- all_window_species %>%
    group_by(dataset, species) %>%
    slice_max(order_by = trap_days_window, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(dataset, species, p_tte_full_ref = p_tte)
  
  community_p_full <- p_tte_reference %>%
    group_by(dataset) %>%
    summarise(mean_p_full = mean(p_tte_full_ref, na.rm = TRUE), .groups = "drop") %>%
    mutate(s_mean_p_full = as.numeric(scale(mean_p_full)))
  
  df %>%
    group_by(dataset, species, ds_sp, protocol, protocol_bin) %>%
    summarise(
      l_trapdays   = median(l_trapdays,   na.rm = TRUE),
      l_nsites     = median(l_nsites,     na.rm = TRUE),
      s_latitude   = median(s_latitude,   na.rm = TRUE),
      s_trap_array = median(s_trap_array, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(richness,        by = c("dataset", "protocol")) %>%
    left_join(community_p_full, by = "dataset") %>%
    mutate(l_n_species = log1p(n_species))
}


############################################################
# 2) MODEL REGISTRY — RICHNESS MSE ONLY
############################################################

model_specs <- list(
  
  # ── Base model: dataset random intercept ──────────────────────────────────
  rich_mse_base = list(
    resp   = "mse_sr_raref",
    family = gaussian(),
    formula = mse_sr_raref ~ 1 + protocol_bin +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | dataset)
  ),
  
  # ── Dataset-specific protocol slope ───────────────────────────────────────
  rich_mse_ds_slope = list(
    resp   = "mse_sr_raref",
    family = gaussian(),
    formula = mse_sr_raref ~ 1 + protocol_bin +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 + protocol_bin || dataset)
  ),
  
  # ── Community mean detectability instead of l_n_species ───────────────────
  # Tests whether assemblage-level detectability explains MSE better than
  # raw richness (portfolio effect).
  rich_mse_mean_p_full = list(
    resp   = "mse_sr_raref",
    family = gaussian(),
    formula = mse_sr_raref ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_mean_p_full + s_latitude + s_trap_array +
      (1 | dataset)
  )
)


############################################################
# 3) GENERIC FITTER
############################################################

fit_from_spec <- function(df_rich, spec_name,
                          iter          = 4000,
                          warmup        = 2000,
                          chains        = 4,
                          cores         = 4,
                          adapt_delta   = 0.98,
                          max_treedepth = 12,
                          seed          = 1) {
  
  sp <- model_specs[[spec_name]]
  stopifnot(!is.null(sp))
  
  # Per-model adapt_delta override
  adapt_delta   <- if (!is.null(sp$adapt_delta))   sp$adapt_delta   else adapt_delta
  max_treedepth <- if (!is.null(sp$max_treedepth)) sp$max_treedepth else max_treedepth
  
  brm(
    formula   = sp$formula,
    data      = df_rich,
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
# 4) DATA PREPARATION
############################################################

# Build species-level summaries to derive dataset-level covariates
df_sum_tmp <- prep_protocol_summaries_spp(all_window_species)

dataset_covs <- df_sum_tmp %>%
  group_by(dataset) %>%
  summarise(
    s_latitude    = median(s_latitude,    na.rm = TRUE),
    s_trap_array  = median(s_trap_array,  na.rm = TRUE),
    l_n_species   = median(l_n_species,   na.rm = TRUE),
    s_mean_p_full = first(s_mean_p_full), # dataset-level; constant within dataset
    .groups = "drop"
  )

df_rich <- all_window_richness %>%
  filter(window_id %in% c("SNAP_EU_CORE", "SNAP_EU_BUFFER")) %>%
  mutate(
    protocol_bin = as.integer(window_id == "SNAP_EU_BUFFER"),
    l_trapdays   = log1p(trap_days_window),
    l_nsites     = log1p(n_sites)
  ) %>%
  left_join(dataset_covs, by = "dataset") %>%
  mutate(
    # mse_sr_raref = d_sr_raref² + sr_raref_var  (already computed in Full1.R)
    mse_sr_raref  = mse_sr_raref,
    # Richness-only datasets (no TTE data) get grand mean (= 0 on scaled)
    s_mean_p_full = if_else(is.na(s_mean_p_full), 0, s_mean_p_full)
  ) %>%
  filter(is.finite(mse_sr_raref))

message("Richness MSE models sample size: N = ", nrow(df_rich))
message("  Datasets: ", n_distinct(df_rich$dataset))
message("  Protocols: ", paste(unique(df_rich$window_id), collapse = ", "))


############################################################
# 5) FIT MODELS
############################################################

fit_rich_mse_base        <- fit_from_spec(df_rich, "rich_mse_base")
pp_check(fit_rich_mse_base, ndraws = 100)

fit_rich_mse_ds_slope    <- fit_from_spec(df_rich, "rich_mse_ds_slope")
pp_check(fit_rich_mse_ds_slope, ndraws = 100)

fit_rich_mse_mean_p_full <- fit_from_spec(df_rich, "rich_mse_mean_p_full")
pp_check(fit_rich_mse_mean_p_full, ndraws = 100)


############################################################
# 6) LOO MODEL COMPARISON
############################################################

loo_rich_mse_base        <- loo(fit_rich_mse_base,        moment_match = TRUE)
loo_rich_mse_ds_slope    <- loo(fit_rich_mse_ds_slope,    moment_match = TRUE)
loo_rich_mse_mean_p_full <- loo(fit_rich_mse_mean_p_full, moment_match = TRUE)

message("\n═══ RICHNESS MSE MODEL SELECTION ═══")
rich_mse_selection <- select_best_model(
  list(
    rich_mse_base        = loo_rich_mse_base,
    rich_mse_ds_slope    = loo_rich_mse_ds_slope,
    rich_mse_mean_p_full = loo_rich_mse_mean_p_full
  )
)


############################################################
# 7) SUMMARIZE RESULTS
############################################################

summary(fit_rich_mse_base)
draws_rich_mse <- as_draws_df(fit_rich_mse_base)

# Protocol effect
mean(draws_rich_mse$b_protocol_bin < 0)  # P(BUFFER reduces MSE for richness)

# Covariate effects
cov_names <- c("l_trapdays","l_nsites","s_latitude","s_trap_array","l_n_species")
summ_cov_rich_mse <- map_df(cov_names, ~{
  col <- paste0("b_", .x)
  d   <- draws_rich_mse[[col]]
  tibble(
    term  = .x,
    mean  = mean(d),
    q025  = quantile(d, 0.025),
    q975  = quantile(d, 0.975),
    p_neg = mean(d < 0),
    p_pos = mean(d > 0)
  )
})
print(summ_cov_rich_mse)

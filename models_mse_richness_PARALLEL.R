############################################################
#   PROTOCOL EVALUATION: MSE-BASED ANALYSIS
#   SPECIES RICHNESS METRIC (mse_sr_raref)
#   PARALLEL FITTING VERSION
#
#   Machine:  16 physical / 22 logical cores
#   Strategy: 3 models × 4 cores = 12 cores simultaneously
#             (3 models fit concurrently in one batch)
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
#   Models (dataset-level, Gaussian family): 3 total
#     - rich_mse_base        : dataset random intercept
#     - rich_mse_ds_slope    : dataset-specific protocol slope
#     - rich_mse_mean_p_full : community detectability instead of richness
#
#   NOTE: absolute deviation richness models (abs_d_sr_raref,
#         prop_sr_full) are in
#         models_abs_richness_PARALLEL.R
#
#   Requires: future, furrr  (install if missing)
#     install.packages(c("future", "furrr"))
#
############################################################

library(dplyr)
library(brms)
library(readr)
library(tibble)
library(posterior)
library(purrr)
library(future)
library(furrr)


############################################################
# PARALLELISM SETTINGS
############################################################

N_WORKERS       <- 3   # 3 models fit simultaneously (one per model)
CORES_PER_MODEL <- 4   # Stan chains per model
# 3 × 4 = 12 cores total — well within 16-core budget

plan(multisession, workers = N_WORKERS)

# Confirm plan
message("Parallel plan: ", class(plan())[1])
message("Workers: ", N_WORKERS, " × ", CORES_PER_MODEL, " cores = ",
        N_WORKERS * CORES_PER_MODEL, " cores total")


############################################################
# 1) DATA PREPARATION
############################################################

# Reuses the same helper as models_abs_richness_PARALLEL.R to build
# dataset-level covariates from the species-level data.
prep_protocol_summaries_spp <- function(all_window_species,
                                        protocol_keep = c("SNAP_EU_CORE","SNAP_EU_BUFFER")) {

  req <- c("dataset","species","window_id",
           "lambda","log_rate","matched_rate",
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
      abs_d_lambda       = abs(d_lambda),
      abs_d_rate         = abs(d_rate),
      abs_d_matched_rate = abs(d_matched_rate)
    ) %>%
    filter(!is.na(protocol),
           is.finite(abs_d_lambda),
           is.finite(abs_d_rate),
           is.finite(abs_d_matched_rate))

  richness <- df %>%
    group_by(dataset, protocol) %>%
    summarise(n_species = n_distinct(species), .groups = "drop")

  lambda_reference <- all_window_species %>%
    group_by(dataset, species) %>%
    slice_max(order_by = trap_days_window, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(dataset, species, lambda_full_ref = lambda)

  community_lambda_full <- lambda_reference %>%
    group_by(dataset) %>%
    summarise(mean_lambda_full = mean(lambda_full_ref, na.rm = TRUE), .groups = "drop") %>%
    mutate(s_mean_lambda_full = as.numeric(scale(mean_lambda_full)))

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
    left_join(community_lambda_full, by = "dataset") %>%
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
  rich_mse_mean_lambda_full = list(
    resp   = "mse_sr_raref",
    family = gaussian(),
    formula = mse_sr_raref ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_mean_lambda_full + s_latitude + s_trap_array +
      (1 | dataset)
  )
)


############################################################
# 3) GENERIC FITTER
############################################################

fit_one <- function(spec_name, df_rich,
                    iter          = 4000,
                    warmup        = 2000,
                    chains        = 4,
                    cores         = CORES_PER_MODEL,
                    adapt_delta   = 0.98,
                    max_treedepth = 12) {

  library(brms)

  sp <- model_specs[[spec_name]]
  stopifnot(!is.null(sp))

  # Per-spec overrides
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
    seed      = 1,
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
    s_mean_lambda_full = first(s_mean_lambda_full), # dataset-level; constant within dataset
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
    s_mean_lambda_full = if_else(is.na(s_mean_lambda_full), 0, s_mean_lambda_full)
  ) %>%
  filter(is.finite(mse_sr_raref))

message("Richness MSE models sample size: N = ", nrow(df_rich))
message("  Datasets: ", n_distinct(df_rich$dataset))
message("  Protocols: ", paste(unique(df_rich$window_id), collapse = ", "))


############################################################
# 5) FIT ALL MODELS IN PARALLEL
############################################################

all_spec_names <- names(model_specs)

message("\n─── Starting parallel fit: ", length(all_spec_names),
        " richness MSE models across ", N_WORKERS, " workers ───")
message("Started at: ", Sys.time())

t_start <- proc.time()

fits_list <- future_map(
  all_spec_names,
  \(nm) fit_one(nm, df_rich),
  .options = furrr_options(seed = TRUE),   # safe RNG across workers
  .progress = TRUE                          # progress bar in console
)

names(fits_list) <- all_spec_names

t_elapsed <- proc.time() - t_start
message("\n─── All models fitted ───")
message("Total elapsed time: ", round(t_elapsed["elapsed"] / 60, 1), " minutes")
message("Finished at: ", Sys.time())


############################################################
# 6) SAVE ALL FITS
############################################################

saveRDS(fits_list, "fits_mse_richness_all.rds")
message("Fits saved to fits_mse_richness_all.rds")


############################################################
# 7) LOAD FITS  (re-entry point — skip §4–6 after first run)
############################################################

fits_list      <- readRDS("fits_mse_richness_all.rds")
all_spec_names <- names(fits_list)


############################################################
# 8) DIAGNOSTICS & POSTERIOR PREDICTIVE CHECKS
############################################################

check_diagnostics <- function(fit, name) {
  s  <- rstan::get_sampler_params(fit$fit, inc_warmup = FALSE)
  nd <- sum(sapply(s, function(x) sum(x[, "divergent__"])))
  rh <- max(brms::rhat(fit), na.rm = TRUE)
  message(sprintf("  %-35s | divergences: %d | max Rhat: %.4f", name, nd, rh))
}

message("\n─── Diagnostics ───")
iwalk(fits_list, check_diagnostics)

message("\n─── Posterior predictive checks ───")
library(ggplot2)

# Helper: run all three checks for one model and print them
pp_check_full <- function(fit, label) {
  message("  ", label)
  print(pp_check(fit, ndraws = 100)                                     + ggtitle(paste(label, "— density")))
  print(pp_check(fit, ndraws = 100, type = "stat",    stat = "mean")   + ggtitle(paste(label, "— mean")))
  print(pp_check(fit, ndraws = 100, type = "stat_2d", stat = c("mean", "sd")) + ggtitle(paste(label, "— mean vs SD")))
}

pp_check_full(fits_list[["rich_mse_base"]], "Richness MSE base")


############################################################
# 9) LOO-CV
############################################################

message("\n─── LOO-CV ───")

# All 3 models fit in one batch — use 3 workers
plan(multisession, workers = length(fits_list))
message("LOO workers: ", length(fits_list))

t_loo <- proc.time()

loo_list <- future_map(
  fits_list,
  function(fit) {
    library(brms)
    loo(fit, moment_match = TRUE)
  },
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)
names(loo_list) <- all_spec_names

message("LOO finished in ", round((proc.time() - t_loo)["elapsed"] / 60, 1), " minutes")

saveRDS(loo_list, "loo_mse_richness_all.rds")


############################################################
# 10) MODEL SELECTION
############################################################

message("\n═══ RICHNESS MSE MODEL SELECTION ═══")
rich_mse_selection <- select_best_model(
  list(
    rich_mse_base        = loo_list[["rich_mse_base"]],
    rich_mse_ds_slope    = loo_list[["rich_mse_ds_slope"]],
    rich_mse_mean_lambda_full = loo_list[["rich_mse_mean_lambda_full"]]
  )
)


############################################################
# 11) SUMMARIZE RESULTS (BASE MODEL)
############################################################

summary(fits_list[["rich_mse_base"]])
draws_rich_mse <- as_draws_df(fits_list[["rich_mse_base"]])

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


############################################################
# 12) RESET PLAN WHEN DONE
############################################################

plan(sequential)
message("Parallel plan reset to sequential.")


############################################################
# 13) CLEAN ENVIRONMENT
#     Keep only the fits list, LOO list, and modelling data.
#     Remove all intermediate/helper objects.
############################################################

rm(list = setdiff(
  ls(),
  c("fits_list", "loo_list", "all_spec_names",
    "df_rich",
    # model selection results
    "rich_mse_selection",
    # summary tables
    "draws_rich_mse", "summ_cov_rich_mse",
    # upstream data (keep for downstream scripts)
    "all_window_species", "all_window_richness", "all_dropped_species",
    "window_grid", "ds_paths", "spp_keep", "anchors")
))
gc()
message("Environment cleaned. Retained: fits_list, loo_list, df_rich, selection results.")

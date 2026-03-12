############################################################
#   PROTOCOL EVALUATION: MSE-BASED ANALYSIS
#   SPECIES-LEVEL DETECTION METRICS (LAMBDA / LRR / MRATE)
#   PARALLEL FITTING VERSION
#
#   Machine:  16 physical / 22 logical cores
#   Strategy: 4 models × 4 cores = 16 cores simultaneously
#             Wall-clock time ≈ ceil(N models / 4) × per-model time
#
#   MSE = Bias² + Variance
#   where:
#     - Bias = (estimate - truth)²
#     - Variance = SE²
#
#   METRICS ANALYZED:
#   1. MSE of daily detection rate (lambda)               (mse_lambda)
#   2. MSE of event rate (raw difference)                  (mse_rate)
#   3. MSE of matched-camera daily detection rate         (mse_matched_rate)
#
#   Models: 4 variants × 3 metrics = 12 total
#     - base / spxprot / guild / richness_xprot
#
#   NOTE: richness MSE models (mse_sr_raref) are in
#         models_mse_richness_PARALLEL.R
#
#   Requires: future, furrr  (install if missing)
#     install.packages(c("future", "furrr"))
#
############################################################

library(dplyr)
library(brms)
library(readr)
library(tidyr)
library(tibble)
library(posterior)
library(purrr)
library(future)
library(furrr)


############################################################
# PARALLELISM SETTINGS
############################################################

N_WORKERS       <- 4   # models fitting simultaneously
CORES_PER_MODEL <- 4   # Stan chains per model
# N_WORKERS × CORES_PER_MODEL should not exceed physical cores (16)

plan(multisession, workers = N_WORKERS)

# Confirm plan
message("Parallel plan: ", class(plan())[1])
message("Workers: ", N_WORKERS, " × ", CORES_PER_MODEL, " cores = ",
        N_WORKERS * CORES_PER_MODEL, " cores total")


############################################################
# 1) DATA PREPARATION
############################################################

prep_mse_summaries <- function(all_window_species,
                                protocol_keep = c("SNAP_EU_CORE","SNAP_EU_BUFFER")) {

  req <- c("dataset","species","window_id",
           "mse_lambda","mse_rate","mse_matched_rate",
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
           is.finite(mse_lambda),
           is.finite(mse_rate),
           is.finite(mse_matched_rate))

  # Aggregate to one row per dataset × species × protocol
  df_sum <- df %>%
    group_by(dataset, species, ds_sp, protocol, protocol_bin) %>%
    summarise(
      mse_lambda       = median(mse_lambda,       na.rm = TRUE),
      mse_rate         = median(mse_rate,           na.rm = TRUE),
      mse_matched_rate = median(mse_matched_rate, na.rm = TRUE),
      l_trapdays       = median(l_trapdays,       na.rm = TRUE),
      l_nsites         = median(l_nsites,         na.rm = TRUE),
      s_latitude       = median(s_latitude,       na.rm = TRUE),
      s_trap_array     = median(s_trap_array,     na.rm = TRUE),
      w_truth          = median(w_truth,          na.rm = TRUE),
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
  # LAMBDA (DAILY DETECTION RATE) MSE MODELS
  # ════════════════════════════════════════════════════════

  mse_lambda_base = list(
    resp   = "mse_lambda",
    family = Gamma(link = "log"),
    formula = mse_lambda ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 | species)
  ),

  mse_lambda_spxprot = list(
    resp   = "mse_lambda",
    family = Gamma(link = "log"),
    formula = mse_lambda ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 + protocol_bin || species)
  ),

  mse_lambda_guild = list(
    resp   = "mse_lambda",
    family = Gamma(link = "log"),
    formula = mse_lambda ~ 1 + protocol_bin * guild_major +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 | species)
  ),

  mse_lambda_richness_xprot = list(
    resp   = "mse_lambda",
    family = Gamma(link = "log"),
    formula = mse_lambda ~ 1 + protocol_bin * l_n_species +
      l_trapdays + l_nsites + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),

  # ════════════════════════════════════════════════════════
  # EVENT RATE MSE MODELS (RAW DIFFERENCE)
  # ════════════════════════════════════════════════════════

  mse_rate_base = list(
    resp   = "mse_rate",
    family = Gamma(link = "log"),
    formula = mse_rate ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 | species)
  ),

  mse_rate_spxprot = list(
    resp   = "mse_rate",
    family = Gamma(link = "log"),
    formula = mse_rate ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 + protocol_bin || species)
  ),

  mse_rate_guild = list(
    resp   = "mse_rate",
    family = Gamma(link = "log"),
    formula = mse_rate ~ 1 + protocol_bin * guild_major +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 | species)
  ),

  mse_rate_richness_xprot = list(
    resp   = "mse_rate",
    family = Gamma(link = "log"),
    formula = mse_rate ~ 1 + protocol_bin * l_n_species +
      l_trapdays + l_nsites + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),

  # ════════════════════════════════════════════════════════
  # MATCHED-CAMERA DETECTION RATE MSE MODELS
  # ════════════════════════════════════════════════════════

  mse_mrate_base = list(
    resp   = "mse_matched_rate",
    family = Gamma(link = "log"),
    formula = mse_matched_rate ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 | species)
  ),

  mse_mrate_spxprot = list(
    resp   = "mse_matched_rate",
    family = Gamma(link = "log"),
    formula = mse_matched_rate ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 + protocol_bin || species)
  ),

  mse_mrate_guild = list(
    resp   = "mse_matched_rate",
    family = Gamma(link = "log"),
    formula = mse_matched_rate ~ 1 + protocol_bin * guild_major +
      l_trapdays + l_nsites + s_latitude + s_trap_array + l_n_species +
      (1 | ds_sp) + (1 | species)
  ),

  mse_mrate_richness_xprot = list(
    resp   = "mse_matched_rate",
    family = Gamma(link = "log"),
    formula = mse_matched_rate ~ 1 + protocol_bin * l_n_species +
      l_trapdays + l_nsites + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  )
)


############################################################
# 3) GENERIC FITTER
############################################################

fit_one <- function(spec_name, df,
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
    data      = df,
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
# 5) DATA PREPARATION
############################################################

df_mse <- prep_mse_summaries(all_window_species)

traits <- read_csv("species_meta.csv", show_col_types = FALSE) %>%
  select(species, guild_major, family)

df_mse <- df_mse %>%
  left_join(traits, by = "species") %>%
  mutate(
    guild_major = factor(guild_major),
    family      = factor(family)
  )

message("MSE detection models sample size: N = ", nrow(df_mse))
message("  Datasets: ", n_distinct(df_mse$dataset))
message("  Species: ",  n_distinct(df_mse$species))
message("  Protocols: ", paste(unique(df_mse$protocol), collapse = ", "))


############################################################
# 6) FIT ALL MODELS IN PARALLEL
############################################################

all_spec_names <- names(model_specs)

message("\n─── Starting parallel fit: ", length(all_spec_names),
        " MSE detection models across ", N_WORKERS, " workers ───")
message("Started at: ", Sys.time())

t_start <- proc.time()

fits_list <- future_map(
  all_spec_names,
  \(nm) fit_one(nm, df_mse),
  .options = furrr_options(seed = TRUE),   # safe RNG across workers
  .progress = TRUE                          # progress bar in console
)

names(fits_list) <- all_spec_names

t_elapsed <- proc.time() - t_start
message("\n─── All models fitted ───")
message("Total elapsed time: ", round(t_elapsed["elapsed"] / 60, 1), " minutes")
message("Finished at: ", Sys.time())


############################################################
# 7) SAVE ALL FITS
############################################################

saveRDS(fits_list, "fits_mse_detection_all.rds")
message("Fits saved to fits_mse_detection_all.rds")


############################################################
# 8) LOAD FITS  (re-entry point — skip §5–7 after first run)
############################################################

fits_list      <- readRDS("fits_mse_detection_all.rds")
all_spec_names <- names(fits_list)


############################################################
# 9) DIAGNOSTICS & POSTERIOR PREDICTIVE CHECKS
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

pp_check_full(fits_list[["mse_lambda_base"]], "Lambda MSE base")
pp_check_full(fits_list[["mse_lrr_base"]],   "LRR MSE base")
pp_check_full(fits_list[["mse_mrate_base"]], "MRate MSE base")


############################################################
# 10) LOO-CV
############################################################

message("\n─── LOO-CV ───")

# Scale up: one worker per model
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

saveRDS(loo_list, "loo_mse_detection_all.rds")


############################################################
# 11) MODEL SELECTION PER METRIC
############################################################

message("\n═══ LAMBDA MSE MODEL SELECTION ═══")
lambda_mse_selection <- select_best_model(
  list(
    lambda_base           = loo_list[["mse_lambda_base"]],
    lambda_spx            = loo_list[["mse_lambda_spxprot"]],
    lambda_guild          = loo_list[["mse_lambda_guild"]],
    lambda_richness_xprot = loo_list[["mse_lambda_richness_xprot"]]
  )
)

message("\n═══ LOG-RATE MSE MODEL SELECTION ═══")
lrr_mse_selection <- select_best_model(
  list(
    lrr_base           = loo_list[["mse_lrr_base"]],
    lrr_spx            = loo_list[["mse_lrr_spxprot"]],
    lrr_guild          = loo_list[["mse_lrr_guild"]],
    lrr_richness_xprot = loo_list[["mse_lrr_richness_xprot"]]
  )
)

message("\n═══ MATCHED RATE MSE MODEL SELECTION ═══")
mrate_mse_selection <- select_best_model(
  list(
    mrate_base           = loo_list[["mse_mrate_base"]],
    mrate_spx            = loo_list[["mse_mrate_spxprot"]],
    mrate_guild          = loo_list[["mse_mrate_guild"]],
    mrate_richness_xprot = loo_list[["mse_mrate_richness_xprot"]]
  )
)


############################################################
# 12) CROSS-METRIC COMPARISON (BASE MODELS)
############################################################

all_mse_effects <- bind_rows(
  compute_effect_summary(fits_list[["mse_lambda_base"]], "mse_lambda_base"),
  compute_effect_summary(fits_list[["mse_lrr_base"]],    "mse_lrr_base"),
  compute_effect_summary(fits_list[["mse_mrate_base"]],  "mse_mrate_base")
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
# 13) RESET PLAN WHEN DONE
############################################################

plan(sequential)
message("Parallel plan reset to sequential.")


############################################################
# 14) CLEAN ENVIRONMENT
#     Keep only the fits list, LOO list, and modelling data.
#     Remove all intermediate/helper objects.
############################################################

rm(list = setdiff(
  ls(),
  c("fits_list", "loo_list", "all_spec_names",
    "df_mse",
    # model selection results
    "lambda_mse_selection", "lrr_mse_selection", "mrate_mse_selection",
    # summary tables
    "all_mse_effects", "protocol_mse_effects", "richness_mse_effects",
    # upstream data (keep for downstream scripts)
    "all_window_species", "all_window_richness", "all_dropped_species",
    "window_grid", "ds_paths", "spp_keep", "anchors")
))
gc()
message("Environment cleaned. Retained: fits_list, loo_list, df_mse, selection results.")

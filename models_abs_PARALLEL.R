############################################################
#   PROTOCOL EVALUATION: ABSOLUTE DEVIATION MODELS
#   PARALLEL FITTING VERSION
#
#   Machine:  16 physical / 22 logical cores
#   Strategy: 4 models × 4 cores = 16 cores simultaneously
#             Wall-clock time ≈ ceil(N models / 4) × per-model time
#
#   Models:
#     - 14 original species-level (TTE / LRR / SCOV)
#     - 6  guild_minor variants (2 schemes × 3 metrics)
#     - 2  abs_d_sr_raref richness models (Gaussian)
#     - 2  prop_sr_full   richness models (Beta) — proportion of FULL spp detected
#     - 4  mean_p_full variants (TTE / LRR / SCOV / raref)
#     - 1  mean_p_full variant for prop_sr_full
#   guild_minor_habitat : carnivore size + ungulate habitat use
#   guild_minor_diet    : carnivore size + ungulate diet (browse/mix/graze)
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

N_WORKERS    <- 4   # models fitting simultaneously
CORES_PER_MODEL <- 4   # Stan chains per model
# N_WORKERS × CORES_PER_MODEL should not exceed physical cores (16)

plan(multisession, workers = N_WORKERS)

# Confirm plan
message("Parallel plan: ", class(plan())[1])
message("Workers: ", N_WORKERS, " × ", CORES_PER_MODEL, " cores = ",
        N_WORKERS * CORES_PER_MODEL, " cores total")


############################################################
# 1) DATA PREPARATION  (unchanged from models_abs_UPDATED.R)
############################################################

prep_protocol_summaries <- function(all_window_species,
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
      w_truth      = log1p(trap_days_full),
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
  
  # Community mean detectability from the longest (best-effort) window per
  # dataset × species — a dataset-level covariate independent of protocol.
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
      abs_d_p_tte       = median(abs_d_p_tte,       na.rm = TRUE),
      abs_d_log_rate    = median(abs_d_log_rate,     na.rm = TRUE),
      abs_d_spatial_cov = median(abs_d_spatial_cov,  na.rm = TRUE),
      l_trapdays        = median(l_trapdays,          na.rm = TRUE),
      l_nsites          = median(l_nsites,            na.rm = TRUE),
      s_latitude        = median(s_latitude,          na.rm = TRUE),
      s_trap_array      = median(s_trap_array,        na.rm = TRUE),
      w_truth           = median(w_truth,             na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(richness,        by = c("dataset", "protocol")) %>%
    left_join(community_p_full, by = "dataset") %>%
    mutate(l_n_species = log1p(n_species))
}


############################################################
# 2) MODEL REGISTRY  (14 models, unchanged from UPDATED)
############################################################

model_specs <- list(
  
  # ── TTE ────────────────────────────────────────────────
  tte_cov_base = list(
    resp = "abs_d_p_tte", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  tte_cov_spxprot = list(
    resp = "abs_d_p_tte", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 + protocol_bin || species)
  ),
  tte_cov_guild = list(
    resp = "abs_d_p_tte", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin * guild_major +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  tte_cov_richness_xprot = list(
    resp = "abs_d_p_tte", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin * l_n_species +
      l_trapdays + l_nsites + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  
  # ── LRR ────────────────────────────────────────────────
  lrr_cov_base = list(
    resp = "abs_d_log_rate", family = lognormal(),
    formula = abs_d_log_rate ~ 1 + protocol_bin +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  lrr_cov_spxprot = list(
    resp = "abs_d_log_rate", family = lognormal(),
    formula = abs_d_log_rate ~ 1 + protocol_bin +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 + protocol_bin || species)
  ),
  lrr_cov_guild = list(
    resp = "abs_d_log_rate", family = lognormal(),
    formula = abs_d_log_rate ~ 1 + protocol_bin * guild_major +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  lrr_cov_richness_xprot = list(
    resp = "abs_d_log_rate", family = lognormal(),
    formula = abs_d_log_rate ~ 1 + protocol_bin * l_n_species +
      l_trapdays + l_nsites + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  
  # ── SCOV ───────────────────────────────────────────────
  scov_cov_base = list(
    resp = "abs_d_spatial_cov", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  scov_cov_spxprot = list(
    resp = "abs_d_spatial_cov", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 + protocol_bin || species)
  ),
  scov_cov_guild = list(
    resp = "abs_d_spatial_cov", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin * guild_major +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  scov_cov_richness_xprot = list(
    resp = "abs_d_spatial_cov", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin * l_n_species +
      l_trapdays + l_nsites + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),

  # ── GUILD_MINOR HABITAT (carnivore size + ungulate habitat use) ─────────────
  tte_cov_guild_minor_hab = list(
    resp = "abs_d_p_tte", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin * guild_minor_habitat +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  lrr_cov_guild_minor_hab = list(
    resp = "abs_d_log_rate", family = lognormal(),
    formula = abs_d_log_rate ~ 1 + protocol_bin * guild_minor_habitat +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  scov_cov_guild_minor_hab = list(
    resp = "abs_d_spatial_cov", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin * guild_minor_habitat +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),

  # ── GUILD_MINOR DIET (carnivore size + ungulate diet browse/mix/graze) ───────
  tte_cov_guild_minor_diet = list(
    resp = "abs_d_p_tte", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin * guild_minor_diet +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  lrr_cov_guild_minor_diet = list(
    resp = "abs_d_log_rate", family = lognormal(),
    formula = abs_d_log_rate ~ 1 + protocol_bin * guild_minor_diet +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  scov_cov_guild_minor_diet = list(
    resp = "abs_d_spatial_cov", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin * guild_minor_diet +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),

  # ── RICHNESS: abs deviation of rarefied richness (Gaussian) ───────────────
  rich_base = list(
    resp = "abs_d_sr_raref", family = gaussian(),
    formula = abs_d_sr_raref ~ 1 + protocol_bin +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | dataset)
  ),
  rich_ds_slope = list(
    resp = "abs_d_sr_raref", family = gaussian(),
    formula = abs_d_sr_raref ~ 1 + protocol_bin +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 + protocol_bin || dataset),
    adapt_delta = 0.99   # 2 divergences on previous run
  ),

  # ── RICHNESS: MSE criterion (bias² + sampling variance) — Gaussian ─────────
  # Direct analogue of mse_p_tte / mse_log_rate / mse_spatial_cov.
  # sr_raref_var is the analytical incidence-rarefaction variance (Colwell 2012),
  # so mse_sr_raref = d_sr_raref² + sr_raref_var captures both bias and imprecision.
  rich_mse_base = list(
    resp = "mse_sr_raref", family = gaussian(),
    formula = mse_sr_raref ~ 1 + protocol_bin +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | dataset)
  ),
  rich_mse_ds_slope = list(
    resp = "mse_sr_raref", family = gaussian(),
    formula = mse_sr_raref ~ 1 + protocol_bin +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 + protocol_bin || dataset)
  ),

  # ── RICHNESS: proportion of FULL species detected (Beta) ──────────────────
  # prop_sr_full = |A ∩ FULL| / |FULL|  ∈ (0, 1)
  # Higher = more of the full-year list recovered; modelled directly (not as deviation).
  # Covariates match the rarefied-richness models but no effort-confound correction
  # needed because the denominator is fixed (FULL species count per dataset).
  rich_prop_base = list(
    resp = "prop_sr_full", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 | dataset)
  ),
  rich_prop_ds_slope = list(
    resp = "prop_sr_full", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin +
      l_trapdays + l_nsites + l_n_species + s_latitude + s_trap_array +
      (1 + protocol_bin || dataset)
  ),

  # ── COMMUNITY MEAN DETECTABILITY (s_mean_p_full) ─────────────────────────────
  # These models swap l_n_species for s_mean_p_full to test whether
  # the average detectability of the local assemblage (from the full reference
  # protocol) explains residual deviation better than species richness alone.
  tte_mean_p_full = list(
    resp = "abs_d_p_tte", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_mean_p_full + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  lrr_mean_p_full = list(
    resp = "abs_d_log_rate", family = lognormal(),
    formula = abs_d_log_rate ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_mean_p_full + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  scov_mean_p_full = list(
    resp = "abs_d_spatial_cov", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_mean_p_full + s_latitude + s_trap_array +
      (1 | ds_sp) + (1 | species)
  ),
  rich_mean_p_full = list(
    resp = "abs_d_sr_raref", family = gaussian(),
    formula = abs_d_sr_raref ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_mean_p_full + s_latitude + s_trap_array +
      (1 | dataset)
  ),
  rich_mse_mean_p_full = list(
    resp = "mse_sr_raref", family = gaussian(),
    formula = mse_sr_raref ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_mean_p_full + s_latitude + s_trap_array +
      (1 | dataset)
  ),
  rich_prop_mean_p_full = list(
    resp = "prop_sr_full", family = Beta(link = "logit"),
    formula = y ~ 1 + protocol_bin +
      l_trapdays + l_nsites + s_mean_p_full + s_latitude + s_trap_array +
      (1 | dataset)
  ),
)


############################################################
# 3) PARALLEL FITTER
############################################################

# Single-model fitting function (called inside each worker)
fit_one <- function(spec_name, df_sum, df_rich,
                    iter        = 4000,
                    warmup      = 2000,
                    chains      = 4,
                    cores       = CORES_PER_MODEL,
                    adapt_delta = 0.98,
                    max_treedepth = 12) {
  
  # Packages must be loaded inside each worker
  library(brms)
  library(dplyr)
  
  sp <- model_specs[[spec_name]]
  stopifnot(!is.null(sp))
  
  # Per-model adapt_delta overrides the function default if specified in spec
  adapt_delta    <- if (!is.null(sp$adapt_delta))    sp$adapt_delta    else adapt_delta
  max_treedepth  <- if (!is.null(sp$max_treedepth))  sp$max_treedepth  else max_treedepth
  
  # Richness models use a different data frame
  df <- if (startsWith(spec_name, "rich_")) df_rich else df_sum
  
  # For Beta-family richness models, drop rows where the response is NA or
  # outside (0,1) before applying the Smithson & Verkuilen transformation
  if (inherits(sp$family, "brmsfamily") && sp$family$family == "beta") {
    resp_vals <- df[[sp$resp]]
    df <- df[!is.na(resp_vals) & is.finite(resp_vals) &
               resp_vals > 0 & resp_vals < 1, , drop = FALSE]
    N    <- nrow(df)
    df$y <- (df[[sp$resp]] * (N - 1) + 0.5) / N
  }
  
  brm(
    formula   = sp$formula,
    data      = df,
    family    = sp$family,
    iter      = iter,
    warmup    = warmup,
    chains    = chains,
    cores     = cores,
    seed      = 42,
    control   = list(adapt_delta   = adapt_delta,
                     max_treedepth = max_treedepth),
    save_pars = save_pars(all = TRUE)
  )
}


############################################################
# 4) DATA PREPARATION
############################################################

df_sum <- prep_protocol_summaries(all_window_species)

traits <- read_csv("species_meta.csv", show_col_types = FALSE) %>%
  select(species, guild_major, family, guild_minor_habitat, guild_minor_diet)

df_sum <- df_sum %>%
  left_join(traits, by = "species") %>%
  mutate(
    guild_major         = factor(guild_major),
    family              = factor(family),
    guild_minor_habitat = factor(guild_minor_habitat),
    guild_minor_diet    = factor(guild_minor_diet)
  )

# Richness data (dataset-level, not species-level)
# latitude, trap_array, n_species not in all_window_richness — join from df_sum
dataset_covs <- df_sum %>%
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
    abs_d_sr_raref = abs(d_sr_raref),
    # mse_sr_raref = bias² + sampling variance (already computed in Full1.R)
    # Keep as-is; Gaussian family, non-negative by construction
    mse_sr_raref   = mse_sr_raref,
    # Richness-only datasets (no TTE data) get grand mean (= 0 on scaled)
    s_mean_p_full  = if_else(is.na(s_mean_p_full), 0, s_mean_p_full)
  ) %>%
  filter(
    is.finite(abs_d_sr_raref),
    # prop_sr_full must be in (0,1) for Beta; NA rows dropped
    # at fit time via the resp-specific filter in fit_one
    TRUE
  )

message("df_sum N = ", nrow(df_sum))
message("df_rich N = ", nrow(df_rich))


############################################################
# 5) FIT ALL MODELS IN PARALLEL BATCHES
############################################################
#
#  N_WORKERS = 4 workers, each running one model (4 Stan chains)
#  14 models → 4 batches of 4, 4, 4, 2  (last batch leaves 2 workers idle)
#
############################################################

all_spec_names <- names(model_specs)

message("\n─── Starting parallel fit: ", length(all_spec_names),
        " models across ", N_WORKERS, " workers ───")
message("Started at: ", Sys.time())

t_start <- proc.time()

fits_list <- future_map(
  all_spec_names,
  \(nm) fit_one(nm, df_sum, df_rich),
  .options = furrr_options(seed = TRUE),   # safe RNG across workers
  .progress = TRUE                          # progress bar in console
)

names(fits_list) <- all_spec_names

t_elapsed <- proc.time() - t_start
message("\n─── All models fitted ───")
message("Total elapsed time: ", round(t_elapsed["elapsed"] / 60, 1), " minutes")
message("Finished at: ", Sys.time())


############################################################
# 6) UNPACK FITS INTO NAMED OBJECTS
############################################################

fit_tte_base            <- fits_list[["tte_cov_base"]]
fit_tte_spx             <- fits_list[["tte_cov_spxprot"]]
fit_tte_guild           <- fits_list[["tte_cov_guild"]]
fit_tte_richness_xprot  <- fits_list[["tte_cov_richness_xprot"]]
fit_tte_guild_minor_hab  <- fits_list[["tte_cov_guild_minor_hab"]]
fit_tte_guild_minor_diet <- fits_list[["tte_cov_guild_minor_diet"]]

fit_lrr_base            <- fits_list[["lrr_cov_base"]]
fit_lrr_spx             <- fits_list[["lrr_cov_spxprot"]]
fit_lrr_guild           <- fits_list[["lrr_cov_guild"]]
fit_lrr_richness_xprot  <- fits_list[["lrr_cov_richness_xprot"]]
fit_lrr_guild_minor_hab  <- fits_list[["lrr_cov_guild_minor_hab"]]
fit_lrr_guild_minor_diet <- fits_list[["lrr_cov_guild_minor_diet"]]

fit_scov_base           <- fits_list[["scov_cov_base"]]
fit_scov_spx            <- fits_list[["scov_cov_spxprot"]]
fit_scov_guild          <- fits_list[["scov_cov_guild"]]
fit_scov_richness_xprot <- fits_list[["scov_cov_richness_xprot"]]
fit_scov_guild_minor_hab  <- fits_list[["scov_cov_guild_minor_hab"]]
fit_scov_guild_minor_diet <- fits_list[["scov_cov_guild_minor_diet"]]

fit_rich_base            <- fits_list[["rich_base"]]
fit_rich_dsxprot         <- fits_list[["rich_ds_slope"]]

fit_rich_mse_base        <- fits_list[["rich_mse_base"]]
fit_rich_mse_dsxprot     <- fits_list[["rich_mse_ds_slope"]]

fit_rich_prop_base       <- fits_list[["rich_prop_base"]]
fit_rich_prop_dsxprot    <- fits_list[["rich_prop_ds_slope"]]

fit_rich_jaccard_base    <- fits_list[["rich_jaccard_base"]]
fit_rich_jaccard_dsxprot <- fits_list[["rich_jaccard_ds_slope"]]

fit_tte_mean_p_full          <- fits_list[["tte_mean_p_full"]]
fit_lrr_mean_p_full          <- fits_list[["lrr_mean_p_full"]]
fit_scov_mean_p_full         <- fits_list[["scov_mean_p_full"]]
fit_rich_mean_p_full         <- fits_list[["rich_mean_p_full"]]
fit_rich_mse_mean_p_full     <- fits_list[["rich_mse_mean_p_full"]]
fit_rich_prop_mean_p_full    <- fits_list[["rich_prop_mean_p_full"]]
fit_rich_jaccard_mean_p_full <- fits_list[["rich_jaccard_mean_p_full"]]


############################################################
# 7) SAVE ALL FITS
############################################################

saveRDS(fits_list, "fits_abs_all.rds")
message("Fits saved to fits_abs_all.rds")


############################################################
# 8) DIAGNOSTICS & POSTERIOR PREDICTIVE CHECKS
############################################################
fits_list <- readRDS("fits_abs_all.rds")

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
  print(pp_check(fit, ndraws = 100)                              + ggtitle(paste(label, "— density")))
  print(pp_check(fit, ndraws = 100, type = "stat",   stat = "mean") + ggtitle(paste(label, "— mean")))
  print(pp_check(fit, ndraws = 100, type = "stat_2d", stat = c("mean", "sd")) + ggtitle(paste(label, "— mean vs SD")))
}

pp_check_full(fits_list[["tte_cov_base"]],  "TTE base")
pp_check_full(fits_list[["lrr_cov_base"]],  "LRR base")
pp_check_full(fits_list[["scov_cov_base"]], "SCOV base")


############################################################
# 9) LOO-CV  (one worker per model — no Stan chains needed)
############################################################


message("\n─── LOO-CV ───")

# Scale up: 14 workers × 1 core each (14 physical cores used)
plan(multisession, workers = length(fits_list))
message("LOO workers: ", length(fits_list))

t_loo <- proc.time()

loo_list <- future_map(
  fits_list,
  function(fit) {
    library(brms)                        # must load inside each worker
    loo(fit, moment_match = TRUE)
  },
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)
names(loo_list) <- all_spec_names

message("LOO finished in ", round((proc.time() - t_loo)["elapsed"] / 60, 1), " minutes")

saveRDS(loo_list, "loo_abs_all.rds")


############################################################
# 10) RESET PLAN WHEN DONE
############################################################

plan(sequential)
message("Parallel plan reset to sequential.")

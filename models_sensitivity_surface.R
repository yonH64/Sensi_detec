# models_sensitivity_surface.R
# ============================================================================
# Sensitivity surface models: how does detection metric deviation from the
# 12-month benchmark vary as a function of window timing (day-of-year),
# window duration, species identity, and environmental context?
#
# Approach:  GAMMs via mgcv::bam() with cyclic splines for day-of-year,
#            tensor product surfaces for timing × duration interaction,
#            and species-specific or shared surfaces depending on model variant.
#
# Input:     sensitivity_species_data.rds   (from prep_sensitivity_data.R)
#            sensitivity_richness_data.rds   (from prep_sensitivity_data.R)
#
# Output:    sensitivity_gam_models.rds      — all fitted models in a named list
#            sensitivity_models_env.RData   — cleaned environment for downstream scripts
#
# Timing:    ~2–3 minutes total on 4 threads (86k species rows, 20k richness rows)
# ============================================================================

library(mgcv)
library(dplyr)
library(tibble)

N_THREADS <- 4
cat("=== Sensitivity surface models ===\n")
cat("Threads:", N_THREADS, "\n\n")

# ── Load data ────────────────────────────────────────────────────────────────

sens_species  <- readRDS("sensitivity_species_data.rds")
sens_richness <- readRDS("sensitivity_richness_data.rds")

cat("Species-level:", nrow(sens_species), "rows,",
    n_distinct(sens_species$species), "species,",
    n_distinct(sens_species$dataset), "datasets\n")
cat("Richness-level:", nrow(sens_richness), "rows,",
    n_distinct(sens_richness$dataset), "datasets\n\n")

# Cyclic knots — day_start runs 1–365 with 7-day steps
cc_knots <- list(day_start = c(0, 365))

# ── Helper: fit and time a bam() model ───────────────────────────────────────

fit_bam <- function(formula, data, family, knots = cc_knots, label = "") {
  cat(sprintf("  Fitting %-45s ... ", label))
  t0 <- proc.time()
  fit <- bam(
    formula  = formula,
    data     = data,
    family   = family,
    method   = "fREML",
    discrete = TRUE,
    nthreads = N_THREADS,
    knots    = knots
  )
  elapsed <- round((proc.time() - t0)[3], 1)
  cat(sprintf("done (%ss, dev.expl = %.1f%%)\n",
              elapsed, 100 * summary(fit)$dev.expl))
  fit
}

# ============================================================================
# PART 1: MODEL COMPARISON — lambda only
# ============================================================================
# Fits M1–M6 for abs_d_lambda to justify the species-specific surface (M6).
# All models include: s_bio4, l_trapdays, l_nsites, s_latitude, s_trap_array,
# and dataset×species random effects s(ds_sp_f, bs = "re").

cat("── Part 1: Model comparison (abs_d_lambda) ──────────────────────\n")

fits_comparison <- list()

# M1: Shared surface + covariates
fits_comparison$M1_base <- fit_bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(12, 8)) +
    s_bio4 + l_trapdays + l_nsites + s_latitude + s_trap_array +
    s(species_f, bs = "re") +
    s(ds_sp_f, bs = "re"),
  data = sens_species, family = Gamma(link = "log"),
  label = "M1: shared surface"
)

# M2: Shared surface + guild-varying seasonality
fits_comparison$M2_guild_season <- fit_bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(12, 8)) +
    s(day_start, bs = "cc", by = guild_major, k = 12) +
    guild_major +
    s_bio4 + l_trapdays + l_nsites + s_latitude + s_trap_array +
    s(species_f, bs = "re") +
    s(ds_sp_f, bs = "re"),
  data = sens_species, family = Gamma(link = "log"),
  label = "M2: + guild-varying seasonality"
)

# M3: Full guild × surface
fits_comparison$M3_guild_surface <- fit_bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(12, 8),
       by = guild_major) +
    guild_major +
    s_bio4 + l_trapdays + l_nsites + s_latitude + s_trap_array +
    s(species_f, bs = "re") +
    s(ds_sp_f, bs = "re"),
  data = sens_species, family = Gamma(link = "log"),
  label = "M3: guild-specific surfaces"
)

# M4: Shared surface + bio4 × duration interaction
fits_comparison$M4_env_x_dur <- fit_bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(12, 8)) +
    s(window_len, by = s_bio4, k = 8) +
    s_bio4 + l_trapdays + l_nsites + s_latitude + s_trap_array +
    s(species_f, bs = "re") +
    s(ds_sp_f, bs = "re"),
  data = sens_species, family = Gamma(link = "log"),
  label = "M4: + bio4 x duration"
)

# M5: Guild surface + bio4 × duration
fits_comparison$M5_guild_env <- fit_bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(12, 8),
       by = guild_major) +
    guild_major +
    s(window_len, by = s_bio4, k = 8) +
    s_bio4 + l_trapdays + l_nsites + s_latitude + s_trap_array +
    s(species_f, bs = "re") +
    s(ds_sp_f, bs = "re"),
  data = sens_species, family = Gamma(link = "log"),
  label = "M5: guild surface + bio4 x dur"
)

# M_hab: Minor guild (habitat) surfaces
fits_comparison$M_hab <- fit_bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6),
       by = guild_minor_habitat) +
    guild_minor_habitat +
    s_bio4 + l_trapdays + l_nsites + s_latitude + s_trap_array +
    s(species_f, bs = "re") +
    s(ds_sp_f, bs = "re"),
  data = sens_species, family = Gamma(link = "log"),
  label = "M_hab: minor habitat surfaces"
)

# M_diet: Minor guild (diet) surfaces
fits_comparison$M_diet <- fit_bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6),
       by = guild_minor_diet) +
    guild_minor_diet +
    s_bio4 + l_trapdays + l_nsites + s_latitude + s_trap_array +
    s(species_f, bs = "re") +
    s(ds_sp_f, bs = "re"),
  data = sens_species, family = Gamma(link = "log"),
  label = "M_diet: minor diet surfaces"
)

# M6: Species-specific surfaces (the best model from pilot)
fits_comparison$M6_species <- fit_bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6),
       by = species_f) +
    species_f +
    s_bio4 + l_trapdays + l_nsites + s_latitude + s_trap_array +
    s(ds_sp_f, bs = "re"),
  data = sens_species, family = Gamma(link = "log"),
  label = "M6: species-specific surfaces"
)

# Print comparison table
cat("\n── Model comparison summary ────────────────────────────────────\n")
comparison_tbl <- tibble(
  model = names(fits_comparison),
  AIC   = sapply(fits_comparison, AIC),
  dev_expl = sapply(fits_comparison, \(m) summary(m)$dev.expl)
) |>
  arrange(AIC) |>
  mutate(delta_AIC = AIC - min(AIC),
         AIC = round(AIC),
         dev_expl = round(dev_expl, 3),
         delta_AIC = round(delta_AIC))
print(comparison_tbl, n = 10)
cat("\n")


# ============================================================================
# PART 2: SPECIES-LEVEL DETECTION MODELS (M6 structure, all metrics)
# ============================================================================
# Fits the winning M6 structure for each detection metric.

cat("── Part 2: Species-level detection models (M6) ────────────────\n")

# Shared formula components (RHS excluding response-specific terms)
# Species surfaces + covariates + random effects
m6_rhs <- ~ te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6),
               by = species_f) +
  species_f +
  s_bio4 + l_trapdays + l_nsites + s_latitude + s_trap_array +
  s(ds_sp_f, bs = "re")

fits_detection <- list()

# 2a. |d_lambda| — TTE daily detection rate deviation (absolute)
fits_detection$abs_d_lambda <- fits_comparison$M6_species  # already fitted above

# 2b. |d_rate| — encounter rate deviation (absolute)
fits_detection$abs_d_rate <- fit_bam(
  update(m6_rhs, abs_d_rate ~ .),
  data = sens_species, family = Gamma(link = "log"),
  label = "abs_d_rate (Gamma)"
)

# 2c. |d_matched_rate| — deconfounded spatial-coverage-based rate (absolute)
fits_detection$abs_d_matched_rate <- fit_bam(
  update(m6_rhs, abs_d_matched_rate ~ .),
  data = sens_species, family = Gamma(link = "log"),
  label = "abs_d_matched_rate (Gamma)"
)

# 2d. d_rate (signed) — bidirectional encounter rate deviation
fits_detection$d_rate_signed <- fit_bam(
  update(m6_rhs, d_rate ~ .),
  data = sens_species, family = gaussian(),
  label = "d_rate signed (Gaussian)"
)

cat("\n")


# ============================================================================
# PART 3: COMMUNITY-LEVEL MODELS (richness & similarity)
# ============================================================================
# Dataset-level models with shared (not species-specific) surfaces.

cat("── Part 3: Community-level models ─────────────────────────────\n")

fits_community <- list()

# 3a. Rarefied richness deviation — Gaussian (can be negative)
fits_community$d_sr_raref <- fit_bam(
  d_sr_raref ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(12, 8)) +
    s_bio4 + l_trapdays + l_nsites + s_latitude +
    s(dataset_f, bs = "re"),
  data = sens_richness, family = gaussian(),
  label = "d_sr_raref (Gaussian)"
)

# 3b. Proportion of full-year species recovered — Beta
sens_rich_beta <- sens_richness |>
  filter(!is.na(prop_sr_full)) |>
  mutate(prop_sr_adj = pmin(pmax(prop_sr_full, 0.001), 0.999))

fits_community$prop_sr_full <- fit_bam(
  prop_sr_adj ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(12, 8)) +
    s_bio4 + l_trapdays + l_nsites + s_latitude +
    s(dataset_f, bs = "re"),
  data = sens_rich_beta, family = betar(link = "logit"),
  label = "prop_sr_full (Beta)"
)

# 3c. Rank correlation of detection rates — Beta on transformed scale
# Filter to windows with >= 5 shared species: Spearman's rho with 3–4 species
# has very few possible values and produces 33% perfect correlations (rho = 1),
# creating boundary mass that violates beta regression assumptions. Filtering
# to >= 5 shared species removes most artefactual boundary inflation and more
# than doubles deviance explained (17% → 37%).
sens_rich_rho <- sens_richness |>
  filter(!is.na(rho_lambda), n_shared_spp >= 5) |>
  mutate(rho_adj = pmin(pmax((rho_lambda + 1) / 2, 0.001), 0.999))

fits_community$rho_lambda <- fit_bam(
  rho_adj ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(12, 8)) +
    s_bio4 + l_trapdays + l_nsites + s_latitude +
    s(dataset_f, bs = "re"),
  data = sens_rich_rho, family = betar(link = "logit"),
  label = "rho_lambda (Beta)"
)

cat("\n")


# ============================================================================
# PART 4: DIAGNOSTICS
# ============================================================================

cat("── Part 4: Diagnostics ────────────────────────────────────────\n\n")

# Summary table of all final models
all_final <- c(fits_detection, fits_community)

diag_tbl <- tibble(
  model       = names(all_final),
  family      = sapply(all_final, \(m) m$family$family),
  link        = sapply(all_final, \(m) m$family$link),
  n_obs       = sapply(all_final, \(m) m$df.null + 1),
  edf_total   = sapply(all_final, \(m) round(sum(m$edf), 1)),
  dev_expl    = sapply(all_final, \(m) round(summary(m)$dev.expl, 3)),
  r_sq        = sapply(all_final, \(m) round(summary(m)$r.sq, 3)),
  AIC         = sapply(all_final, \(m) round(AIC(m)))
)
print(diag_tbl, n = 20)

# Concurvity check for the main model (abs_d_lambda)
cat("\n── Concurvity (abs_d_lambda, worst case) ──\n")
cc <- concurvity(fits_detection$abs_d_lambda, full = TRUE)
print(round(cc, 3))

cat("\n")


# ============================================================================
# PART 5: SAVE
# ============================================================================

all_models <- list(
  comparison = fits_comparison,
  detection  = fits_detection,
  community  = fits_community
)

saveRDS(all_models, "sensitivity_gam_models.rds")
cat("Saved all models to sensitivity_gam_models.rds\n")
cat(sprintf("Total models: %d comparison + %d detection + %d community = %d\n",
            length(fits_comparison), length(fits_detection),
            length(fits_community),
            length(fits_comparison) + length(fits_detection) + length(fits_community)))

# Clean up intermediate objects and save working environment
rm(fits_comparison, fits_detection, fits_community, all_final,
   sens_rich_beta, sens_rich_rho,
   m6_rhs, comparison_tbl, diag_tbl, cc_knots, N_THREADS, fit_bam)
save.image("sensitivity_models_env.RData")
cat("Saved environment to sensitivity_models_env.RData\nDone.\n")

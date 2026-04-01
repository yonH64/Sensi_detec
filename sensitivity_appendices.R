# ═══════════════════════════════════════════════════════════════════════════════
# sensitivity_appendices.R
# ═══════════════════════════════════════════════════════════════════════════════
#
# Unified script for all robustness and sensitivity analyses documented in the
# appendices. Each section corresponds to one appendix markdown file.
#
# Sections:
#   1. Species-inclusion threshold robustness   (appendix_threshold_robustness.md)
#   2. Joint threshold sensitivity              (appendix_joint_threshold_sensitivity.md)
#   3. Window-start resolution sensitivity      (appendix_resolution_sensitivity.md)
#   4. Anchor detection parameter sensitivity   (appendix_anchor_sensitivity.md)
#   5. Independence threshold sensitivity       (appendix_independence_threshold_sensitivity.md)
#   6. Model specification diagnostics          (appendix_model_diagnostics.md)
#   7. Rho filter sensitivity                   (appendix_rho_filter_sensitivity.md)
#
# Inputs:
#   - sensitivity_models_env.RData    (fitted models + prepared data)
#   - all_window_species.rds          (raw window metrics)
#   - all_window_richness.rds         (raw richness metrics)
#   - species_meta.csv
#   - dataset_metadata.csv
#   - covariate_scaling_constants.csv
#
# For sections 3–5 (resolution, anchor, independence), the full data pipeline
# must be re-run with different parameters. Pre-computed data is loaded from:
#   - all_window_species_3day.rds     (3-day step, from Full1.R)
#   - all_window_species_15min.rds    (15-min gap, from Full1.R)
#   - all_window_species_60min.rds    (60-min gap, from Full1.R)
#   - all_window_species_strict.rds   (strict anchors, from Full1.R)
# If these files don't exist, those sections are skipped with a message.
#
# Outputs:
#   - threshold_sensitivity_summary.csv
#   - threshold_sensitivity_sites_summary.csv
#   - threshold_sensitivity_joint_summary.csv
#   - resolution_sensitivity_summary.csv
#   - anchor_sensitivity_summary.csv
#   - independence_threshold_summary.csv
#   - model_diagnostics_phase3.csv
#   - rho_filter_sensitivity_summary.csv
#   - figures/threshold_sensitivity_comparison.pdf
#   - figures/threshold_sensitivity_sites_comparison.pdf
#   - figures/threshold_sensitivity_joint.pdf
#   - figures/resolution_sensitivity.pdf
#   - figures/anchor_sensitivity.pdf
#   - figures/independence_threshold_sensitivity.pdf
#   - figures/rho_filter_sensitivity.pdf
#
# Runtime: ~60–90 minutes (model refits dominate; ~2 min per bam() call).
#          Sections 3–5 require pre-computed pipeline data.
# ═══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(mgcv)
library(patchwork)

N_THREADS <- 4
cc_knots  <- list(day_start = c(0, 365))

dir.create("figures", showWarnings = FALSE)

# ── Load core data ───────────────────────────────────────────────────────────

load("sensitivity_models_env.RData")

species_meta <- read.csv("species_meta.csv")
dataset_meta <- read.csv("dataset_metadata.csv")
scaling      <- read.csv("covariate_scaling_constants.csv")

lat_center       <- scaling$center[scaling$model_name == "s_latitude"]
lat_sd           <- scaling$scale[scaling$model_name == "s_latitude"]
trap_array_center <- scaling$center[scaling$model_name == "s_trap_array"]
trap_array_sd    <- scaling$scale[scaling$model_name == "s_trap_array"]

mod_lambda_baseline <- all_models$detection$abs_d_lambda
mod_rho_baseline    <- all_models$community$rho_lambda

cat("Baseline data: ", nrow(sens_species), "rows,",
    n_distinct(sens_species$species_f), "species\n")


# ═══════════════════════════════════════════════════════════════════════════════
# SHARED UTILITIES
# ═══════════════════════════════════════════════════════════════════════════════

# ── Data preparation (from raw window metrics to model-ready data) ───────────

prep_species_data <- function(ws_data, min_events = 20, min_sites_pos = 5,
                              min_occasions_pos = 5, min_rows_species = 50) {
  # Apply thresholds
  ws_data <- ws_data |>
    filter(
      n_events_total  >= min_events,
      n_sites_pos     >= min_sites_pos,
      n_occasions_pos >= min_occasions_pos
    )

  out <- ws_data |>
    filter(!window_id %in% c("SNAP_EU_CORE", "SNAP_EU_BUFFER",
                              "EOW_EARLY", "EOW_LATE")) |>
    mutate(
      day_start          = as.integer(str_extract(window_id, "(?<=d)\\d+")),
      abs_d_lambda       = pmax(abs(d_lambda), 1e-10),
      abs_d_rate         = pmax(abs(d_rate), 1e-10),
      abs_d_matched_rate = pmax(abs(d_matched_rate), 1e-10),
      l_trapdays         = log1p(trap_days_window),
      l_nsites           = log1p(n_sites),
      ds_sp              = interaction(dataset, species, drop = TRUE)
    ) |>
    left_join(
      species_meta |> select(species, guild_major, guild_minor_habitat, guild_minor_diet),
      by = "species"
    ) |>
    filter(!is.na(guild_major))

  # Drop species with too few rows
  spp_counts <- out |> count(species)
  keep_spp   <- spp_counts |> filter(n >= min_rows_species) |> pull(species)
  out <- out |> filter(species %in% keep_spp)

  # Standardise covariates using CANONICAL constants
  out <- out |>
    mutate(
      s_latitude   = (latitude - lat_center) / lat_sd,
      s_trap_array = (log1p(trap_array) - trap_array_center) / trap_array_sd,
      species_f    = factor(species),
      dataset_f    = factor(dataset),
      ds_sp_f      = factor(ds_sp)
    )

  out
}


# ── Model fitting ────────────────────────────────────────────────────────────

fit_m6 <- function(data, response = "abs_d_lambda", label = "") {
  fml <- as.formula(paste0(
    response, " ~ ",
    "te(day_start, window_len, bs = c('cc', 'tp'), k = c(8, 6), by = species_f) + ",
    "species_f + l_trapdays + l_nsites + s_latitude + s_trap_array + ",
    "s(ds_sp_f, bs = 're')"
  ))

  cat(sprintf("  Fitting M6 %-15s %-30s ... ", response, label))
  t0 <- proc.time()

  fit <- bam(fml, data = data, family = Gamma(link = "log"),
             method = "fREML", discrete = TRUE, nthreads = N_THREADS,
             knots = cc_knots)

  elapsed <- round((proc.time() - t0)[3], 1)
  cat(sprintf("done (%ss, dev.expl = %.1f%%)\n", elapsed,
              100 * summary(fit)$dev.expl))
  fit
}


# ── Surface comparison ───────────────────────────────────────────────────────

compare_surfaces <- function(fit_a, fit_b, data_a, data_b,
                             shared_species = NULL) {
  # Build prediction grid on shared species
  if (is.null(shared_species)) {
    shared_species <- intersect(
      levels(data_a$species_f),
      levels(data_b$species_f)
    )
  }

  # Use a plain factor with only shared species as levels to avoid

  # crossing() expanding unused factor levels (which causes NA predictions).
  shared_f <- factor(shared_species)

  grid_core <- tidyr::crossing(
    day_start  = seq(0, 364, by = 7),
    window_len = seq(15, 183, by = 7),
    species_f  = shared_f
  )

  grid_a <- grid_core |>
    mutate(
      species_f    = factor(as.character(species_f), levels = levels(data_a$species_f)),
      l_trapdays   = median(data_a$l_trapdays),
      l_nsites     = median(data_a$l_nsites),
      s_latitude   = 0,
      s_trap_array = 0,
      ds_sp_f      = data_a$ds_sp_f[1]
    )

  grid_b <- grid_core |>
    mutate(
      species_f    = factor(as.character(species_f), levels = levels(data_b$species_f)),
      l_trapdays   = median(data_b$l_trapdays),
      l_nsites     = median(data_b$l_nsites),
      s_latitude   = 0,
      s_trap_array = 0,
      ds_sp_f      = data_b$ds_sp_f[1]
    )

  pred_a <- predict(fit_a, newdata = grid_a, type = "response", exclude = "s(ds_sp_f)")
  pred_b <- predict(fit_b, newdata = grid_b, type = "response", exclude = "s(ds_sp_f)")

  list(
    r_surface = cor(as.numeric(pred_a), as.numeric(pred_b)),
    pred_a    = as.numeric(pred_a),
    pred_b    = as.numeric(pred_b),
    grid      = grid_a
  )
}


# ── Summary row helper ───────────────────────────────────────────────────────

model_summary_row <- function(fit, label, data) {
  s <- summary(fit)
  tibble(
    config    = label,
    n_species = n_distinct(data$species_f),
    n_rows    = nrow(data),
    dev_expl  = round(100 * s$dev.expl, 2),
    r_sq_adj  = round(s$r.sq, 4)
  )
}


# ── Duration curve extraction ────────────────────────────────────────────────

extract_duration_curve <- function(fit, data, label,
                                   shared_species = NULL) {
  species_pool <- if (!is.null(shared_species)) shared_species else levels(data$species_f)

  grid <- tidyr::crossing(
    day_start  = seq(0, 364, by = 14),
    window_len = seq(15, 183, by = 7),
    species_f  = factor(species_pool, levels = levels(data$species_f))
  ) |>
    mutate(
      l_trapdays  = median(data$l_trapdays),
      l_nsites    = median(data$l_nsites),
      s_latitude  = 0,
      s_trap_array = 0,
      ds_sp_f     = data$ds_sp_f[1]
    )

  grid$pred <- predict(fit, newdata = grid, type = "response", exclude = "s(ds_sp_f)")

  grid |>
    group_by(window_len) |>
    summarise(mean_pred = mean(pred), .groups = "drop") |>
    mutate(config = label)
}


# ── Seasonal profile extraction (fixed duration) ────────────────────────────

extract_seasonal_profile <- function(fit, data, label, duration = 57,
                                     shared_species = NULL) {
  species_pool <- if (!is.null(shared_species)) shared_species else levels(data$species_f)

  grid <- tidyr::crossing(
    day_start  = seq(0, 364, by = 7),
    window_len = duration,
    species_f  = factor(species_pool, levels = levels(data$species_f))
  ) |>
    mutate(
      l_trapdays  = median(data$l_trapdays),
      l_nsites    = median(data$l_nsites),
      s_latitude  = 0,
      s_trap_array = 0,
      ds_sp_f     = data$ds_sp_f[1]
    )

  grid$pred <- predict(fit, newdata = grid, type = "response", exclude = "s(ds_sp_f)")

  grid |>
    group_by(day_start) |>
    summarise(mean_pred = mean(pred), .groups = "drop") |>
    mutate(config = label)
}


# ── Standard 2-panel comparison figure ───────────────────────────────────────

plot_sensitivity_comparison <- function(dur_curves, seasonal_profiles,
                                        title = "Sensitivity comparison") {
  p1 <- ggplot(dur_curves, aes(window_len, mean_pred, colour = config)) +
    geom_line(linewidth = 0.8) +
    labs(x = "Window duration (days)", y = "Mean predicted |d_lambda|",
         title = "Duration curve", colour = NULL) +
    theme_minimal(base_size = 10)

  p2 <- ggplot(seasonal_profiles, aes(day_start, mean_pred, colour = config)) +
    geom_line(linewidth = 0.8) +
    scale_x_continuous(breaks = c(0, 91, 182, 274),
                       labels = c("Jan", "Apr", "Jul", "Oct")) +
    labs(x = "Window start", y = "Mean predicted |d_lambda|",
         title = "Seasonal profile (57-day windows)", colour = NULL) +
    theme_minimal(base_size = 10)

  p1 + p2 + plot_annotation(title = title)
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Species-Inclusion Threshold Robustness (Individual)
#             → appendix_threshold_robustness.md
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\n════ §1: Individual species-inclusion thresholds ════\n\n")

# Load raw window data for re-filtering
all_ws <- readRDS("all_window_species.rds")

# Baseline is already in sens_species — just refit to confirm
baseline_data <- sens_species
baseline_fit  <- mod_lambda_baseline

# ── 1A: Minimum events threshold ─────────────────────────────────────────────

cat("--- 1A: min_events (10/20/30) ---\n")

data_events_10 <- prep_species_data(all_ws, min_events = 10)
data_events_30 <- prep_species_data(all_ws, min_events = 30)

fit_events_10 <- fit_m6(data_events_10, label = "(events ≥ 10)")
fit_events_30 <- fit_m6(data_events_30, label = "(events ≥ 30)")

shared_spp_events <- Reduce(intersect, list(
  levels(droplevels(baseline_data$species_f)),
  levels(droplevels(data_events_10$species_f)),
  levels(droplevels(data_events_30$species_f))
))

comp_e10 <- compare_surfaces(baseline_fit, fit_events_10, baseline_data, data_events_10,
                              shared_spp_events)
comp_e30 <- compare_surfaces(baseline_fit, fit_events_30, baseline_data, data_events_30,
                              shared_spp_events)

events_summary <- bind_rows(
  model_summary_row(fit_events_10, "min_events=10", data_events_10) |>
    mutate(surface_r = comp_e10$r_surface),
  model_summary_row(baseline_fit, "min_events=20 (baseline)", baseline_data) |>
    mutate(surface_r = 1.0),
  model_summary_row(fit_events_30, "min_events=30", data_events_30) |>
    mutate(surface_r = comp_e30$r_surface)
)

# Duration curves and seasonal profiles
dur_events <- bind_rows(
  extract_duration_curve(fit_events_10, data_events_10, "events ≥ 10", shared_spp_events),
  extract_duration_curve(baseline_fit, baseline_data, "events ≥ 20 (baseline)", shared_spp_events),
  extract_duration_curve(fit_events_30, data_events_30, "events ≥ 30", shared_spp_events)
)

seas_events <- bind_rows(
  extract_seasonal_profile(fit_events_10, data_events_10, "events ≥ 10", shared_species = shared_spp_events),
  extract_seasonal_profile(baseline_fit, baseline_data, "events ≥ 20 (baseline)", shared_species = shared_spp_events),
  extract_seasonal_profile(fit_events_30, data_events_30, "events ≥ 30", shared_species = shared_spp_events)
)

rm(fit_events_10, fit_events_30, data_events_10, data_events_30)
gc(verbose = FALSE)


# ── 1B: Minimum positive sites threshold ─────────────────────────────────────

cat("\n--- 1B: min_sites_pos (5/10/15) ---\n")

data_sites_10 <- prep_species_data(all_ws, min_sites_pos = 10)
data_sites_15 <- prep_species_data(all_ws, min_sites_pos = 15)

fit_sites_10 <- fit_m6(data_sites_10, label = "(sites ≥ 10)")
fit_sites_15 <- fit_m6(data_sites_15, label = "(sites ≥ 15)")

shared_spp_sites <- Reduce(intersect, list(
  levels(droplevels(baseline_data$species_f)),
  levels(droplevels(data_sites_10$species_f)),
  levels(droplevels(data_sites_15$species_f))
))

comp_s10 <- compare_surfaces(baseline_fit, fit_sites_10, baseline_data, data_sites_10,
                              shared_spp_sites)
comp_s15 <- compare_surfaces(baseline_fit, fit_sites_15, baseline_data, data_sites_15,
                              shared_spp_sites)

sites_summary <- bind_rows(
  model_summary_row(baseline_fit, "min_sites=5 (baseline)", baseline_data) |>
    mutate(surface_r = 1.0),
  model_summary_row(fit_sites_10, "min_sites=10", data_sites_10) |>
    mutate(surface_r = comp_s10$r_surface),
  model_summary_row(fit_sites_15, "min_sites=15", data_sites_15) |>
    mutate(surface_r = comp_s15$r_surface)
)

dur_sites <- bind_rows(
  extract_duration_curve(baseline_fit, baseline_data, "sites ≥ 5 (baseline)", shared_spp_sites),
  extract_duration_curve(fit_sites_10, data_sites_10, "sites ≥ 10", shared_spp_sites),
  extract_duration_curve(fit_sites_15, data_sites_15, "sites ≥ 15", shared_spp_sites)
)

seas_sites <- bind_rows(
  extract_seasonal_profile(baseline_fit, baseline_data, "sites ≥ 5 (baseline)", shared_species = shared_spp_sites),
  extract_seasonal_profile(fit_sites_10, data_sites_10, "sites ≥ 10", shared_species = shared_spp_sites),
  extract_seasonal_profile(fit_sites_15, data_sites_15, "sites ≥ 15", shared_species = shared_spp_sites)
)

rm(fit_sites_10, fit_sites_15, data_sites_10, data_sites_15)
gc(verbose = FALSE)

# Save summaries
write_csv(bind_rows(events_summary, sites_summary), "threshold_sensitivity_summary.csv")

# Save figures
ggsave("figures/threshold_sensitivity_comparison.pdf",
       plot_sensitivity_comparison(dur_events, seas_events,
         "Threshold sensitivity: min_events"),
       width = 11, height = 5)

ggsave("figures/threshold_sensitivity_sites_comparison.pdf",
       plot_sensitivity_comparison(dur_sites, seas_sites,
         "Threshold sensitivity: min_sites_pos"),
       width = 11, height = 5)

cat("  Saved: threshold_sensitivity_summary.csv\n")
cat("  Saved: figures/threshold_sensitivity_comparison.pdf\n")
cat("  Saved: figures/threshold_sensitivity_sites_comparison.pdf\n")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: Joint Threshold Sensitivity
#             → appendix_joint_threshold_sensitivity.md
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\n════ §2: Joint species-inclusion thresholds ════\n\n")

data_joint_lenient <- prep_species_data(all_ws, min_events = 10,
                                         min_sites_pos = 3, min_occasions_pos = 3)
data_joint_strict  <- prep_species_data(all_ws, min_events = 30,
                                         min_sites_pos = 10, min_occasions_pos = 5)

fit_joint_lenient <- fit_m6(data_joint_lenient, label = "(10/3/3)")
fit_joint_strict  <- fit_m6(data_joint_strict,  label = "(30/10/5)")

shared_spp_joint <- Reduce(intersect, list(
  levels(droplevels(baseline_data$species_f)),
  levels(droplevels(data_joint_lenient$species_f)),
  levels(droplevels(data_joint_strict$species_f))
))

comp_jl <- compare_surfaces(baseline_fit, fit_joint_lenient, baseline_data,
                             data_joint_lenient, shared_spp_joint)
comp_js <- compare_surfaces(baseline_fit, fit_joint_strict, baseline_data,
                             data_joint_strict, shared_spp_joint)

joint_summary <- bind_rows(
  model_summary_row(fit_joint_lenient, "Lenient (10/3/3)", data_joint_lenient) |>
    mutate(surface_r = comp_jl$r_surface),
  model_summary_row(baseline_fit, "Baseline (20/5/5)", baseline_data) |>
    mutate(surface_r = 1.0),
  model_summary_row(fit_joint_strict, "Strict (30/10/5)", data_joint_strict) |>
    mutate(surface_r = comp_js$r_surface)
)

write_csv(joint_summary, "threshold_sensitivity_joint_summary.csv")

dur_joint <- bind_rows(
  extract_duration_curve(fit_joint_lenient, data_joint_lenient, "Lenient (10/3/3)", shared_spp_joint),
  extract_duration_curve(baseline_fit, baseline_data, "Baseline (20/5/5)", shared_spp_joint),
  extract_duration_curve(fit_joint_strict, data_joint_strict, "Strict (30/10/5)", shared_spp_joint)
)
seas_joint <- bind_rows(
  extract_seasonal_profile(fit_joint_lenient, data_joint_lenient, "Lenient (10/3/3)", shared_species = shared_spp_joint),
  extract_seasonal_profile(baseline_fit, baseline_data, "Baseline (20/5/5)", shared_species = shared_spp_joint),
  extract_seasonal_profile(fit_joint_strict, data_joint_strict, "Strict (30/10/5)", shared_species = shared_spp_joint)
)

ggsave("figures/threshold_sensitivity_joint.pdf",
       plot_sensitivity_comparison(dur_joint, seas_joint,
         "Joint threshold sensitivity"),
       width = 11, height = 5)

rm(fit_joint_lenient, fit_joint_strict, data_joint_lenient, data_joint_strict)
gc(verbose = FALSE)

cat("  Saved: threshold_sensitivity_joint_summary.csv\n")
cat("  Saved: figures/threshold_sensitivity_joint.pdf\n")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: Window-Start Resolution Sensitivity
#             → appendix_resolution_sensitivity.md
#
# Requires pre-computed data: all_window_species_3day.rds (from Full1.R with
# step = 3 in make_window_template). If not available, section is skipped.
# The 14-day resolution is obtained by subsetting the baseline 7-day data.
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\n════ §3: Window-start resolution ════\n\n")

if (file.exists("all_window_species_3day.rds")) {
  all_ws_3day <- readRDS("all_window_species_3day.rds")
  data_3day <- prep_species_data(all_ws_3day)

  # 14-day: subset every other start position from baseline
  # day_start starts at 1 (not 0), so use modulo match on the first position
  first_day <- min(baseline_data$day_start)
  data_14day <- baseline_data |>
    filter((day_start - first_day) %% 14 == 0)
  data_14day$ds_sp_f <- droplevels(data_14day$ds_sp_f)

  fit_3day  <- fit_m6(data_3day,  label = "(3-day step)")
  fit_14day <- fit_m6(data_14day, label = "(14-day step)")

  shared_spp_res <- Reduce(intersect, list(
    levels(droplevels(baseline_data$species_f)),
    levels(droplevels(data_3day$species_f)),
    levels(droplevels(data_14day$species_f))
  ))

  comp_3d  <- compare_surfaces(baseline_fit, fit_3day, baseline_data, data_3day,
                                shared_spp_res)
  comp_14d <- compare_surfaces(baseline_fit, fit_14day, baseline_data, data_14day,
                                shared_spp_res)

  res_summary <- bind_rows(
    model_summary_row(fit_3day, "3-day step", data_3day) |>
      mutate(surface_r = comp_3d$r_surface),
    model_summary_row(baseline_fit, "7-day step (baseline)", baseline_data) |>
      mutate(surface_r = 1.0),
    model_summary_row(fit_14day, "14-day step", data_14day) |>
      mutate(surface_r = comp_14d$r_surface)
  )

  write_csv(res_summary, "resolution_sensitivity_summary.csv")

  dur_res <- bind_rows(
    extract_duration_curve(fit_3day, data_3day, "3-day", shared_spp_res),
    extract_duration_curve(baseline_fit, baseline_data, "7-day (baseline)", shared_spp_res),
    extract_duration_curve(fit_14day, data_14day, "14-day", shared_spp_res)
  )
  seas_res <- bind_rows(
    extract_seasonal_profile(fit_3day, data_3day, "3-day", shared_species = shared_spp_res),
    extract_seasonal_profile(baseline_fit, baseline_data, "7-day (baseline)", shared_species = shared_spp_res),
    extract_seasonal_profile(fit_14day, data_14day, "14-day", shared_species = shared_spp_res)
  )

  ggsave("figures/resolution_sensitivity.pdf",
         plot_sensitivity_comparison(dur_res, seas_res,
           "Resolution sensitivity: 3/7/14-day step"),
         width = 11, height = 5)

  cat("  Saved: resolution_sensitivity_summary.csv\n")
  cat("  Saved: figures/resolution_sensitivity.pdf\n")

  rm(all_ws_3day, data_3day, data_14day, fit_3day, fit_14day)
  gc(verbose = FALSE)

} else {
  cat("  SKIPPED: all_window_species_3day.rds not found.\n")
  cat("  To generate: re-run Full1.R with step = 3 in make_window_template().\n")
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Anchor Detection Parameter Sensitivity
#             → appendix_anchor_sensitivity.md
#
# Requires pre-computed data: all_window_species_strict.rds (from Full1.R
# with strict anchor parameters). Relaxed config = baseline (no binding
# parameters), so only strict is tested.
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\n════ §4: Anchor detection parameters ════\n\n")

if (file.exists("all_window_species_strict.rds")) {
  all_ws_strict <- readRDS("all_window_species_strict.rds")
  data_strict <- prep_species_data(all_ws_strict)

  fit_strict <- fit_m6(data_strict, label = "(strict anchors)")

  shared_spp_anchor <- intersect(
    levels(droplevels(baseline_data$species_f)),
    levels(droplevels(data_strict$species_f))
  )

  comp_strict <- compare_surfaces(baseline_fit, fit_strict, baseline_data,
                                   data_strict, shared_spp_anchor)

  anchor_summary <- bind_rows(
    tibble(config = "Relaxed", n_species = NA, n_rows = NA,
           dev_expl = NA, r_sq_adj = NA, surface_r = 1.000,
           note = "Identical to baseline"),
    model_summary_row(baseline_fit, "Current (baseline)", baseline_data) |>
      mutate(surface_r = 1.0, note = ""),
    model_summary_row(fit_strict, "Strict", data_strict) |>
      mutate(surface_r = comp_strict$r_surface, note = "")
  )

  write_csv(anchor_summary, "anchor_sensitivity_summary.csv")

  dur_anchor <- bind_rows(
    extract_duration_curve(baseline_fit, baseline_data, "Current (baseline)", shared_spp_anchor),
    extract_duration_curve(fit_strict, data_strict, "Strict", shared_spp_anchor)
  )
  seas_anchor <- bind_rows(
    extract_seasonal_profile(baseline_fit, baseline_data, "Current (baseline)", shared_species = shared_spp_anchor),
    extract_seasonal_profile(fit_strict, data_strict, "Strict", shared_species = shared_spp_anchor)
  )

  ggsave("figures/anchor_sensitivity.pdf",
         plot_sensitivity_comparison(dur_anchor, seas_anchor,
           "Anchor parameter sensitivity: Current vs Strict"),
         width = 11, height = 5)

  cat("  Saved: anchor_sensitivity_summary.csv\n")
  cat("  Saved: figures/anchor_sensitivity.pdf\n")

  rm(all_ws_strict, data_strict, fit_strict)
  gc(verbose = FALSE)

} else {
  cat("  SKIPPED: all_window_species_strict.rds not found.\n")
  cat("  To generate: re-run Full1.R with strict anchor parameters.\n")
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5: Independence Threshold Sensitivity
#             → appendix_independence_threshold_sensitivity.md
#
# Requires: all_window_species_15min.rds, all_window_species_60min.rds
# (from Full1.R with independence_gap = 15 or 60).
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\n════ §5: Independence threshold ════\n\n")

has_15 <- file.exists("all_window_species_15min.rds")
has_60 <- file.exists("all_window_species_60min.rds")

if (has_15 && has_60) {
  all_ws_15 <- readRDS("all_window_species_15min.rds")
  all_ws_60 <- readRDS("all_window_species_60min.rds")

  data_15 <- prep_species_data(all_ws_15)
  data_60 <- prep_species_data(all_ws_60)

  fit_15 <- fit_m6(data_15, label = "(15-min gap)")
  fit_60 <- fit_m6(data_60, label = "(60-min gap)")

  shared_spp_indep <- Reduce(intersect, list(
    levels(droplevels(baseline_data$species_f)),
    levels(droplevels(data_15$species_f)),
    levels(droplevels(data_60$species_f))
  ))

  comp_15 <- compare_surfaces(baseline_fit, fit_15, baseline_data, data_15,
                               shared_spp_indep)
  comp_60 <- compare_surfaces(baseline_fit, fit_60, baseline_data, data_60,
                               shared_spp_indep)

  indep_summary <- bind_rows(
    model_summary_row(fit_15, "15-min", data_15) |> mutate(surface_r = comp_15$r_surface),
    model_summary_row(baseline_fit, "30-min (baseline)", baseline_data) |> mutate(surface_r = 1.0),
    model_summary_row(fit_60, "60-min", data_60) |> mutate(surface_r = comp_60$r_surface)
  )

  # Extract parametric coefficients for comparison
  coef_table <- bind_rows(
    broom::tidy(fit_15, parametric = TRUE) |> filter(str_detect(term, "^(l_|s_)")) |> mutate(config = "15-min"),
    broom::tidy(baseline_fit, parametric = TRUE) |> filter(str_detect(term, "^(l_|s_)")) |> mutate(config = "30-min"),
    broom::tidy(fit_60, parametric = TRUE) |> filter(str_detect(term, "^(l_|s_)")) |> mutate(config = "60-min")
  ) |>
    select(config, term, estimate, p.value) |>
    pivot_wider(names_from = config, values_from = c(estimate, p.value))

  write_csv(indep_summary, "independence_threshold_summary.csv")

  dur_indep <- bind_rows(
    extract_duration_curve(fit_15, data_15, "15-min", shared_spp_indep),
    extract_duration_curve(baseline_fit, baseline_data, "30-min (baseline)", shared_spp_indep),
    extract_duration_curve(fit_60, data_60, "60-min", shared_spp_indep)
  )
  seas_indep <- bind_rows(
    extract_seasonal_profile(fit_15, data_15, "15-min", shared_species = shared_spp_indep),
    extract_seasonal_profile(baseline_fit, baseline_data, "30-min (baseline)", shared_species = shared_spp_indep),
    extract_seasonal_profile(fit_60, data_60, "60-min", shared_species = shared_spp_indep)
  )

  ggsave("figures/independence_threshold_sensitivity.pdf",
         plot_sensitivity_comparison(dur_indep, seas_indep,
           "Independence threshold: 15/30/60-min gap"),
         width = 11, height = 5)

  cat("  Saved: independence_threshold_summary.csv\n")
  cat("  Saved: figures/independence_threshold_sensitivity.pdf\n")

  rm(all_ws_15, all_ws_60, data_15, data_60, fit_15, fit_60)
  gc(verbose = FALSE)

} else {
  cat("  SKIPPED: pre-computed pipeline data not found.\n")
  cat("  Need: all_window_species_15min.rds, all_window_species_60min.rds\n")
  cat("  To generate: re-run Full1.R with independence_gap = 15 and 60.\n")
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6: Model Specification Diagnostics
#             → appendix_model_diagnostics.md
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\n════ §6: Model specification diagnostics ════\n\n")

diag_results <- list()

# ── 6A: k-doubling ───────────────────────────────────────────────────────────

cat("--- 6A: Basis dimension (k-doubling) ---\n")

cat("  Fitting M6 with k = (16, 12) ... ")
t0 <- proc.time()
fit_k16 <- bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(16, 12),
       by = species_f) +
    species_f + l_trapdays + l_nsites + s_latitude + s_trap_array +
    s(ds_sp_f, bs = "re"),
  data = baseline_data,
  family = Gamma(link = "log"),
  method = "fREML", discrete = TRUE, nthreads = N_THREADS,
  knots = cc_knots
)
elapsed <- round((proc.time() - t0)[3], 1)
s_k16 <- summary(fit_k16)
cat(sprintf("done (%ss, dev.expl = %.1f%%)\n", elapsed, 100 * s_k16$dev.expl))

comp_k16 <- compare_surfaces(baseline_fit, fit_k16, baseline_data, baseline_data)

diag_results$k_doubling <- tibble(
  test = "k-doubling (8,6 → 16,12)",
  dev_expl_baseline = round(100 * summary(baseline_fit)$dev.expl, 2),
  dev_expl_variant  = round(100 * s_k16$dev.expl, 2),
  aic_baseline = round(AIC(baseline_fit)),
  aic_variant  = round(AIC(fit_k16)),
  surface_r    = round(comp_k16$r_surface, 4)
)

rm(fit_k16); gc(verbose = FALSE)


# ── 6B: Gamma vs Tweedie ────────────────────────────────────────────────────

cat("\n--- 6B: Response distribution (Gamma vs Tweedie) ---\n")

cat("  Fitting M6 with Tweedie family ... ")
t0 <- proc.time()
fit_tw <- bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6),
       by = species_f) +
    species_f + l_trapdays + l_nsites + s_latitude + s_trap_array +
    s(ds_sp_f, bs = "re"),
  data = baseline_data,
  family = tw(),
  method = "fREML", discrete = TRUE, nthreads = N_THREADS,
  knots = cc_knots
)
elapsed <- round((proc.time() - t0)[3], 1)
tw_p <- fit_tw$family$getTheta(TRUE)
cat(sprintf("done (%ss, p = %.3f)\n", elapsed, tw_p))

comp_tw <- compare_surfaces(baseline_fit, fit_tw, baseline_data, baseline_data)

diag_results$tweedie <- tibble(
  test = "Gamma vs Tweedie",
  dev_expl_variant = round(100 * summary(fit_tw)$dev.expl, 2),
  estimated_p = round(tw_p, 3),
  surface_r = round(comp_tw$r_surface, 4)
)

rm(fit_tw); gc(verbose = FALSE)


# ── 6C: Nested random effects ───────────────────────────────────────────────

cat("\n--- 6C: Nested random effects (base-dataset + slice) ---\n")

# Create base-dataset × species factor
baseline_data <- baseline_data |>
  mutate(
    base_dataset = sub("_slice\\d+$", "", as.character(dataset)),
    base_ds_sp   = interaction(base_dataset, species, drop = TRUE),
    base_ds_sp_f = factor(base_ds_sp)
  )

cat("  RE levels: ds_sp_f =", nlevels(baseline_data$ds_sp_f),
    "| base_ds_sp_f =", nlevels(baseline_data$base_ds_sp_f), "\n")

cat("  Fitting M6 with nested RE ... ")
t0 <- proc.time()
fit_nested <- bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6),
       by = species_f) +
    species_f + l_trapdays + l_nsites + s_latitude + s_trap_array +
    s(base_ds_sp_f, bs = "re") + s(ds_sp_f, bs = "re"),
  data = baseline_data,
  family = Gamma(link = "log"),
  method = "fREML", discrete = TRUE, nthreads = N_THREADS,
  knots = cc_knots
)
elapsed <- round((proc.time() - t0)[3], 1)
cat(sprintf("done (%ss)\n", elapsed))

comp_nested <- compare_surfaces(baseline_fit, fit_nested, baseline_data, baseline_data)

# Extract RE SDs
vc_nested <- gam.vcomp(fit_nested)
base_sd  <- vc_nested$vc["s(base_ds_sp_f)", "std.dev"]
slice_sd <- vc_nested$vc["s(ds_sp_f)", "std.dev"]

diag_results$nested_re <- tibble(
  test = "Nested RE (base + slice)",
  dev_expl_variant = round(100 * summary(fit_nested)$dev.expl, 2),
  re_base_sd  = round(base_sd, 3),
  re_slice_sd = round(slice_sd, 3),
  delta_aic   = round(AIC(fit_nested) - AIC(baseline_fit)),
  surface_r   = round(comp_nested$r_surface, 4)
)

cat(sprintf("  Base-dataset SD: %.3f | Slice SD: %.3f | ΔAIC: %d | r: %.3f\n",
            base_sd, slice_sd, diag_results$nested_re$delta_aic,
            comp_nested$r_surface))

rm(fit_nested); gc(verbose = FALSE)


# ── 6D: Smooth BIO4 (historical — BIO4 was subsequently dropped) ────────────
#
# BIO4 was removed from the final models due to collinearity with latitude
# (r = 0.77) and LOO instability (sign-flip when BE-Leuven removed). This
# diagnostic tested whether a nonlinear BIO4 term improved fit *before* the
# decision to remove it. Retained for completeness.
#
# The baseline model in the current environment does NOT contain s_bio4.
# To reproduce this diagnostic, s_bio4 must be temporarily added.

cat("\n--- 6D: Smooth BIO4 (historical check) ---\n")

if ("s_bio4" %in% names(baseline_data)) {
  cat("  s_bio4 found in data; fitting linear and smooth variants ... \n")

  cat("  Fitting with linear s_bio4 ... ")
  t0 <- proc.time()
  fit_bio4_lin <- bam(
    abs_d_lambda ~
      te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6),
         by = species_f) +
      species_f + s_bio4 + l_trapdays + l_nsites + s_latitude + s_trap_array +
      s(ds_sp_f, bs = "re"),
    data = baseline_data,
    family = Gamma(link = "log"),
    method = "fREML", discrete = TRUE, nthreads = N_THREADS,
    knots = cc_knots
  )
  elapsed <- round((proc.time() - t0)[3], 1)
  cat(sprintf("done (%ss)\n", elapsed))

  cat("  Fitting with smooth s(s_bio4, k=5) ... ")
  t0 <- proc.time()
  fit_bio4_smooth <- bam(
    abs_d_lambda ~
      te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6),
         by = species_f) +
      species_f + s(s_bio4, k = 5) + l_trapdays + l_nsites + s_latitude + s_trap_array +
      s(ds_sp_f, bs = "re"),
    data = baseline_data,
    family = Gamma(link = "log"),
    method = "fREML", discrete = TRUE, nthreads = N_THREADS,
    knots = cc_knots
  )
  elapsed <- round((proc.time() - t0)[3], 1)
  cat(sprintf("done (%ss)\n", elapsed))

  # Extract edf of smooth term
  sm_table <- summary(fit_bio4_smooth)$s.table
  bio4_edf <- sm_table[grep("s_bio4", rownames(sm_table)), "edf"]

  diag_results$smooth_bio4 <- tibble(
    test = "Smooth BIO4 (s(s_bio4, k=5))",
    bio4_edf = round(bio4_edf, 2),
    delta_aic = round(AIC(fit_bio4_smooth) - AIC(fit_bio4_lin)),
    note = "BIO4 subsequently dropped from all models"
  )

  cat(sprintf("  Smooth edf: %.2f | ΔAIC: %d\n",
              bio4_edf, diag_results$smooth_bio4$delta_aic))

  rm(fit_bio4_lin, fit_bio4_smooth)
  gc(verbose = FALSE)

} else {
  cat("  s_bio4 not in current data (BIO4 already removed from pipeline).\n")
  cat("  To reproduce: add s_bio4 = scale(bio4_temp_seasonality) to the data.\n")

  diag_results$smooth_bio4 <- tibble(
    test = "Smooth BIO4 (historical)",
    bio4_edf = 1.00,
    delta_aic = 0,
    note = "BIO4 removed; smooth collapsed to edf=1 (linear) when tested"
  )
}


# Save diagnostics in clean tabular format
diag_summary <- tibble(
  check = c("k-doubling", "Gamma vs Tweedie", "Nested RE"),
  baseline_dev_expl = rep(diag_results$k_doubling$dev_expl_baseline, 3),
  variant_dev_expl  = c(
    diag_results$k_doubling$dev_expl_variant,
    diag_results$tweedie$dev_expl_variant,
    diag_results$nested_re$dev_expl_variant
  ),
  surface_r = c(
    diag_results$k_doubling$surface_r,
    diag_results$tweedie$surface_r,
    diag_results$nested_re$surface_r
  ),
  notes = c(
    paste0("ΔAIC = ", diag_results$k_doubling$aic_variant - diag_results$k_doubling$aic_baseline),
    paste0("p = ", diag_results$tweedie$estimated_p),
    paste0("base_ds SD=", diag_results$nested_re$re_base_sd,
           ", slice SD=", diag_results$nested_re$re_slice_sd)
  )
)

write_csv(diag_summary, "model_diagnostics_phase3.csv")
cat("\n  Saved: model_diagnostics_phase3.csv\n")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 7: Rho Filter Sensitivity
#             → appendix_rho_filter_sensitivity.md
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\n════ §7: Rho filter sensitivity ════\n\n")

# The rho_lambda model uses beta regression on rank correlations between
# sub-window and full-year detection rate rankings. Windows with few shared
# species produce artefactual boundary values (rho ≈ 1).

# Load richness data (contains rho_lambda and n_species_shared)
if (!exists("sens_richness")) {
  sens_richness <- readRDS("sensitivity_richness_data.rds")
}

# Use the n_shared_spp column from the richness data (count of species present
# in both window and FULL). This column is computed during the data pipeline.
if ("n_shared_spp" %in% names(sens_richness)) {
  sens_richness <- sens_richness |> mutate(n_shared = n_shared_spp)
} else if ("n_shared_spp" %in% names(all_window_richness)) {
  shared_counts <- all_window_richness |>
    filter(!window_id %in% c("SNAP_EU_CORE", "SNAP_EU_BUFFER", "EOW_EARLY", "EOW_LATE")) |>
    select(dataset, window_id, n_shared = n_shared_spp)
  sens_richness <- sens_richness |>
    left_join(shared_counts, by = c("dataset", "window_id"))
} else {
  # Fallback: count from all_window_species (approximate)
  shared_counts <- all_window_species |>
    filter(!window_id %in% c("SNAP_EU_CORE", "SNAP_EU_BUFFER", "EOW_EARLY", "EOW_LATE")) |>
    group_by(dataset, window_id) |>
    summarise(n_shared = n(), .groups = "drop")
  sens_richness <- sens_richness |>
    left_join(shared_counts, by = c("dataset", "window_id"))
}

# Test cutoffs — retain fits for surface correlation computation
rho_cutoffs <- c(3, 4, 5, 6, 8)
rho_results <- list()
rho_fits    <- list()
rho_datasets <- list()

for (cutoff in rho_cutoffs) {
  cat(sprintf("  Cutoff ≥ %d shared species ... ", cutoff))

  rho_data <- sens_richness |>
    filter(!is.na(rho_lambda), n_shared >= cutoff) |>
    # Rescale from [-1,1] to [0,1] then squeeze — matches models_sensitivity_surface.R
    mutate(rho_adj = pmin(pmax((rho_lambda + 1) / 2, 0.001), 0.999))

  n_obs <- nrow(rho_data)
  n_ds  <- n_distinct(rho_data$dataset)
  pct_boundary <- mean(rho_data$rho_lambda > 0.999) * 100

  if (n_obs < 100) {
    cat(sprintf("too few obs (%d), skipping\n", n_obs))
    next
  }

  t0 <- proc.time()
  fit_rho <- tryCatch(
    bam(
      rho_adj ~
        te(day_start, window_len, bs = c("cc", "tp"), k = c(12, 8)) +
        l_trapdays + l_nsites + s_latitude +
        s(dataset_f, bs = "re"),
      data = rho_data,
      family = betar(link = "logit"),
      method = "fREML", discrete = TRUE, nthreads = N_THREADS,
      knots = cc_knots
    ),
    error = function(e) {
      cat(sprintf("ERROR: %s\n", conditionMessage(e)))
      NULL
    }
  )
  elapsed <- round((proc.time() - t0)[3], 1)

  if (!is.null(fit_rho)) {
    s_rho <- summary(fit_rho)
    cat(sprintf("done (%ss, dev.expl = %.1f%%, n = %d, %d datasets, boundary = %.1f%%)\n",
                elapsed, 100 * s_rho$dev.expl, n_obs, n_ds, pct_boundary))

    rho_results[[as.character(cutoff)]] <- tibble(
      min_shared   = cutoff,
      n_obs        = n_obs,
      n_datasets   = n_ds,
      pct_boundary = round(pct_boundary, 1),
      dev_expl     = round(100 * s_rho$dev.expl, 1)
    )
    rho_fits[[as.character(cutoff)]]    <- fit_rho
    rho_datasets[[as.character(cutoff)]] <- rho_data
  }

  rm(fit_rho, rho_data)
  gc(verbose = FALSE)
}

# Surface correlations vs baseline (≥5)
baseline_key <- "5"
rho_surface_r <- numeric(length(rho_results))
names(rho_surface_r) <- names(rho_results)

if (baseline_key %in% names(rho_fits)) {
  fit_base_rho  <- rho_fits[[baseline_key]]
  data_base_rho <- rho_datasets[[baseline_key]]

  pred_grid <- tidyr::crossing(
    day_start  = seq(0, 364, by = 7),
    window_len = seq(15, 183, by = 7)
  ) |>
    mutate(
      l_trapdays = median(data_base_rho$l_trapdays),
      l_nsites   = median(data_base_rho$l_nsites),
      s_latitude = 0,
      dataset_f  = data_base_rho$dataset_f[1]
    )

  pred_base <- predict(fit_base_rho, newdata = pred_grid, type = "response",
                       exclude = "s(dataset_f)")

  for (key in names(rho_fits)) {
    if (key == baseline_key) {
      rho_surface_r[key] <- 1.0
    } else {
      grid_k <- pred_grid |>
        mutate(dataset_f = rho_datasets[[key]]$dataset_f[1])
      pred_k <- predict(rho_fits[[key]], newdata = grid_k, type = "response",
                        exclude = "s(dataset_f)")
      rho_surface_r[key] <- round(cor(as.numeric(pred_base), as.numeric(pred_k)), 3)
    }
  }
}

rho_summary <- bind_rows(rho_results) |>
  mutate(surface_r = rho_surface_r[as.character(min_shared)])

rm(rho_fits, rho_datasets)
gc(verbose = FALSE)
write_csv(rho_summary, "rho_filter_sensitivity_summary.csv")

# Figure: boundary mass vs deviance explained
p_rho1 <- ggplot(rho_summary, aes(min_shared, pct_boundary)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  labs(x = "Min shared species", y = "% at boundary (rho ≈ 1)",
       title = "Boundary mass by cutoff") +
  theme_minimal(base_size = 10)

p_rho2 <- ggplot(rho_summary, aes(min_shared, dev_expl)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  labs(x = "Min shared species", y = "Deviance explained (%)",
       title = "Model fit by cutoff") +
  theme_minimal(base_size = 10)

p_rho3 <- ggplot(rho_summary, aes(min_shared, n_obs)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  labs(x = "Min shared species", y = "Observations",
       title = "Data retention") +
  theme_minimal(base_size = 10)

p_rho4 <- ggplot(rho_summary, aes(min_shared, n_datasets)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  labs(x = "Min shared species", y = "Datasets retained",
       title = "Geographic coverage") +
  theme_minimal(base_size = 10)

ggsave("figures/rho_filter_sensitivity.pdf",
       (p_rho1 + p_rho2) / (p_rho3 + p_rho4) +
         plot_annotation(title = "Rho filter sensitivity"),
       width = 10, height = 8)

cat("  Saved: rho_filter_sensitivity_summary.csv\n")
cat("  Saved: figures/rho_filter_sensitivity.pdf\n")


# ═══════════════════════════════════════════════════════════════════════════════
# DONE
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\n═══════════════════════════════════════════════════════\n")
cat("sensitivity_appendices.R complete.\n")
cat("═══════════════════════════════════════════════════════\n\n")

cat("Output files:\n")
for (f in c("threshold_sensitivity_summary.csv",
            "threshold_sensitivity_joint_summary.csv",
            "resolution_sensitivity_summary.csv",
            "anchor_sensitivity_summary.csv",
            "independence_threshold_summary.csv",
            "model_diagnostics_phase3.csv",
            "rho_filter_sensitivity_summary.csv")) {
  cat(sprintf("  %s %s\n", ifelse(file.exists(f), "✓", "✗"), f))
}

cat("\nFigures:\n")
for (f in c("figures/threshold_sensitivity_comparison.pdf",
            "figures/threshold_sensitivity_sites_comparison.pdf",
            "figures/threshold_sensitivity_joint.pdf",
            "figures/resolution_sensitivity.pdf",
            "figures/anchor_sensitivity.pdf",
            "figures/independence_threshold_sensitivity.pdf",
            "figures/rho_filter_sensitivity.pdf")) {
  cat(sprintf("  %s %s\n", ifelse(file.exists(f), "✓", "✗"), f))
}

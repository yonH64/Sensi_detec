# ====================================================================
# appendix_study_design.R — Study design covariates analysis
# ====================================================================
#
# Tests whether camera-trap study design parameters (placement strategy,
# grid spacing, camera height, trail use) predict deviation magnitude
# beyond the baseline M6 sensitivity surface model.
#
# Requires:
#   - sensitivity_models_env.RData (fitted models + prepared data)
#   - dataset_meta.xlsx           (study design metadata)
#
# Outputs:
#   - appendix_study_design.md
#   - figures/study_design_spacing_leverage.pdf
#   - figures/study_design_nsites_gradient.pdf
#   - study_design_summary.csv
#
# NOTE: This script will need re-running once dataset_meta.xlsx is
#       complete for all datasets and the final dataset is added.
#       Imputed rows are flagged and sensitivity to their removal is tested.
# ====================================================================

library(tidyverse)
library(mgcv)
library(readxl)

# ── 0. Load data and models ──────────────────────────────────────────

load("sensitivity_models_env.RData")

mod_lambda_baseline <- all_models$comparison$M6_species
# sens_species is loaded from the environment

META_PATH <- "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets/dataset_meta.xlsx"
stopifnot(file.exists(META_PATH))

# ── 1. Load and clean design metadata ────────────────────────────────

dm_raw <- read_excel(META_PATH)

dm_clean <- dm_raw |>
  slice(-1) |>  # row 1 is a description header

  mutate(across(where(is.character), \(x) {
    x <- trimws(x)
    x[x == "" | x == " "] <- NA_character_
    x
  })) |>
  select(dataset, target, placement, grid_spacing, bait, trail, height)

# Parse numeric columns
dm_clean <- dm_clean |>
  mutate(
    grid_spacing_km = parse_number(grid_spacing),
    height_cm = case_when(
      str_detect(height, "and") ~ {
        vals <- as.numeric(str_extract_all(height, "\\d+")[[1]])
        mean(vals)
      },
      !is.na(height) ~ as.numeric(height),
      TRUE ~ NA_real_
    )
  )

# Fix SE-grimso "40 and 80" → 60 (vectorised case_when doesn't handle multi-value extraction)
dm_clean$height_cm[dm_clean$dataset == "SE-grimso"] <- 60

# ── 2. Map dataset names between xlsx and model data ─────────────────
# The xlsx uses base names; the model data uses sliced names and spatial splits.

dm_clean <- dm_clean |>
  mutate(dataset = case_when(
    dataset == "FR-monpellier" ~ "FR-montpellier",  # typo in xlsx
    TRUE ~ dataset
  ))

# Expand rows for spatial splits (same design applies to both clusters/splits)
mica_row   <- dm_clean |> filter(dataset == "NL-MICA")
grimso_row <- dm_clean |> filter(dataset == "SE-grimso")
bredvik_row <- dm_clean |> filter(dataset == "SE-bredvik")

dm_expanded <- dm_clean |>
  filter(!dataset %in% c("NL-MICA", "SE-grimso", "SE-bredvik",
                          "SI-strunjan", "SI-oljka")) |>
  bind_rows(
    mica_row |> mutate(dataset = "NL-MICA_C1"),
    mica_row |> mutate(dataset = "NL-MICA_C2"),
    grimso_row |> mutate(dataset = "SE-grimso-high"),
    grimso_row |> mutate(dataset = "SE-grimso-low"),
    bredvik_row |> mutate(dataset = "SE-Jarnashalvon")
  )

# ── 3. Impute missing design metadata ────────────────────────────────
# 5 datasets lack design info. Impute plausible values and flag them.

dm_imputed <- dm_expanded |>
  mutate(
    imputed = dataset %in% c("BE-Leuven", "GE-Mueritz_NP", "NO-nina",
                              "SP-donana", "PL-kampinos_NP",
                              "NL-MICA_C1", "NL-MICA_C2"),
    target = coalesce(target, "multiple"),
    bait   = coalesce(bait, "no"),
    placement = case_when(
      !is.na(placement) ~ placement,
      dataset == "BE-Leuven"     ~ "systematic",
      dataset == "GE-Mueritz_NP" ~ "random",
      dataset == "NO-nina"       ~ "random",
      dataset == "SP-donana"     ~ "random",
      dataset == "PL-kampinos_NP" ~ "systematic",
      dataset %in% c("NL-MICA_C1", "NL-MICA_C2") ~ "systematic",
    ),
    trail = case_when(
      !is.na(trail) ~ trail,
      dataset %in% c("BE-Leuven", "GE-Mueritz_NP", "NO-nina",
                      "SP-donana", "NL-MICA_C1", "NL-MICA_C2") ~ "off",
      TRUE ~ trail
    ),
    height_cm = case_when(
      !is.na(height_cm) ~ height_cm,
      dataset == "BE-Leuven"     ~ 50,
      dataset == "GE-Mueritz_NP" ~ 50,
      dataset == "NO-nina"       ~ 100,
      dataset == "SP-donana"     ~ 50,
      dataset == "PL-kampinos_NP" ~ 50,
      dataset %in% c("NL-MICA_C1", "NL-MICA_C2") ~ 50,
    ),
    grid_spacing_km = case_when(
      !is.na(grid_spacing_km) ~ grid_spacing_km,
      dataset == "BE-Leuven"     ~ 0.07,
      dataset == "GE-Mueritz_NP" ~ 1,
      dataset == "NO-nina"       ~ 1,
      dataset == "SP-donana"     ~ 1,
      dataset == "FR-montpellier" ~ 1,
      dataset %in% c("NL-MICA_C1", "NL-MICA_C2") ~ 1,
    )
  )

# ── 4. Join to model data ────────────────────────────────────────────

sens_design <- sens_species |>
  mutate(base_dataset = str_remove(dataset, "_slice\\d+$")) |>
  left_join(
    dm_imputed |> select(dataset, placement, grid_spacing_km,
                         trail, height_cm, imputed),
    by = c("base_dataset" = "dataset")
  ) |>
  mutate(
    # Collapse rare "mixed" levels into majority category
    placement_f = factor(case_when(
      placement == "mixed" ~ "random",
      TRUE ~ placement
    )),
    trail_f = factor(case_when(
      trail == "mixed" ~ "off",
      TRUE ~ trail
    )),
    s_height  = as.numeric(scale(height_cm)),
    s_spacing = as.numeric(scale(log(grid_spacing_km))),
    species_f = factor(species),
    ds_sp_f   = factor(paste(dataset, species, sep = "_"))
  )

stopifnot(sum(is.na(sens_design$placement_f)) == 0)

cat("Design data joined:", nrow(sens_design), "rows,",
    n_distinct(sens_design$base_dataset), "datasets,",
    n_distinct(sens_design$species), "species\n")

# ── 5. Coverage and contrast summary ─────────────────────────────────

cat("\n=== Design covariate contrast (observation-level) ===\n\n")

cat("placement:\n")
print(sens_design |> count(placement_f) |> mutate(pct = round(100 * n / sum(n), 1)))

cat("\ntrail:\n")
print(sens_design |> count(trail_f) |> mutate(pct = round(100 * n / sum(n), 1)))

cat("\nheight_cm:\n")
print(sens_design |> count(height_cm) |> mutate(pct = round(100 * n / sum(n), 1)))

cat("\ngrid_spacing_km:\n")
print(summary(sens_design$grid_spacing_km))

# ── 6. Collinearity check ────────────────────────────────────────────

cat("\n=== Correlation between design and existing covariates ===\n")
cor_mat <- sens_design |>
  distinct(base_dataset, .keep_all = TRUE) |>
  select(s_latitude, s_trap_array, s_height, s_spacing) |>
  cor(use = "complete.obs")
print(round(cor_mat, 2))

# ====================================================================
# § A: Full-data models (with RE)
# ====================================================================

cat("\n\n=== § A: Design covariates with dataset×species RE ===\n")
cat("Fitting M6 + all 4 design covariates (with RE)...\n")

t0 <- Sys.time()
mod_D1_withRE <- bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6), by = species_f) +
    species_f +
    l_trapdays + l_nsites + s_latitude + s_trap_array +
    placement_f + trail_f + s_height + s_spacing +
    s(ds_sp_f, bs = "re"),
  data = sens_design,
  family = Gamma(link = "log"),
  knots = list(day_start = c(0, 365)),
  discrete = TRUE, nthreads = 4
)
cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

cat("  Dev.expl:", round(summary(mod_D1_withRE)$dev.expl * 100, 2), "%\n")
cat("  ΔAIC vs M6 baseline:", round(AIC(mod_D1_withRE) - AIC(mod_lambda_baseline), 1), "\n")

# Extract design coefficients
coefs_withRE <- as.data.frame(summary(mod_D1_withRE)$p.table) |>
  tibble::rownames_to_column("term") |>
  filter(str_detect(term, "placement|trail|s_height|s_spacing"))
cat("\nDesign coefficients (with RE):\n")
print(coefs_withRE |> mutate(across(where(is.numeric), \(x) signif(x, 3))))

# ====================================================================
# § B: No-RE models — isolate between-dataset signal
# ====================================================================

cat("\n\n=== § B: No-RE models ===\n")

# Baseline: no design, no RE
cat("Fitting baseline (no design, no RE)...\n")
t0 <- Sys.time()
mod_noRE_base <- bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6), by = species_f) +
    species_f + l_trapdays + l_nsites + s_latitude + s_trap_array,
  data = sens_design,
  family = Gamma(link = "log"),
  knots = list(day_start = c(0, 365)),
  discrete = TRUE, nthreads = 4
)
cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

# Full design: no RE
cat("Fitting + all design (no RE)...\n")
t0 <- Sys.time()
mod_noRE_design <- bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6), by = species_f) +
    species_f + l_trapdays + l_nsites + s_latitude + s_trap_array +
    placement_f + trail_f + s_height + s_spacing,
  data = sens_design,
  family = Gamma(link = "log"),
  knots = list(day_start = c(0, 365)),
  discrete = TRUE, nthreads = 4
)
cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

# Individual covariates
fit_noRE_single <- function(extra, label) {
  f <- as.formula(paste(
    "abs_d_lambda ~",
    "te(day_start, window_len, bs = c('cc','tp'), k = c(8,6), by = species_f) +",
    "species_f + l_trapdays + l_nsites + s_latitude + s_trap_array +",
    extra
  ))
  m <- bam(f, data = sens_design, family = Gamma(link = "log"),
           knots = list(day_start = c(0, 365)), discrete = TRUE, nthreads = 4)
  s <- summary(m)
  coefs <- as.data.frame(s$p.table)
  coefs$term <- rownames(coefs)
  design_row <- coefs |> filter(str_detect(term, "placement|trail|s_height|s_spacing"))
  tibble(
    model = label,
    dev_expl = round(s$dev.expl * 100, 2),
    delta_aic = round(AIC(m) - AIC(mod_noRE_base), 0),
    beta = round(design_row$Estimate[1], 3),
    p_value = signif(design_row$`Pr(>|t|)`[1], 3)
  )
}

cat("Fitting individual design covariates (no RE)...\n")
individual_results <- bind_rows(
  tibble(model = "baseline (no design)", dev_expl = round(summary(mod_noRE_base)$dev.expl * 100, 2),
         delta_aic = 0L, beta = NA_real_, p_value = NA_real_),
  fit_noRE_single("placement_f", "+ placement"),
  fit_noRE_single("s_spacing",   "+ spacing"),
  fit_noRE_single("s_height",    "+ height"),
  fit_noRE_single("trail_f",     "+ trail"),
  tibble(model = "+ all 4 design", dev_expl = round(summary(mod_noRE_design)$dev.expl * 100, 2),
         delta_aic = round(AIC(mod_noRE_design) - AIC(mod_noRE_base), 0),
         beta = NA_real_, p_value = NA_real_)
)

cat("\nIndividual design covariate contributions (no RE, full data):\n")
print(individual_results)

# ====================================================================
# § C: Dataset-only RE + design covariates — variance decomposition
# ====================================================================

cat("\n\n=== § C: Dataset-only RE ===\n")

sens_design$dataset_f <- factor(sens_design$dataset)

cat("Fitting dataset-only RE, no design...\n")
t0 <- Sys.time()
mod_dsRE_base <- bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6), by = species_f) +
    species_f + l_trapdays + l_nsites + s_latitude + s_trap_array +
    s(dataset_f, bs = "re"),
  data = sens_design,
  family = Gamma(link = "log"),
  knots = list(day_start = c(0, 365)),
  discrete = TRUE, nthreads = 4
)
cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

cat("Fitting dataset-only RE + all design...\n")
t0 <- Sys.time()
mod_dsRE_design <- bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6), by = species_f) +
    species_f + l_trapdays + l_nsites + s_latitude + s_trap_array +
    placement_f + trail_f + s_height + s_spacing +
    s(dataset_f, bs = "re"),
  data = sens_design,
  family = Gamma(link = "log"),
  knots = list(day_start = c(0, 365)),
  discrete = TRUE, nthreads = 4
)
cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

vc_base <- gam.vcomp(mod_dsRE_base)
vc_design <- gam.vcomp(mod_dsRE_design)
re_sd_base   <- vc_base$vc["s(dataset_f)", "std.dev"]
re_sd_design <- vc_design$vc["s(dataset_f)", "std.dev"]

cat("\nDataset RE SD without design:", round(re_sd_base, 3), "\n")
cat("Dataset RE SD with design:   ", round(re_sd_design, 3), "\n")
cat("RE variance reduction:       ", round((1 - re_sd_design^2 / re_sd_base^2) * 100, 1), "%\n")

# ====================================================================
# § D: Leverage diagnostics — BE-Leuven and SI-serknica
# ====================================================================

cat("\n\n=== § D: Leverage diagnostics ===\n")

# D1: Coefficient stability without BE-Leuven
cat("Fitting no-RE model excluding BE-Leuven...\n")
sens_noBE <- sens_design |> filter(base_dataset != "BE-Leuven")

mod_noBE <- bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6),
       by = factor(species, levels = levels(sens_noBE$species_f))) +
    factor(species, levels = levels(sens_noBE$species_f)) +
    l_trapdays + l_nsites + s_latitude + s_trap_array +
    placement_f + trail_f + s_height + s_spacing,
  data = sens_noBE |> mutate(species_f = factor(species)),
  family = Gamma(link = "log"),
  knots = list(day_start = c(0, 365)),
  discrete = TRUE, nthreads = 4
)

coefs_noBE <- as.data.frame(summary(mod_noBE)$p.table) |>
  tibble::rownames_to_column("term") |>
  filter(str_detect(term, "placement|trail|s_height|s_spacing"))

cat("\nCoefficients excluding BE-Leuven:\n")
print(coefs_noBE |> mutate(across(where(is.numeric), \(x) signif(x, 3))))

# D2: Verified-only (exclude all imputed datasets)
cat("\nFitting verified-only model (no imputed, no RE)...\n")
sens_verified <- sens_design |>
  filter(!imputed) |>
  mutate(species_f = factor(species))

mod_verified <- bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6), by = species_f) +
    species_f + l_trapdays + l_nsites + s_latitude + s_trap_array +
    placement_f + trail_f + s_height + s_spacing,
  data = sens_verified,
  family = Gamma(link = "log"),
  knots = list(day_start = c(0, 365)),
  discrete = TRUE, nthreads = 4
)

coefs_verified <- as.data.frame(summary(mod_verified)$p.table) |>
  tibble::rownames_to_column("term") |>
  filter(str_detect(term, "placement|trail|s_height|s_spacing"))

cat("\nCoefficients (verified-only):\n")
print(coefs_verified |> mutate(across(where(is.numeric), \(x) signif(x, 3))))

# D3: Placement without SI-serknica
cat("\nFitting verified-only, excluding SI-serknica...\n")
sens_v_noserk <- sens_verified |>
  filter(base_dataset != "SI-serknica") |>
  mutate(species_f = factor(species))

mod_noserk <- bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6), by = species_f) +
    species_f + l_trapdays + l_nsites + s_latitude + s_trap_array +
    placement_f,
  data = sens_v_noserk,
  family = Gamma(link = "log"),
  knots = list(day_start = c(0, 365)),
  discrete = TRUE, nthreads = 4
)

coef_noserk <- as.data.frame(summary(mod_noserk)$p.table) |>
  tibble::rownames_to_column("term") |>
  filter(term == "placement_fsystematic")

cat("\nPlacement coefficient without SI-serknica:\n")
print(coef_noserk |> mutate(across(where(is.numeric), \(x) signif(x, 3))))

# ====================================================================
# § E: l_nsites scaling and residual gradient
# ====================================================================

cat("\n\n=== § E: l_nsites scaling diagnostic ===\n")

# Residuals from the no-RE baseline
sens_design$resid_noRE <- residuals(mod_noRE_base, type = "deviance")

# Residual gradient (resid_183d − resid_15d) by dataset
resid_gradient <- sens_design |>
  filter(window_len %in% c(15, 183)) |>
  group_by(base_dataset, window_len) |>
  summarise(mean_resid = mean(resid_noRE), .groups = "drop") |>
  pivot_wider(names_from = window_len, values_from = mean_resid,
              names_prefix = "resid_")

# l_nsites scaling by dataset
nsites_scaling <- sens_design |>
  group_by(base_dataset, window_len) |>
  summarise(mean_l_nsites = mean(l_nsites), .groups = "drop") |>
  filter(window_len %in% c(15, 183)) |>
  pivot_wider(names_from = window_len, values_from = mean_l_nsites,
              names_prefix = "lns_") |>
  mutate(l_nsites_range = lns_183 - lns_15)

resid_gradient <- resid_gradient |>
  left_join(nsites_scaling |> select(base_dataset, l_nsites_range),
            by = "base_dataset") |>
  mutate(gradient = resid_183 - resid_15)

r_gradient <- cor(resid_gradient$gradient, resid_gradient$l_nsites_range,
                  use = "complete.obs")
cat("Correlation between residual gradient and l_nsites scaling: r =",
    round(r_gradient, 3), "\n")

# Species-level l_nsites sensitivity
nsites_by_species <- sens_design |>
  filter(window_len == 57) |>
  group_by(species) |>
  filter(n() >= 20) |>
  summarise(
    cor_nsites_dev = cor(l_nsites, abs_d_lambda),
    mean_abs_d = mean(abs_d_lambda),
    n = n(),
    .groups = "drop"
  ) |>
  arrange(cor_nsites_dev)

cat("\nl_nsites × abs_d_lambda correlation by species (57d windows):\n")
print(nsites_by_species |> mutate(across(where(is.numeric), \(x) round(x, 3))))

# ====================================================================
# § F: Surface sensitivity — does any of this change predictions?
# ====================================================================

cat("\n\n=== § F: Surface sensitivity ===\n")

cat("Fitting l_nsites × window_len interaction model...\n")
t0 <- Sys.time()
mod_nsites_int <- bam(
  abs_d_lambda ~
    te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6), by = species_f) +
    species_f +
    l_trapdays + l_nsites + s_latitude + s_trap_array +
    l_nsites:window_len +
    s(ds_sp_f, bs = "re"),
  data = sens_design,
  family = Gamma(link = "log"),
  knots = list(day_start = c(0, 365)),
  discrete = TRUE, nthreads = 4
)
cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
cat("  ΔAIC vs baseline:", round(AIC(mod_nsites_int) - AIC(mod_lambda_baseline), 0), "\n")

# Compute surface correlations across all species
all_spp <- levels(sens_design$species_f)
surface_cors <- tibble()

for (sp in all_spp) {
  valid_level <- sens_design |>
    filter(species == sp) |>
    slice(1) |>
    pull(ds_sp_f) |>
    as.character()

  pg <- tidyr::crossing(
    day_start = seq(1, 358, by = 14),
    window_len = seq(15, 183, by = 14),
    species_f = factor(sp, levels = levels(sens_design$species_f)),
    l_trapdays = mean(sens_design$l_trapdays),
    l_nsites = mean(sens_design$l_nsites),
    s_latitude = 0,
    s_trap_array = 0,
    ds_sp_f = factor(valid_level, levels = levels(sens_design$ds_sp_f))
  )

  p1 <- predict(mod_lambda_baseline, newdata = pg, type = "response",
                exclude = "s(ds_sp_f)")
  p2 <- predict(mod_nsites_int, newdata = pg, type = "response",
                exclude = "s(ds_sp_f)")

  surface_cors <- bind_rows(surface_cors, tibble(species = sp, r = cor(p1, p2)))
}

cat("\nSurface correlation (M6 vs l_nsites interaction model):\n")
cat("  Median r:", round(median(surface_cors$r), 4), "\n")
cat("  Min r:   ", round(min(surface_cors$r), 4), "\n")

# ====================================================================
# § G: Figures
# ====================================================================

cat("\n\n=== § G: Generating figures ===\n")
dir.create("figures", showWarnings = FALSE)

# Figure 1: Spacing leverage — BE-Leuven
sens_design$resid_noRE <- residuals(mod_noRE_base, type = "deviance")

ds_resid <- sens_design |>
  group_by(base_dataset, placement_f, grid_spacing_km, imputed) |>
  summarise(mean_resid = mean(resid_noRE), n = n(), .groups = "drop")

p1 <- ggplot(ds_resid, aes(log(grid_spacing_km), mean_resid)) +
  geom_point(aes(color = placement_f, size = n)) +
  geom_text(aes(label = str_remove(base_dataset, "^[A-Z]{2}-")),
            size = 2, vjust = -1, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, color = "gray30", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "log(Grid spacing, km)", y = "Mean deviance residual",
       color = "Placement", size = "N obs") +
  theme_minimal(base_size = 10)

ggsave("figures/study_design_spacing_leverage.pdf", p1, width = 7, height = 5)
cat("  Saved figures/study_design_spacing_leverage.pdf\n")

# Figure 2: Residual gradient vs nsites scaling
p2 <- ggplot(resid_gradient, aes(l_nsites_range, gradient)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue", linewidth = 0.7) +
  geom_text(aes(label = str_remove(base_dataset, "^[A-Z]{2}-")),
            size = 2, vjust = -0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "l_nsites range (15d \u2192 183d)",
       y = "Residual gradient (resid_183d \u2212 resid_15d)") +
  theme_minimal(base_size = 10)

ggsave("figures/study_design_nsites_gradient.pdf", p2, width = 7, height = 5)
cat("  Saved figures/study_design_nsites_gradient.pdf\n")

# ====================================================================
# § H: Summary CSV
# ====================================================================

# Coefficient comparison table
coefs_full <- as.data.frame(summary(mod_noRE_design)$p.table) |>
  tibble::rownames_to_column("term") |>
  filter(str_detect(term, "placement|trail|s_height|s_spacing|s_latitude")) |>
  transmute(term, est_full = Estimate, p_full = `Pr(>|t|)`)

coefs_noBE_tbl <- coefs_noBE |>
  transmute(term, est_excl_BE = Estimate, p_excl_BE = `Pr(>|t|)`)

coefs_ver_tbl <- coefs_verified |>
  transmute(term, est_verified = Estimate, p_verified = `Pr(>|t|)`)

summary_tbl <- coefs_full |>
  left_join(coefs_noBE_tbl, by = "term") |>
  left_join(coefs_ver_tbl, by = "term")

write_csv(summary_tbl, "study_design_summary.csv")
cat("\n  Saved study_design_summary.csv\n")

cat("\n=== appendix_study_design.R complete ===\n")

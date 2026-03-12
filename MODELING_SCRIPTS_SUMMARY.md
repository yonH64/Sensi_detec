# Modeling Scripts Summary

## Overview

Four modeling scripts fit Bayesian hierarchical models comparing SNAP_EU_CORE vs SNAP_EU_BUFFER protocols. Two complementary approaches (absolute deviation and MSE), each split into detection-level and richness-level analyses. All use parallel fitting via `future`/`furrr`.

| Script | Response Type | Metrics | Family |
|--------|--------------|---------|--------|
| `models_abs_detection_PARALLEL.R` | \|estimate − truth\| | Lambda, log-rate, matched rate | Lognormal |
| `models_abs_richness_PARALLEL.R` | \|estimate − truth\|, proportions | Rarefied richness, prop_sr_full | Gaussian / Beta |
| `models_mse_detection_PARALLEL.R` | Bias² + Variance | Lambda, log-rate, matched rate | Gamma(log) |
| `models_mse_richness_PARALLEL.R` | Bias² + Variance | Rarefied richness | Gaussian |

---

## Absolute Deviation — Detection (`models_abs_detection_PARALLEL.R`)

### Model Variants per Metric (LAMBDA / RATE / MRATE)

| Name | Structure | Purpose |
|------|-----------|---------|
| `*_cov_base` | protocol + covariates + `(1\|ds_sp)` + `(1\|species)` | Baseline protocol effect |
| `*_cov_spxprot` | + `(1 + protocol_bin \|\| species)` | Species-specific protocol slopes |
| `*_cov_guild` | + `protocol_bin * guild_major` | Guild × protocol interaction |
| `*_cov_richness_xprot` | + `protocol_bin * l_n_species` | Richness moderates protocol effect |
| `*_cov_guild_minor_hab` | + `protocol_bin * guild_minor_habitat` | Habitat guild interaction |
| `*_cov_guild_minor_diet` | + `protocol_bin * guild_minor_diet` | Diet guild interaction |
| `*_mean_lambda_full` | `s_mean_lambda_full` replaces `l_n_species` | Community detectability vs richness |

**Total: 21 detection models** (7 variants × 3 metrics)

### Data Preparation

`prep_protocol_summaries()` filters to CORE/BUFFER, computes covariates (`l_trapdays`, `l_nsites`, `s_latitude`, `s_trap_array`), takes species × dataset × protocol medians, joins richness counts and community mean detectability (`mean_lambda_full`).

### Metric Definitions

- **LAMBDA** (`abs_d_lambda`): Daily detection rate from TTE, effort-free. Lognormal family.
- **Rate** (`abs_d_rate`): Raw encounter rate difference. Lognormal family.
- **MRATE** (`abs_d_matched_rate`): Daily detection rate derived from spatial coverage on matched cameras: `-log(1 - spatial_cov_m) / L`. Fully deconfounded for both denominator and accumulation effects. Lognormal family.

---

## Absolute Deviation — Richness (`models_abs_richness_PARALLEL.R`)

### Response Variables
- `abs_d_sr_raref` — absolute deviation of rarefied richness (Gaussian)
- `prop_sr_full` — proportion of FULL species detected (Beta)

### Model Variants per Response
- Base: `protocol_bin` + covariates + `(1 | dataset)`
- Dataset-slope: + `(1 + protocol_bin || dataset)`
- `mean_lambda_full` variants for each

**Total: 6 richness models**

---

## MSE — Detection (`models_mse_detection_PARALLEL.R`)

Same structure as absolute detection models but with MSE response: `(estimate − truth)² + SE²`. All use `Gamma(link = "log")`. Captures both bias and variance of estimates.

Includes base, spxprot, guild, and richness×prot variants per metric (LAMBDA / RATE / MRATE).

**Total: 12 MSE detection models** (4 variants × 3 metrics)

---

## MSE — Richness (`models_mse_richness_PARALLEL.R`)

Response: `mse_sr_raref = d_sr_raref² + sr_raref_var` (analytical rarefaction variance).

Base, dataset-slope, and `mean_lambda_full` variants.

**Total: 3 MSE richness models**

---

## Results Scripts

| Script | Input | Output |
|--------|-------|--------|
| `questions_outputs.R` | `fits_abs_detection_all.rds`, `fits_abs_richness_all.rds`, LOO objects | Q1–Q8 summaries, CSV files |
| `questions_outputs_mse.R` | `fits_mse_detection_all.rds`, `fits_mse_richness_all.rds`, LOO objects | Q1m–Q7m summaries, CSV files |

### Questions Addressed (Q1–Q8, absolute deviation)
- **Q1**: Protocol effect on absolute deviation (base models)
- **Q2**: Best model per metric (LOO-CV)
- **Q3**: Species-specific protocol variation (spxprot)
- **Q4**: Guild differences (guild models)
- **Q5**: Richness moderation (richness×prot)
- **Q6**: Covariate effects (best models, fixed effects)
- **Q7**: Random effect SDs (grouping structure variance)
- **Q8**: Community detectability vs. species richness

---

## Parallelism

All PARALLEL scripts use `future::plan(multisession, workers = 4)` with `CORES_PER_MODEL = 4` (16 cores total). Each model runs 4 chains × 4000 iterations (2000 warmup), `adapt_delta = 0.98`, `max_treedepth = 12`.

---

## Superseded Scripts

| Old Script | Replaced By |
|------------|-------------|
| `models_abs_UPDATED.R` | `models_abs_detection_PARALLEL.R` + `models_abs_richness_PARALLEL.R` |
| `models_mse_UPDATED.R`, `models_mse.R` | `models_mse_detection_PARALLEL.R` + `models_mse_richness_PARALLEL.R` |
| `models_abs.R`, `modeling.R` | Earlier iterations |

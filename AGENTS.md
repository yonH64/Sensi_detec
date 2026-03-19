# Sensitivity of Detection — Project Memory

## Project Overview

PhD project evaluating how camera-trap temporal sampling window design affects species detection metrics. Uses a **sensitivity surface approach**: a sliding window grid (25 durations [15–183 days] × 53 start positions, stepped every 7 days across the full year) quantifies detection metric deviation from a 12-month benchmark as a function of window timing, duration, species identity, and environmental context. The Snapshot Europe CORE (61d, Sep 1–Oct 31) and BUFFER (89d, Aug 18–Nov 14) protocols are evaluated as specific points on this surface, not as the primary comparison. Models are fitted as GAMMs (`mgcv::bam()`) with cyclic splines for day-of-year circularity and species-specific 2D tensor product surfaces.

## Key Data Files

| File | Description |
|------|-------------|
| `all_window_species.rds` | Core dataset: detection metrics per species × window × dataset. Primary deviation columns: `d_lambda`, `d_rate`, plus matched-camera columns (see Metric Inventory below). |
| `all_window_richness.rds` | Community-level richness metrics per window × dataset. |
| `all_dropped_species.rds` | Species that failed minimum viable sample thresholds (`min_events`, `min_sites_pos`, `min_occasions_pos`). |
| `species_meta.csv` | Species taxonomy: `species`, `guild_major`, `guild_minor_habitat`, `guild_minor_diet`. 63 species. |
| `dataset_metadata.csv` | Dataset-level metrics + environmental covariates + camera setup metadata. Written by `dataset_report.R`. |
| `dataset_meta.xlsx` | (in Datasets/) Provider, country, camera setup metadata per dataset. Joined by `dataset_report.R`. |
| `wrapped*.rds` (~8–9 MB each) | Intermediate outputs from `dataset_wrapper()`. |
| `sensitivity_species_data.rds` | Prepared species-level modeling data for sensitivity surface analysis. 222,748 rows, 29 species × window × 35 datasets. |
| `sensitivity_richness_data.rds` | Prepared richness-level modeling data. 46,375 rows (window × 35 datasets including sliced datasets). |
| `sensitivity_gam_models.rds` | All fitted GAM models via `mgcv::bam()`. Contains: 8 model comparison variants for lambda (M1–M6 + M_hab + M_diet), M6-structure detection models for 4 metrics, and 3 community models (richness, proportion, rank correlation). |
| `sensitivity_models_env.RData` | Full R environment saved after model fitting. Contains `all_models`, `sens_species`, `sens_richness`, and `all_window_species`/`all_window_richness` (minus intermediate fitting objects). Loaded by `sensitivity_results.R` and `sensitivity_figures.R`. |
| `interannual_cv_lambda_full.csv` | Inter-annual CV of lambda_full for 33 species × site combinations across 4 multi-year sites. Used to quantify benchmark noise floor. |
| `benchmark_lambda_180_270.rds` | TTE lambda computed for 180d and 270d windows at every 7-day start position across all 23 dataset-slices. Used for benchmark robustness check. |
| `robustness_check_data.rds` | Existing sub-windows matched to their closest-centred 180d and 270d benchmarks, with deviations computed relative to each benchmark. |
| `robustness_benchmark_summary.csv` | Summary table: shape correlations, mean deviations, and sign statistics by window length × benchmark duration. |
| `BIG_PICTURE_Pager_SamplingWindow.Rmd` | Original pager (Option A framing: deviation from 365d benchmark). |
| `BIG_PICTURE_Pager_OptionB.Rmd` | Pager reframed as "seasonal detection variation" with application-specific guidance (standardised monitoring, occupancy, trend detection, community comparison). |
| `BIG_PICTURE_Pager_OptionC.Rmd` | Pager with dual-target approach: models both bias (|d|) and MSE (bias² + variance) as parallel surfaces. Preferred framing. |

### Old brms model files (from superseded protocol comparison, kept for reference)

| File | Description |
|------|-------------|
| `fits_abs_all.rds` (559 MB) | All fitted `brms` model objects (absolute deviation + richness models). |
| `loo_abs_all.rds` (165 KB) | LOO-CV results for model comparison. |
| `fits_env_covariate_comparison.rds` | 15 brms models comparing 5 environmental covariates × 3 detection metrics. |
| `loo_env_covariate_comparison.rds` | LOO-CV results for the 15 env covariate models. |
| `fits_env_bio4_quadratic.rds` | 3 brms models testing quadratic temperature seasonality. |
| `fits_env_bio4_bio15_combined.rds` | 3 brms models with both temp seasonality + precip seasonality. |

## Script Execution Order

### 1. Data Pipeline
- **`helpers.R`** — **Canonical shared utilities.** Source this at the top of any standalone script. Contains: `resolve_dataset_path`, `dataset_output_dir`, `read_csv_src`, `set_camtrap_cols`, `parse_ts_safe`, `parse_num_safe`, `pad_bbox_safe`, `daily_deployment_effort_fast` (vectorised, preferred), `daily_deployment_effort` (legacy, used by `dataset_wrapper.R`), `monthly_deployment_effort`, `taxa_filter`, `make_window_template`.
- **`Full1.R`** — Main orchestration script. Sources `helpers.R` at the top. Defines `add_snapshot_europe_windows`, `compute_full_effort`, `split_spatial_by_distance`, `find_anchors` (renamed from `slice_auto_seasons`), `build_window_metrics_fast1`, and `dataset_wrapper1`. Runs the full pipeline (anchor detection → window metrics → saving results).
- **`build_window_metrics.R`** — Calculates `p_tte`, `log_rate`, `spatial_cov` per species × window. Called by `dataset_wrapper()`.
- **`dataset_wrapper.R`** — Batch-processes all datasets through the window metric pipeline.

### 2. Sensitivity Surface Analysis (current)
- **`prep_sensitivity_data.R`** — Prepares species-level and richness-level data for the sensitivity surface analysis. Joins species traits, environmental covariates (via `join_env_covariates()` with slice-name matching), and standardises covariates. Outputs `sensitivity_species_data.rds` and `sensitivity_richness_data.rds`.
- **`models_sensitivity_surface.R`** — All GAM/GAMM models for the sensitivity surface. Uses `mgcv::bam()` with cyclic splines and species-specific tensor product surfaces. Fits 8 model comparison variants (M1–M6 + M_hab + M_diet) for lambda, plus M6 for all other detection metrics, plus 3 community models. Outputs `sensitivity_gam_models.rds` and `sensitivity_models_env.RData` (cleaned environment for downstream scripts).
- **`sensitivity_results.R`** — Extracts derived quantities from fitted models. Loads `sensitivity_models_env.RData`; creates model shorthand aliases (`mod_lambda`, `mod_rate`, `mod_matched`, `mod_rate_signed`, `mod_richness`, `mod_prop`, `mod_rho`). Answers Q1–Q8 (duration effect, seasonal profiles, full surface predictions, species/guild variation, BIO4 effect, protocol evaluation, optimal timing, richness recovery). Q9 computes the benchmark noise floor (inter-annual CV of lambda_full from multi-year sites, SNR = |d_lambda| / SD_benchmark per observation). Outputs Q1–Q9 CSV files + `model_comparison_table.csv`. **Note:** Uses `tidyr::crossing()` explicitly to avoid namespace conflict with `igraph::crossing()`.
- **`sensitivity_figures.R`** — Publication figures (9 main + 1 supplementary PDF): main surface, guild/species surfaces, duration curves, signed rate surface, richness surfaces, protocol evaluation, model comparison, benchmark noise floor (4-panel: CV by species, overall SNR curve, % below noise floor, species-specific SNR curves), and Fig S1 benchmark robustness check (3-panel: shape correlations across benchmarks, deviation magnitude convergence, seasonal profiles at 29d/85d for 180d/270d/365d benchmarks). Loads `sensitivity_models_env.RData` + Q CSV files + `robustness_check_data.rds` + `robustness_benchmark_summary.csv`.

### 2b. Model Fitting — OLD protocol comparison (superseded by sensitivity surface)
- **`models_abs_detection_PARALLEL.R`** — Absolute deviation models for 3 detection metrics. All lognormal family.
- **`models_abs_richness_PARALLEL.R`** — Absolute deviation models for richness metrics. Gaussian/Beta families.
- **`models_mse_detection_PARALLEL.R`** — MSE models for 3 detection metrics. Gamma(log) family.
- **`models_mse_richness_PARALLEL.R`** — MSE model for rarefied richness.

### 3. Results & Interpretation
- **`sensitivity_results.R`** — See §2 above.
- **`questions_outputs.R`** — **SUPERSEDED.** Old protocol-comparison results extraction (brms, CORE vs BUFFER). Kept for reference only. All old Q CSV files have been removed.
- **`questions_outputs_mse.R`** — **SUPERSEDED.** Same structure for MSE model results.

### 4. Visualization & Reporting
- **`Full_visu.R`** — Main visualization script: dropped-window checks, calendar heatmaps, metric distributions, protocol comparisons.
- **`dataset_report.R`** — Unified per-dataset diagnostic report. Produces one PDF per dataset with: dashboard (metrics + map + effort sparkline), effort timeline with anchor slices & Snapshot Europe windows, species × month detection heatmap, per-slice richness maps & accumulation curves, individual species panels, rejected species table. Also extracts WorldClim + MODIS environmental covariates and joins `dataset_meta.xlsx`. Stores enriched metrics in `dataset_overview_metrics` and saves `dataset_metadata.csv`.
- **`visualize_suppression_effect.R`** — Publication figures for the latitude effect. **Note:** Originally built around the suppression effect narrative (old `d_spatial_cov` metric); that pattern does not hold for the deconfounded metrics. Script needs updating to reflect the direct latitude effect.
- **`portfolio_effect_comparison.R`** — Portfolio effect analysis and figures. **Note:** Dataset-level portfolio effect (species richness predicting approximation quality) is absent for the deconfounded metrics; script may need revision.


## Metric Inventory

### Detection metrics (per species × window × dataset)

**Question answered:** "How well does this sub-window capture species-level detection compared to the full year?"

#### Core metrics (computed per window, all cameras)

| Column | Formula | Description |
|--------|---------|-------------|
| `lambda` | `n_first_detections / total_exposure` | Daily detection rate from TTE (time-to-first-event per camera). Effort-free. |
| `rate` | `n_events / trap_days` | Encounter rate (events per trap-day). Effort-free. |
| `log_rate` | `log(rate + 1e-6)` | Log-transformed encounter rate. |
| `spatial_cov` | `n_sites_pos / n_sites` | Proportion of active cameras that detected the species. **Confounded** with window length (see below). Retained as a descriptive statistic; not used for deviations. |

#### Deviations from FULL benchmark (all cameras)

| Column | Formula | Description |
|--------|---------|-------------|
| `d_lambda` | `lambda_window − lambda_full` | TTE daily detection rate deviation. **Clean** — no structural confound. |
| `d_rate` | `rate_window − rate_full` | Raw encounter rate deviation. **Clean** — no structural confound. Replaced `d_log_rate` which had a log-amplification confound at low baseline rates (see Structural Confounds for rationale). |
| `mse_lambda` | `d_lambda² + lambda_se²` | MSE criterion (bias² + variance). |
| `mse_rate` | `d_rate² + rate_se²` | MSE criterion. |

Note: `d_lambda` and `d_rate` are moderately correlated (r ≈ 0.21) but capture partly different aspects of detection (TTE-derived vs event-based).

#### Matched-camera metrics (cameras active in BOTH window and FULL)

The matched-camera approach restricts all computations to cameras that were active in both the sub-window and the FULL 12-month window. This fixes the **denominator confound** (different camera pools) and enables fair spatial coverage comparisons.

| Column | Formula | Description |
|--------|---------|-------------|
| `n_matched` | `length(intersect(cams_window, cams_full))` | Number of shared cameras. |
| `spatial_cov_m` | `n_detected_shared / n_matched` | Spatial coverage on matched cameras (window side). |
| `spatial_cov_full_m` | `n_detected_shared_full / n_matched` | Spatial coverage on matched cameras (FULL side). |
| `matched_rate` | `-log(1 − spatial_cov_m) / L` | Daily detection rate derived from matched spatial_cov. Deconfounds **both** denominator and accumulation effects. |
| `matched_rate_full` | `-log(1 − spatial_cov_full_m) / 365` | Same, FULL side. |
| `d_matched_rate` | `matched_rate − matched_rate_full` | **Fully deconfounded spatial-coverage-based detection rate deviation.** |
| `rate_m` / `rate_full_m` | `events_shared / effort_shared` | Encounter rate on matched cameras. |
| `d_rate_m` | `rate_m − rate_full_m` | Matched encounter rate deviation (raw difference). Controls for camera composition. |
| `mse_matched_rate` | `d_matched_rate² + matched_rate_se²` | MSE for matched detection rate. |
| `mse_rate_m` | `d_rate_m² + rate_m_se²` | MSE for matched encounter rate. |
| `cam_det_ids` | (list column) | Camera IDs that detected this species in this window. |
| `cam_rate_tbl` | (list column) | Per-camera events/effort/rate table. |

### Richness metrics (per window × dataset)

**Question answered:** "How well does this sub-window capture community-level diversity compared to the full year?"

| Column | Formula | Description |
|--------|---------|-------------|
| `sr_obs` | count of species passing thresholds | Observed richness. |
| `sr_raref` | Incidence-based rarefaction (Colwell 2012) | Richness standardised to common sampling intensity. |
| `d_sr_raref` | `sr_raref_window − sr_raref_full` | Rarefied richness deviation. |
| `prop_sr_full` | &#124;spp_window ∩ spp_full&#124; / &#124;spp_full&#124; | Proportion of FULL species recovered. |
| `rho_lambda` | Spearman ρ across species | Rank preservation of detection rate. |
| `rho_log_rate` | Spearman ρ across species | Rank preservation of encounter rate. |
| `rho_spatial_cov` | Spearman ρ across species | Rank preservation of spatial coverage. |

### Design / complementary variables

| Column | Description |
|--------|-------------|
| `trap_days_window` / `trap_days_full` | Total trap-days in window / FULL. |
| `n_sites` / `n_sites_full` | Active cameras in window / FULL. |
| `n_events_total` | Independent events for this species. |
| `n_sites_pos` | Cameras detecting this species. |
| `n_occasions_pos` | Occasion periods with detections. |
| `latitude` | Centroid latitude of camera array. |
| `trap_array` | Max pairwise camera distance (km). |
| `window_len` | Window duration (days). |

### Environmental covariates (per dataset, from `dataset_report.R`)

| Column | Source | Description |
|--------|--------|-------------|
| `bio1_mean_temp` | WorldClim v2.1 | Mean annual temperature (°C × 10). |
| `bio4_temp_seasonality` | WorldClim v2.1 | Temperature seasonality (SD × 100). **Best environmental predictor of detection deviation** — outperforms raw latitude by LOO-CV. |
| `bio5_tmax_warmest` | WorldClim v2.1 | Max temperature of warmest month. |
| `bio6_tmin_coldest` | WorldClim v2.1 | Min temperature of coldest month (winter severity). |
| `bio12_annual_precip` | WorldClim v2.1 | Annual precipitation (mm). |
| `bio15_precip_seasonality` | WorldClim v2.1 | Precipitation seasonality (CV). |
| `ndvi_amplitude` | MODIS MOD13Q1 | Annual NDVI max − min (vegetation phenology amplitude). |
| `ndvi_cv` | MODIS MOD13Q1 | Annual NDVI coefficient of variation. |
| `daylength_range_hr` | Computed from latitude | Longest − shortest day (hours). Purely latitudinal. |

### Camera setup metadata (per dataset, from `dataset_meta.xlsx`)

| Column | Description |
|--------|-------------|
| `provider` | Data provider name. |
| `country` | Country of camera array. |
| `target` | Target species focus (single/multi-species). |
| `bait` | Whether bait was used. |
| `trail` | Camera placement (on/off trail). |
| `spacing` | Camera spacing (km). |
| `height` | Camera height (cm). |


## Structural Confounds Identified

### `p_tte` was structurally confounded with window length (REMOVED)

`p_tte = 1 − exp(−λ × L)` is mathematically guaranteed to be smaller for shorter windows, even if the daily rate λ is identical. The deviation `d_p_tte = p_tte_window − p_tte_full` was 75–79% negative purely from this mathematical scaling. **Replaced by `lambda`** (= `−log(1 − p_tte) / L`), which extracts the window-length-independent daily rate.

### `spatial_cov` has two confounds

1. **Denominator confound**: The camera pool (`n_sites`) changes between windows — sub-windows have fewer active cameras than FULL. This can push `spatial_cov_full` down (dilution by cameras that never detect), partially masking the accumulation effect.
2. **Accumulation confound**: More observation time → more cameras detect the species. This pushes `spatial_cov_full` up relative to sub-windows.

These two confounds act in **opposite directions** and partially cancel in the unmatched `d_spatial_cov` (76% negative). Matching cameras fixes confound #1, exposing confound #2 fully (93% negative). The `-log(1-scov)/L` transformation (→ `matched_rate`) removes both confounds.

**Raw `d_spatial_cov` and `d_spatial_cov_m` are no longer computed.** Use `d_matched_rate` for a deconfounded spatial-coverage-based detection metric.

### `d_log_rate` had a log-amplification confound (REPLACED by `d_rate`)

The previous metric `d_log_rate = log(rate_window + 1e-6) − log(rate_full + 1e-6)` amplified small absolute differences at low baseline encounter rates. For the same absolute rate difference (0.005 events/trap-day), `d_log_rate` was ~5.5× larger at a baseline of 0.02 than at 0.10. Since northern sites have systematically lower encounter rates (median ~0.02 vs ~0.04 south), this produced a **positive** latitude coefficient for |d_log_rate| — the opposite direction from lambda and matched rate.

**Replaced by `d_rate = rate_window − rate_full`**, the raw absolute difference in encounter rate. This has no floor amplification (r = +0.55 with baseline rate, same direction as d_lambda), and produces a negative latitude coefficient consistent with the other two metrics.

### 12-month benchmark defensibility

The 12-month benchmark is defensible as the **operational reference** ("what you'd get with year-round monitoring"), not as ecological "truth." Inter-annual variability in the benchmark is genuine ecological variation, not measurement error. Within each year-slice, sub-windows and FULL share the same ecological conditions, so the within-year deviation is clean.

#### Inter-annual CV of lambda_full (benchmark noise floor)

Quantified from the 4 multi-year sites (BE-Leuven 5 yr, NO-evenstadlia 2, SI-serknica 2, SP-donana 4): 33 species × site combinations.

| Statistic | Value |
|-----------|-------|
| Median CV | 23% |
| IQR | 13–34% |
| Range | 3.8% (Ursus arctos, SI-serknica) – 87% (Cervus elaphus, NO-evenstadlia, 2 yr only) |

CV is weakly correlated with detectability (r = −0.30 with log lambda_full): rare species tend to be noisier, but with substantial scatter.

#### Signal-to-noise ratio (deviation vs benchmark noise)

SNR = |d_lambda| / inter-annual SD(lambda_full). Benchmark noise is not a practical limitation for most of the surface:

| Window length | Median SNR | % observations with SNR < 1 |
|---------------|-----------|------------------------------|
| 15 d | 90 | 0% |
| 57 d | 25 | 0.03% |
| 85 d | 14 | 0.5% |
| 120 d | 8.2 | 2.3% |
| 155 d | 5.0 | 7.9% |
| 183 d | 3.4 | 13.7% |

For the handful of high-CV species at specific sites (Alces alces at NO-evenstadlia, Sciurus vulgaris at SI-serknica), deviations at long windows (≥85 d) fall below benchmark noise for 10–16% of observations. For most species, the SNR = 1 threshold is not reached even at 120 d.

**"Deviation smaller than benchmark noise" could serve as a principled stopping rule** for the diminishing-returns duration question: the point beyond which lengthening the window doesn't improve the signal relative to inter-annual fluctuation.

Saved as `interannual_cv_lambda_full.csv`.

#### Benchmark robustness check (180d and 270d centred benchmarks)

TTE lambda was recomputed for 180-day and 270-day windows at every 7-day start position across all 23 dataset-slices, then matched to each existing sub-window by closest centre. The seasonal shape of the deviation surface was compared across benchmarks.

**Seasonal shape correlation (|d_lambda| profile across start positions, 365d vs alternative):**

| Sub-window length | 365d vs 180d | 365d vs 270d |
|-------------------|-------------|-------------|
| 15 d | 0.995 | 0.998 |
| 29 d | 0.983 | 0.995 |
| 57 d | 0.970 | 0.851 |
| 85 d | 0.933 | 0.602 |
| 120 d | 0.879 | 0.705 |

The 180d benchmark preserves the surface shape very well (r > 0.88 for all window lengths). The 270d benchmark degrades for longer sub-windows because there is little contrast left (e.g., 120d window vs 270d benchmark). This degradation is mechanical, not evidence of benchmark-dependence.

**Conclusion:** The sensitivity surface shape is robust to benchmark choice. The autumn rut spike and the duration smoothing curve appear regardless of whether the benchmark is 180, 270, or 365 days.

Saved as `benchmark_lambda_180_270.rds`, `robustness_check_data.rds`, `robustness_benchmark_summary.csv`. Figure: `figures/FigS1_robustness_benchmark.pdf` (generated by `sensitivity_figures.R`; replaces earlier ad-hoc `robustness_benchmark_comparison.pdf`).


## Model Structure — Sensitivity Surface (current)

### Approach: GAMMs via `mgcv::bam()`

Species-specific 2D tensor product surfaces model how detection metric deviation varies jointly as a function of window start day (circular) and window duration, with species-level surfaces, environmental covariates, effort controls, and dataset × species random effects.

### Best model formula (species-level detection, M6)

```r
abs_d_metric ~
  te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6), by = species_f) +
  species_f +
  s_bio4 + l_trapdays + l_nsites + s_latitude + s_trap_array +
  s(ds_sp_f, bs = "re")
family = Gamma(link = "log")
knots = list(day_start = c(0, 365))
```

### Model comparison (lambda metric, 222,748 rows, 29 species)

| Model | Structure | AIC | Dev. explained |
|-------|-----------|-----|----------------|
| M1 | Shared surface + covariates | -2,217,140 | 84.6% |
| M2 | + guild-varying seasonality | -2,230,664 | 85.4% |
| M3 | Full guild × surface | -2,232,413 | 85.5% |
| M4 | + bio4 × duration | -2,220,211 | 84.7% |
| M5 | Guild surface + bio4 × duration | -2,235,716 | 85.7% |
| M_hab | Minor habitat guild surfaces (9 levels) | -2,235,802 | 85.7% |
| M_diet | Minor diet guild surfaces (10 levels) | -2,239,747 | 86.0% |
| **M6** | **Species-specific surfaces (29 species)** | **-2,258,933** | **87.1%** |

Species-specific surfaces (M6) dominate: ΔAIC ≈ -19,186 over the best guild model (M_diet). The penalization automatically regularises rare species.

### Performance across metrics (M6 structure)

| Response | Family | N | Dev. explained |
|----------|--------|---|----------------|
| `abs_d_lambda` | Gamma(log) | 222,749 | 87.1% |
| `abs_d_matched_rate` | Gamma(log) | 222,749 | 86.9% |
| `abs_d_rate` | Gamma(log) | 222,749 | 71.4% |
| `d_rate` (signed) | Gaussian | 222,749 | 19.0% |
| `prop_sr_full` | Beta(logit) | 46,376 | 77.3% |
| `d_sr_raref` | Gaussian | 46,376 | 24.7% |
| `rho_lambda` | Beta(logit) | 23,944 | 36.3% |

### Key design choices
- **Cyclic spline** `bs = "cc"` for `day_start` with `knots = list(day_start = c(0, 365))`: handles year-wrapping windows (~19% of data).
- **Tensor product** `te()`: allows the timing × duration interaction to be non-separable (short windows are more timing-sensitive than long ones).
- **`by = species_f`**: separate surface per species (29 species). Penalization controls complexity: species with few observations get nearly flat surfaces.
- **`discrete = TRUE, nthreads = 4`** in `bam()`: fast fitting for the species-level models.
- **Window range**: 15–183 days (6 months) in 7-day steps, covering 25 duration levels. Extended from original 120-day cap to characterize the full plateau region.

### Richness models (dataset-level, shared surface)

```r
d_sr_raref ~ te(day_start, window_len, bs = c("cc", "tp"), k = c(12, 8)) +
  s_bio4 + l_trapdays + l_nsites + s_latitude + s(dataset_f, bs = "re")
family = gaussian()
```

## Research Questions — Sensitivity Surface (current)

| Q# | Question | Approach |
|----|----------|----------|
| Q1 | How does deviation depend on window duration? | Marginal effect of `window_len` from the tensor surface |
| Q2 | How does deviation depend on window timing (season)? | Marginal effect of `day_start` from the tensor surface |
| Q3 | How do duration and timing interact? | Full 2D tensor surface; contrast short vs long windows across seasons |
| Q4 | How do species traits modulate the surface? | Species-specific surfaces (M6); aggregate by guild for interpretation |
| Q5 | How does temperature seasonality modulate the surface? | `s_bio4` coefficient; optional `s(window_len, by = s_bio4)` interaction |
| Q6 | Where do named protocols sit on the surface? | Derived: predict at CORE (day 244, 61d), BUFFER (day 230, 89d), EOW_EARLY (day 214, 60d), EOW_LATE (day 274, 60d) |
| Q7 | What is the optimal window design? | Derived: find timing that minimises predicted deviation at given durations |
| Q8 | How does species richness recovery vary across the surface? | Richness models (prop_sr_full, d_sr_raref, rho_lambda) |
| Q9 | At what duration does deviation fall below benchmark noise? | SNR = \|d_lambda\| / inter-annual SD(lambda_full); species-specific noise-floor thresholds from multi-year sites |

### Q output files (all from `sensitivity_results.R`)

| File | Content |
|------|---------|
| `Q1_duration_effect.csv` | Mean/median/IQR predicted deviation by window_len × guild |
| `Q2_seasonal_profiles.csv` | Mean deviation by day_start at fixed durations (15, 30, 60, 90, 120d) |
| `Q3_surface_predictions.csv` | Full 2D surface predictions per metric × guild |
| `Q4_species_guild_surfaces.csv` | Species-level deviation by window_len × season (29 species) |
| `Q5_bio4_effect.csv` | BIO4 parametric coefficients per metric |
| `Q6_protocol_evaluation.csv` | Predicted deviation for 4 named protocols × guild × metric |
| `Q7_optimal_timing.csv` | Best/worst timing per species at key durations |
| `Q8_richness_surface.csv` | Richness surface predictions (d_sr_raref, prop_sr_full, rho_lambda) |
| `Q9_noise_floor.csv` | SNR by window duration (overall) + per-species noise-floor summary at long windows |
| `model_comparison_table.csv` | AIC/deviance comparison for M1–M6 + M_hab + M_diet models |

### Figure output files (all from `sensitivity_figures.R`)

| File | Content |
|------|---------|
| `Fig1_sensitivity_surface.pdf` | Main 2D surface (all species pooled, |d_lambda|) with protocol positions |
| `Fig2_guild_surfaces.pdf` | Guild-specific surfaces (5 guilds excl. Insectivore) |
| `Fig3_species_surfaces.pdf` | Focal species surfaces (6 species, re-predicted from M6) |
| `Fig4_duration_curves.pdf` | Duration-deviation curves by season |
| `Fig5_signed_rate_surface.pdf` | Signed encounter rate surface (bidirectional) |
| `Fig6_richness_surfaces.pdf` | Three-panel: species recovered, rank preservation, richness deviation |
| `Fig7_protocol_evaluation.pdf` | Protocol × guild × metric comparison (CORE/BUFFER/EOW) |
| `Fig8_model_comparison.pdf` | AIC ladder (M1–M6 + M_hab + M_diet) |
| `figures/Fig9_benchmark_noise_floor.pdf` | 4-panel: CV by species, overall SNR, % below noise floor, species-specific SNR curves |
| `figures/FigS1_robustness_benchmark.pdf` | Supplementary: 3-panel benchmark robustness check (shape correlations, deviation magnitude, seasonal profiles for 180d/270d/365d) |
| `figures/residual_diagnostics_by_duration.pdf` | Residual diagnostics stratified by window duration bin (15–60d, 61–120d, 121–183d): residual vs fitted, QQ plots, density |

## Research Questions — OLD protocol comparison (superseded)

| Q# | Question | Source Models |
|----|----------|---------------|
| Q1 | Does protocol affect absolute deviation? | Base models |
| Q2 | Which model fits best per metric? | LOO-CV comparison |
| Q3 | Does effect vary by species? | spxprot models |
| Q4 | Does effect differ by guild? | Guild models |
| Q5 | Does richness moderate protocol effect? | Richness×prot models |
| Q6 | Which covariates drive absolute deviation? | Best models, fixed effects |
| Q7 | How much variance in grouping structure? | Random effect SDs |
| Q8 | Does community detectability explain deviation better than richness? | mean_lambda_full models |

## Key Findings

### Sensitivity surface findings (current)

- **The sensitivity surface is dominated by species identity**: Species-specific 2D surfaces (M6) dramatically outperform guild-level (ΔAIC ≈ -19,186) and shared surfaces. Species with strong seasonal activity patterns (ungulate rut, hibernation) have much more pronounced surfaces than generalists.
- **Ungulate autumn rut is the primary driver of high deviations**: Cervus elaphus and Sus scrofa show intense deviation spikes for windows centered on Sep–Nov, driven by rutting-season detection rate inflation. Capreolus capreolus shows both spring and autumn peaks. Vulpes vulpes and Lepus europaeus have near-flat surfaces.
- **Duration effect is steep then plateauing**: Deviation drops ~80% from 15→60 days, then diminishes. Overall mean |d_lambda| goes from 0.026 (15d) → 0.005 (60d) → 0.002 (120d) → 0.001 (183d).
- **Short windows are most timing-sensitive**: A 15-day window on the rut captures a very different picture than in summer. Windows >75 days are relatively robust to timing choice.
- **Temperature seasonality (BIO4) is significant** as a linear environmental predictor for most detection metrics (β ranges from −0.16 to −0.23 on the Gamma link scale, p < 0.007 for abs_d_rate and abs_d_matched_rate; p = 0.069 marginal for abs_d_lambda). Higher seasonality → smaller deviations after accounting for species and window design.
- **Signed encounter rate deviation is bidirectional**: d_rate is ~59% positive / 41% negative. Summer/autumn windows overestimate encounter rate; winter/spring underestimate. d_lambda and d_matched_rate are ~99% and 95% positive respectively.
- **Richness recovery favours spring**: `prop_sr_full` (proportion of full-year species recovered) peaks for spring-centered windows (best start ≈ day 99) and increases monotonically with duration. Predicted recovery reaches ~70% at 64 days.
- **Rank preservation is uniformly high**: rho_lambda (rank correlation of detection rates between window and FULL) reaches ~0.96 at 64 days across the predicted surface (model restricted to windows with ≥5 shared species, N = 12,028; see methods_note_rho_filter.md). The rut inflates all ungulate rates proportionally, preserving ranking even though absolute values are biased.
- **Snapshot Europe CORE sits on the edge of the autumn hot zone**: Both CORE (day 244, 61d) and BUFFER (day 230, 89d) are positioned in the region of elevated deviation due to the autumn activity spike. At 60d, optimal start ≈ day 344 (early December) for lambda, day 309 (early November) for encounter rate. EOW_EARLY (day 214, 60d) and EOW_LATE (day 274, 60d) split the autumn window at Oct 1; for lambda and matched_rate, all three 60-day autumn protocols show similar predicted deviation — the rut signal spans the full autumn. BUFFER benefits from its longer duration (89d), not its timing.

### Previous findings (from old protocol comparison, retained for reference)

- **Latitude effect (direct, no suppression)**: With the clean metrics (`d_lambda`, `d_rate`, `d_matched_rate`), latitude is a direct, significant predictor of deviation magnitude in simple bivariate regressions. All three metrics show consistent negative slopes (higher latitude → smaller deviations), consistent with compressed phenology reducing seasonal detection fluctuations. Species richness is **not significant** for any metric and adds negligible explanatory power when included alongside latitude.
- **Suppression effect was an artifact of spatial coverage confounds**: The previously reported suppression pattern (latitude non-significant alone but significant controlling for richness) was based on the old `d_spatial_cov` metric, which had denominator and accumulation confounds. Once these were resolved via `d_lambda` and `d_matched_rate`, the suppression mechanism disappears — latitude acts directly, and the portfolio effect at the dataset level is absent.

### Environmental covariate comparison (replacing latitude)

Five candidate environmental covariates were tested as replacements for raw latitude in the old brms model structure (15 models = 5 covariates × 3 metrics):

| Covariate | r with latitude | LOO rank (lambda / rate / mrate) |
|-----------|-----------------|----------------------------------|
| Temp seasonality (BIO4) | 0.69 | **#1 / #2 / #1** |
| NDVI amplitude | 0.36 | #2 / #1 / #2 |
| Latitude | 1.00 | #3 / #3 / #4 |
| Annual precip (BIO12) | −0.33 | #4 / #5 / #3 |
| Precip seasonality (BIO15) | −0.25 | #5 / #4 / #5 |

**Temperature seasonality (BIO4) is the best-performing environmental predictor** across metrics. This finding is confirmed in the GAM sensitivity surface models (Q5), where BIO4 remains significant for encounter rate and matched rate (β ≈ −0.21 to −0.23, p < 0.007), with a marginal effect for lambda (β = −0.16, p = 0.069).

**Note on `dataset_metadata.csv`**: BE-Leuven is missing from this file. Its environmental covariates were extracted manually: BIO4 = 559.2, BIO12 = 774, BIO15 = 11.2, NDVI amplitude = 0.449. The `centroid_lat` for BE-Leuven IS present in the file (50.80). The `join_env_covariates()` function in `prep_sensitivity_data.R` handles both the BE-Leuven manual fix and the slice-name matching for all covariates including `centroid_lat`.

### Slovenian sites drive mid-seasonality band deviations

The hump-shaped appearance in the BIO4 sensitivity surface (mid-seasonality band having highest deviations) is **not a systematic quadratic relationship** — it is driven by 5 Slovenian datasets (SI-eow, SI-senji, SI-abnik, SI-serknica × 2) which have both intermediate temperature seasonality (BIO4 ≈ 664–724) and anomalously high deviations. The two non-Slovenian datasets in the same band (GE-BFNP_201920, NO-nina) have low deviations consistent with the overall declining trend.

The Slovenian signal is dominated by **ungulate species with strong rutting-season detection spikes**: red deer (*Cervus elaphus*) shows d_lambda up to 0.10 in autumn windows at Slovenian sites vs near-zero elsewhere, and roe deer (*Capreolus capreolus*) shows a similar pattern. SI-eow has the highest mean |d_lambda| of any dataset (0.043), roughly 4× higher than non-Slovenian sites at similar BIO4 values.

## Known Issues and Bug Fixes

### `centroid_lat` was not propagated to sliced datasets in richness data (FIXED 2026-03-11)

`join_env_covariates()` in `prep_sensitivity_data.R` originally included `centroid_lat` only in the first-pass join (exact dataset name match) but not in the second-pass (slice-name matching). This caused all 13 sliced datasets (BE-Leuven ×5, NO-evenstadlia ×2, SI-serknica ×2, SP-donana ×4) to have NA `s_latitude` in the richness data. These 11,050 rows were silently dropped by `bam()`, resulting in community models fitted on ~9k rows instead of ~20k. The species-level data was unaffected because it uses `latitude` from the raw window data rather than `centroid_lat` from the metadata join.

**Fix**: `join_env_covariates()` now carries `centroid_lat` alongside the environmental covariates through all three join passes.

### `igraph::crossing()` masks `tidyr::crossing()` (FIXED 2026-03-11)

If `igraph` is loaded, its `crossing()` function masks `tidyr::crossing()`, causing Q6 protocol evaluation to fail. `sensitivity_results.R` now uses `tidyr::crossing()` explicitly.

### rho_lambda model boundary inflation (FIXED 2026-03-12)

The original rho_lambda community model included all non-NA observations, but ~33% of rho values were exactly 1 — artefactual perfect correlations from windows sharing only 3–4 species (Spearman's rho with 3 species has very few discrete possible values). This boundary mass violated beta regression assumptions. Filtering to windows with ≥5 shared species (currently N = 12,028 after dataset expansion) removed most artefactual boundary inflation (dev. explained = 40.4%). See `methods_note_rho_filter.md` for manuscript text.

### `abs_d_rate` has 1 exact zero requiring offset (NOTED 2026-03-18)

With the extended 183-day window range, 1 observation out of 222,748 has `abs_d_rate == 0`. A `pmax(..., 1e-10)` offset is applied before Gamma regression. This was not needed with the original 120-day range.

### SP-donana deployment end-dates corrected (FIXED 2026-03-17)

The original `deployments.csv` for SP-donana had ~50 deployment end-dates that extended beyond actual camera operation. Replaced with `deployments_new.csv` which trims these end-dates. Old file backed up as `deployments_backup_20260316.csv`. The correction changed metric values throughout (all d_rate values, ~91% of d_lambda, ~68% of d_matched_rate) but did not change which species pass thresholds (same row counts per slice). Mean absolute deviations increased slightly (e.g., |d_lambda| from 0.0134 → 0.0145 for SP-donana). Full pipeline re-run, prep data regenerated, and all models re-fitted.

## Parallelism Settings
- `N_THREADS = 4` for `bam()` fitting
- Old brms models used `N_WORKERS = 4` × `CORES_PER_MODEL = 4` = 16 cores

## Superseded Scripts (kept for reference)
- `models_abs_UPDATED.R` → replaced by `models_abs_detection_PARALLEL.R` + `models_abs_richness_PARALLEL.R`
- `models_mse_UPDATED.R`, `models_mse.R` → replaced by `models_mse_detection.R` + `models_mse_richness.R`
- `models_abs.R`, `modeling.R`, `parallel.R`, `parallel_after.R` → older versions
- `Full.R`, `Full0.R` → replaced by `Full1.R`
- `Full_visu0.R` → replaced by `Full_visu.R`
- `dataset_overview.R`, `dataset_overview_v2.R` → replaced by `dataset_report.R`
- `dataset-slices_overview.R` → replaced by `dataset_report.R`

## `set_camtrap_cols()` usage pattern

As of the v2 refactor, `set_camtrap_cols()` returns a named list instead of
injecting into the caller's frame via `list2env`. All active callers unpack it
explicitly:

```r
cols <- set_camtrap_cols(obs, deploy)
list2env(cols, envir = environment())
# sci_col, date_col, depl_id_obs, depl_id_deploy,
# loc_col, lat_col, lon_col, start_col, end_col
# are now available as local variables
```

Superseded scripts (Full.R, Full0.R, parallel.R, parallel_after.R,
diag_12m_slice.R, diag_12m_slice2.R) still use the old `set_camtrap_cols(obs, deploy)`
calling style and have not been updated.

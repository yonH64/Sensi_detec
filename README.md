# Sensitivity of Camera-Trap Sampling Window Design on Detection Metrics

How does the choice of temporal sampling window affect what camera traps tell us about wildlife? This project builds a **sensitivity surface** — a systematic map of how detection metric deviation from a 12-month benchmark varies as a function of window **timing** and **duration** — across 26 European camera-trap datasets (35 dataset-slices) and 29 mammal species. A benchmark noise floor analysis (inter-annual variability of the 12-month reference) provides a species-specific stopping rule for the diminishing-returns question: *how long is long enough?*

![Dataset locations](figures/dataset_map.png)

------------------------------------------------------------------------

## Research Question

Standardised camera-trap protocols like Snapshot Europe prescribe fixed temporal windows (e.g., 61 days starting September 1). But how sensitive are the resulting detection estimates to this specific choice? Rather than comparing two protocols, we ask: **how well does this short window approximate what year-round monitoring would give me?**

We construct a sliding window grid (25 durations [15–183 days] × 53 start positions, stepped weekly across the year) and model the resulting \~223,000 species × window × dataset observations using species-specific GAM surfaces.

------------------------------------------------------------------------

## Key Results

### The sensitivity surface

The main result: a 2D heatmap of predicted absolute deviation in daily detection rate (\|Δλ\|), averaged across species. Bright colours = high deviation. Cyan diamonds mark the Snapshot Europe CORE and BUFFER windows, plus the ENETWILD EOW split.

![Main sensitivity surface](figures/fig1_sensitivity_surface.png)

**Short windows centred on autumn show the largest deviations.** The signal is driven by ungulate rutting-season activity spikes — windows capturing the rut inflate detection rates well above the annual average.

### Species-level variation

Species identity dominates the surface (ΔAIC ≈ −19,200 over guild-level models). Six focal species illustrate the range:

![Species-specific surfaces](figures/fig3_species_surfaces.png)

*Cervus elaphus* and *Sus scrofa* show intense autumn hotspots (rut), while *Vulpes vulpes* and *Lepus europaeus* have near-flat surfaces — their detection is relatively stable year-round.

### Duration effect

Deviation drops steeply with window length, then plateaus. \~80% of the reduction occurs by 60 days. The seasonal sensitivity is strongest for short windows.

![Duration curves by season](figures/fig4_duration_curves.png)

### Signed encounter rate deviation

Not all bias goes one direction. Summer/autumn windows tend to **overestimate** encounter rates (red), while winter/spring windows **underestimate** them (blue):

![Signed rate surface](figures/fig5_signed_rate_surface.png)

### Community-level recovery

Species richness recovery (proportion of full-year species detected) favours spring-started windows and increases monotonically with duration. Rank preservation of detection rates (Spearman ρ) is uniformly high (0.91–0.98) — even biased windows preserve species ordering.

![Richness surfaces](figures/fig6_richness_surfaces.png)

### Protocol evaluation

The Snapshot Europe CORE and BUFFER windows sit on the edge of the autumn deviation hotspot. BUFFER benefits primarily from its longer duration (89d vs 61d), not its timing.

![Protocol evaluation](figures/fig7_protocol_evaluation.png)

### Guild-level surfaces

Deviation surfaces broken down by major taxonomic guild — ungulates dominate the autumn signal:

![Guild surfaces](figures/fig2_guild_surfaces.png)

### Benchmark noise floor & signal-to-noise ratio

The 12-month benchmark itself fluctuates year-to-year (median inter-annual CV = 23% across 33 species × site combinations from 4 multi-year sites). A sub-window's deviation is only informative if it exceeds this inter-annual noise. We define **SNR = |Δλ| / SD(λ_full)** — when SNR < 1, the deviation is smaller than what natural year-to-year variation produces.

| Window duration | Median SNR | % observations below noise floor |
|-----------------|-----------|----------------------------------|
| 15 d            | 86        | 0%                               |
| 60 d            | ~22       | 0.1%                             |
| 90 d            | ~13       | 0.6%                             |
| 120 d           | 8.1       | 2.3%                             |

For most species, the noise floor is not reached even at 120 days — the sensitivity surface carries real signal well beyond typical protocol durations. The stopping rule is genuinely **species-specific**: a handful of high-CV species (e.g., *Alces alces*, *Sciurus vulgaris* at specific sites) approach SNR < 1 at long windows, while most species remain well above.

### Benchmark robustness

The sensitivity surface shape is **robust to benchmark choice**. Recomputing deviations against 180-day and 270-day centred benchmarks preserves the seasonal profile (shape correlation r > 0.88 for all window lengths vs 180d; r > 0.96 for short windows vs 270d). The autumn rut spike and the duration smoothing curve appear regardless of benchmark duration. See supplementary Fig S1 (`figures/FigS1_robustness_benchmark.pdf`).

### Robustness summary

We tested the sensitivity of the main results to every analytical choice point — data pipeline parameters, species-inclusion thresholds, and model specification — by varying each in turn and comparing the predicted surface shape against the baseline. The table below summarises all robustness checks.

| Check | What was varied | Configurations tested | Surface *r* vs baseline | Dev. expl. range | Conclusion |
|-------|-----------------|----------------------|-------------------------|-----------------|------------|
| **Species thresholds** (individual) | min_events, min_sites, min_occasions | 7 variants across 3 thresholds | 0.925–1.000 | 85.8–87.0% | Robust; occasions threshold non-binding |
| **Species thresholds** (joint) | All 3 thresholds simultaneously | Lenient (10/3/3), Baseline (20/5/5), Strict (30/10/5) | 0.970–1.000 | 87.1–87.2% | Robust; no interaction effects |
| **Window-start resolution** | Step size between start positions | 3-day, 7-day, 14-day | 0.979–1.000 | 87.1% (all) | 7-day step adequate |
| **Anchor detection** | Annual slice detection parameters (5 params) | Relaxed, Current, Strict | 0.980–1.000 | 87.0–87.1% | Relaxed = Current; strict drops 3 slices, minimal effect |
| **Independence gap** | Minimum time between independent events | 15-min, 30-min, 60-min | 0.9998–1.000 | 87.1% (all) | Functionally identical |
| **Basis dimension** | Tensor product k values | k = (8,6), k = (16,12) | 0.996 | 87.1–87.3% | Baseline adequate |
| **Response family** | Gamma vs Tweedie | Gamma, Tweedie (estimated *p* = 1.99) | 0.997 | — | Gamma validated |
| **Random effects** | Flat vs nested (base-dataset + slice) | 254 levels vs 183 + 254 levels | 0.985 | — | Slice-level RE sufficient |
| **BIO4 nonlinearity** | Linear vs smooth | Linear, s(bio4, k = 5) → edf = 1.00 | — | — | No nonlinear evidence |
| **Rho filter cutoff** | Min shared species for rank correlation model | ≥3, ≥4, ≥5, ≥6, ≥8 | 0.911–0.953 | 19.0–55.3% | ≥5 balances boundary mass vs coverage |
| **Benchmark duration** | Reference period for deviation computation | 180d, 270d, 365d | 0.879–0.998 (shape *r*) | — | Surface shape preserved |

Across all checks, predicted surface correlations with the baseline exceed 0.92, and in most cases exceed 0.97. No analytical choice point changes the qualitative conclusions. Full details for each check are in the corresponding appendix documents.

------------------------------------------------------------------------

## Methods Overview

| Stage | Script | Description |
|------------------|---------------------|---------------------------------|
| **Data pipeline** | `Full1.R` | Processes raw CamTrap DP archives → detection metrics per species × window × dataset |
| **Data preparation** | `prep_sensitivity_data.R` | Joins traits, environmental covariates, standardises for modelling |
| **Model fitting** | `models_sensitivity_surface.R` | GAMMs via `mgcv::bam()` with cyclic splines and species-specific tensor product surfaces |
| **Results extraction** | `sensitivity_results.R` | Derives Q1–Q9 outputs (duration effect, seasonal profiles, optimal timing, protocol evaluation, benchmark noise floor, etc.) |
| **Figures** | `sensitivity_figures.R` | 9 main + 1 supplementary publication figures |

The best model (M6) fits species-specific 2D surfaces over day-of-year (cyclic) × duration, with temperature seasonality (BIO4), effort controls, and dataset × species random effects. Gamma(log) family for absolute deviations. Deviance explained: **87.1%** for the primary metric (\|Δλ\|).

Full analytical details: [`PAGER_ANALYTICAL_WORKFLOW.md`](PAGER_ANALYTICAL_WORKFLOW.md)

------------------------------------------------------------------------

## Supplementary & Appendix Documents

### Robustness & sensitivity checks

| Document | Description |
|----------|-------------|
| [`appendix_threshold_robustness.md`](appendix_threshold_robustness.md) | Sensitivity to individual species-inclusion thresholds (`min_events`, `min_sites_pos`, `min_occasions_pos`) |
| [`appendix_joint_threshold_sensitivity.md`](appendix_joint_threshold_sensitivity.md) | Sensitivity to all three species-inclusion thresholds varied simultaneously |
| [`appendix_resolution_sensitivity.md`](appendix_resolution_sensitivity.md) | Sensitivity to window-start step size (3-day, 7-day, 14-day) |
| [`appendix_anchor_sensitivity.md`](appendix_anchor_sensitivity.md) | Sensitivity to annual slice detection parameters (`find_anchors()`) |
| [`appendix_independence_threshold_sensitivity.md`](appendix_independence_threshold_sensitivity.md) | Sensitivity to the independence gap for defining detection events (15, 30, 60 min) |
| [`appendix_model_diagnostics.md`](appendix_model_diagnostics.md) | Model specification diagnostics: basis adequacy (k-doubling), Gamma vs Tweedie, nested random effects, nonlinear BIO4 |
| [`appendix_rho_filter_sensitivity.md`](appendix_rho_filter_sensitivity.md) | Sensitivity of the rank correlation (ρ) model to the minimum shared species filter (≥3 to ≥8) |
| [`appendix_benchmark_robustness.md`](appendix_benchmark_robustness.md) | Sensitivity of the surface shape to benchmark duration (180d, 270d, 365d) |

### Methods & analytical notes

| Document | Description |
|----------|-------------|
| [`appendix_structural_confounds.md`](appendix_structural_confounds.md) | Three structural confounds identified in detection metrics and their resolutions |
| [`appendix_model_comparison.md`](appendix_model_comparison.md) | Model selection: M1–M6 + guild variants (AIC comparison, surface structure rationale) |
| [`appendix_rho_filter.md`](appendix_rho_filter.md) | Spearman's ρ boundary inflation and the ≥5 shared species filter |
| [`methods_note_rho_filter.md`](methods_note_rho_filter.md) | Concise methods text for the ρ filter (for manuscript insertion) |
| [`discussion_benchmark_and_framing.md`](discussion_benchmark_and_framing.md) | Discussion: benchmark interpretation and alternative research framings |
| [`PAGER_ANALYTICAL_WORKFLOW.md`](PAGER_ANALYTICAL_WORKFLOW.md) | Full analytical workflow documentation |

------------------------------------------------------------------------

## Key Scripts

| Script | Role |
|-----------------------------------------|-------------------------------|
| `helpers.R` | Shared utility functions (sourced by all other scripts) |
| `Full1.R` | Main data pipeline orchestration |
| `prep_sensitivity_data.R` | Model data preparation |
| `models_sensitivity_surface.R` | All GAM model fitting |
| `sensitivity_results.R` | Research question outputs (Q1–Q9 CSVs) |
| `sensitivity_figures.R` | Publication figures (Figs 1–9, Fig S1) |
| `dataset_report.R` | Per-dataset diagnostic reports + environmental covariate extraction |

------------------------------------------------------------------------

## Data

Raw camera-trap datasets are not included in this repository. Model outputs (`.rds`, `.RData`) and generated results (`Q*.csv`, `Fig*.pdf`) are excluded via `.gitignore` but can be regenerated by running the pipeline.

**Species**: 29 species passing data thresholds across ≥3 datasets (of 63 in the full taxonomy). **Datasets**: 26 camera-trap arrays across 10 European countries, yielding 35 dataset-slices after temporal splitting. **Observations**: 222,748 species × window × dataset rows for species-level models; 46,375 for community-level models.

------------------------------------------------------------------------

*Last updated: 19 March 2026*

# Sensitivity of Camera-Trap Sampling Window Design on Detection Metrics

How does the choice of temporal sampling window affect what camera traps tell us about wildlife? This project builds a **sensitivity surface** — a systematic map of how detection metric deviation from a 12-month benchmark varies as a function of window **timing** and **duration** — across 26 European camera-trap datasets (35 dataset-slices) and 26 mammal species. A benchmark noise floor analysis (inter-annual variability of the 12-month reference) provides a species-specific stopping rule for the diminishing-returns question: *how long is long enough?*

![Dataset locations](figures/dataset_map.png)

------------------------------------------------------------------------

## Research Question

Standardised camera-trap protocols like Snapshot Europe prescribe fixed temporal windows (e.g., 61 days starting September 1). But how sensitive are the resulting detection estimates to this specific choice? Rather than comparing two protocols, we ask: **how well does this short window approximate what year-round monitoring would give me?**

We construct a sliding window grid (16 durations × 53 start positions, stepped weekly across the year) and model the resulting \~124,000 species × window × dataset observations using species-specific GAM surfaces.

------------------------------------------------------------------------

## Key Results

### The sensitivity surface

The main result: a 2D heatmap of predicted absolute deviation in daily detection rate (\|Δλ\|), averaged across species. Bright colours = high deviation. Cyan diamonds mark the Snapshot Europe CORE and BUFFER windows, plus the ENETWILD EOW split.

![Main sensitivity surface](figures/fig1_sensitivity_surface.png)

**Short windows centred on autumn show the largest deviations.** The signal is driven by ungulate rutting-season activity spikes — windows capturing the rut inflate detection rates well above the annual average.

### Species-level variation

Species identity dominates the surface (ΔAIC ≈ −10,600 over guild-level models). Six focal species illustrate the range:

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

------------------------------------------------------------------------

## Methods Overview

| Stage | Script | Description |
|------------------|---------------------|---------------------------------|
| **Data pipeline** | `Full1.R` | Processes raw CamTrap DP archives → detection metrics per species × window × dataset |
| **Data preparation** | `prep_sensitivity_data.R` | Joins traits, environmental covariates, standardises for modelling |
| **Model fitting** | `models_sensitivity_surface.R` | GAMMs via `mgcv::bam()` with cyclic splines and species-specific tensor product surfaces |
| **Results extraction** | `sensitivity_results.R` | Derives Q1–Q9 outputs (duration effect, seasonal profiles, optimal timing, protocol evaluation, benchmark noise floor, etc.) |
| **Figures** | `sensitivity_figures.R` | 9 main + 1 supplementary publication figures |

The best model (M6) fits species-specific 2D surfaces over day-of-year (cyclic) × duration, with temperature seasonality (BIO4), effort controls, and dataset × species random effects. Gamma(log) family for absolute deviations. Deviance explained: **86.8%** for the primary metric (\|Δλ\|).

Full analytical details: [`PAGER_ANALYTICAL_WORKFLOW.md`](PAGER_ANALYTICAL_WORKFLOW.md)

------------------------------------------------------------------------

## Supplementary & Appendix Documents

| Document | Description |
|----------|-------------|
| [`appendix_threshold_robustness.md`](appendix_threshold_robustness.md) | Sensitivity of results to species-inclusion thresholds (`min_events`, `min_sites_pos`, `min_occasions_pos`) |
| [`appendix_benchmark_robustness.md`](appendix_benchmark_robustness.md) | Sensitivity of the surface shape to benchmark duration (180d, 270d, 365d) |
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

**Species**: 26 species passing data thresholds across ≥3 datasets (of 63 in the full taxonomy). **Datasets**: 26 camera-trap arrays across 10 European countries, yielding 35 dataset-slices after temporal splitting. **Observations**: 123,902 species × window × dataset rows for species-level models; \~30,000 for community-level models.

------------------------------------------------------------------------

*Last updated: 18 March 2026*

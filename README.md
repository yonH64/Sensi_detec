# Sensitivity of Camera-Trap Sampling Window Design on Detection Metrics

How does the choice of temporal sampling window affect what camera traps tell us about wildlife? This project builds a **sensitivity surface** — a systematic map of how detection metric deviation from a 12-month benchmark varies as a function of window **timing** and **duration** — across 23 European camera-trap datasets and 24 mammal species.

![Dataset locations](figures/dataset_map.png)

------------------------------------------------------------------------

## Research Question

Standardised camera-trap protocols like Snapshot Europe prescribe fixed temporal windows (e.g., 61 days starting September 1). But how sensitive are the resulting detection estimates to this specific choice? Rather than comparing two protocols, we ask: **how well does this short window approximate what year-round monitoring would give me?**

We construct a sliding window grid (16 durations × 53 start positions, stepped weekly across the year) and model the resulting \~86,000 species × window × dataset observations using species-specific GAM surfaces.

------------------------------------------------------------------------

## Key Results

### The sensitivity surface

The main result: a 2D heatmap of predicted absolute deviation in daily detection rate (\|Δλ\|), averaged across species. Bright colours = high deviation. Cyan diamonds mark the Snapshot Europe CORE and BUFFER windows, plus the ENETWILD EOW split.

![Main sensitivity surface](figures/fig1_sensitivity_surface.png)

**Short windows centred on autumn show the largest deviations.** The signal is driven by ungulate rutting-season activity spikes — windows capturing the rut inflate detection rates well above the annual average.

### Species-level variation

Species identity dominates the surface (ΔAIC ≈ −9,000 over guild-level models). Six focal species illustrate the range:

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

------------------------------------------------------------------------

## Methods Overview

| Stage | Script | Description |
|------------------|---------------------|---------------------------------|
| **Data pipeline** | `Full1.R` | Processes raw CamTrap DP archives → detection metrics per species × window × dataset |
| **Data preparation** | `prep_sensitivity_data.R` | Joins traits, environmental covariates, standardises for modelling |
| **Model fitting** | `models_sensitivity_surface.R` | GAMMs via `mgcv::bam()` with cyclic splines and species-specific tensor product surfaces |
| **Results extraction** | `sensitivity_results.R` | Derives Q1–Q8 outputs (duration effect, seasonal profiles, optimal timing, protocol evaluation, etc.) |
| **Figures** | `sensitivity_figures.R` | 8 publication figures |

The best model (M6) fits species-specific 2D surfaces over day-of-year (cyclic) × duration, with temperature seasonality (BIO4), effort controls, and dataset × species random effects. Gamma(log) family for absolute deviations. Deviance explained: **89.8%** for the primary metric (\|Δλ\|).

Full analytical details: [`PAGER_ANALYTICAL_WORKFLOW.md`](PAGER_ANALYTICAL_WORKFLOW.md)

------------------------------------------------------------------------

## Key Scripts

| Script | Role |
|-----------------------------------------|-------------------------------|
| `helpers.R` | Shared utility functions (sourced by all other scripts) |
| `Full1.R` | Main data pipeline orchestration |
| `prep_sensitivity_data.R` | Model data preparation |
| `models_sensitivity_surface.R` | All GAM model fitting |
| `sensitivity_results.R` | Research question outputs (Q1–Q8 CSVs) |
| `sensitivity_figures.R` | Publication figures |
| `dataset_report.R` | Per-dataset diagnostic reports + environmental covariate extraction |

------------------------------------------------------------------------

## Data

Raw camera-trap datasets are not included in this repository. Model outputs (`.rds`, `.RData`) and generated results (`Q*.csv`, `Fig*.pdf`) are excluded via `.gitignore` but can be regenerated by running the pipeline.

**Species**: 24 species passing data thresholds across ≥3 datasets (of 63 in the full taxonomy). **Datasets**: 23 camera-trap arrays across 10 European countries, yielding 24 dataset-slices after temporal splitting. **Observations**: 86,108 species × window × dataset rows for species-level models; \~20,000 for community-level models.

------------------------------------------------------------------------

*Last updated: March 2026*

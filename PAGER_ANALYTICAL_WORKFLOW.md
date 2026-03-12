---
output:
  pdf_document: default
  html_document: default
---
# Analytical Workflow — Sensitivity of Detection Protocol Study

## Research Question

How does the choice of camera-trap temporal sampling window affect species-level detection metrics and community-level richness estimates? We quantify this using a **sensitivity surface approach**: a systematic grid of sliding windows (varying in duration and timing) is compared against a full 12-month reference period across multiple European camera-trap datasets. The Snapshot Europe **CORE** (61 days, Sep 1–Oct 31) and **BUFFER** (89 days, Aug 18–Nov 14) protocols are evaluated as specific points on this surface, not as the primary comparison.

---

## 1. Data Pipeline (`Full1.R`)

### 1.1 Input Data

- Multiple camera-trap datasets stored as CamTrap DP archives (deployments.csv + observations.csv).
- Each dataset may span multiple years. Multi-year datasets are split into non-overlapping 12-month "anchor" periods (slices) using an automated scoring algorithm that maximises trapping effort and detection coverage. Multiple slices from the same dataset are treated as independent dataset-slices.
- Observations are filtered to a curated list of ~61 European terrestrial mammals (lagomorphs, insectivores, carnivores, ungulates, rodents, and one primate).

### 1.2 Window Grid

A systematic grid of ~850 sampling windows is constructed, varying by:
- **Start day-of-year** (every 7 days across the calendar year, 53 positions)
- **Window length** (15 to 120 days in 7-day steps, 16 durations)

Four named protocol windows are added:
| Window | Start | End | Length | Description |
|--------|-------|-----|--------|-------------|
| SNAP_EU_CORE | Sep 1 (DOY 244) | Oct 31 (DOY 304) | 61 d | Snapshot Europe core window |
| SNAP_EU_BUFFER | Aug 18 (DOY 230) | Nov 14 (DOY 318) | 89 d | Snapshot Europe extended window |
| EOW_EARLY | Aug 2 (DOY 214) | Sep 30 (DOY 273) | 60 d | EOW/ENETWILD pre-rut half |
| EOW_LATE | Oct 1 (DOY 274) | Nov 29 (DOY 333) | 60 d | EOW/ENETWILD post-rut half |

A **FULL** window (365 days) serves as the benchmark for all deviation metrics.

### 1.3 Metric Computation

For each dataset-slice, for each window in the grid, and for each species that passes minimum-data thresholds, the pipeline computes:

**Species-level detection metrics** (computed per species × window × dataset):
- **Daily detection rate** (`lambda`): Derived from a time-to-first-detection exponential model across cameras. Estimated λ (hazard rate) from first-detection times, with censoring for cameras where the species was not detected. Effort-free and window-length-independent.
- **Event rate** (`rate`): Total independent events / total trap-days. Independence threshold: 30-minute gap per camera × species.
- **Spatial detection coverage** (`spatial_cov`): Proportion of cameras with ≥1 detection. Retained as a descriptive statistic; not used for deviations directly (structurally confounded with window length).
- **Matched-camera detection rate** (`matched_rate`): `-log(1 - spatial_cov_m) / L`, where `spatial_cov_m` is spatial coverage computed only on cameras active in both the sub-window and the FULL period. Fully deconfounded for both denominator and accumulation effects.
- **Standard errors** for each metric (delta method for TTE; camera-level empirical SE for rate; derived SE for matched rate).

**Community-level richness metrics** (computed per window × dataset):
- **Observed richness** (`sr_obs`): Count of species passing thresholds.
- **Rarefied richness** (`sr_raref`): Incidence-based rarefaction (Colwell et al. 2012) to a common number of sampled camera × occasion units across windows. Analytical variance is computed.
- **Proportion of FULL species detected** (`prop_sr_full`): |A ∩ FULL| / |FULL|.
- **Rank correlation** (`rho_lambda`): Spearman rank correlation of species detection rates between window and FULL.

**Deviation metrics** (deltas to the FULL benchmark):
- For each species-level metric: `d_X = X_window − X_FULL` and `abs_d_X = |d_X|`.
- MSE criterion: `mse_X = (X_window − X_FULL)² + SE²_window` (bias² + variance).
- For richness: `d_sr_raref`, `mse_sr_raref = d_sr_raref² + sr_raref_var`.

### 1.4 Minimum-Data Thresholds

A species × window combination is retained only if:
- ≥ 20 independent events
- ≥ 5 positive occasions (days with detections)
- ≥ 5 camera sites with detections

Species × window combinations failing these thresholds are logged separately.

### 1.5 Species Traits

Species are classified into functional guilds via a separate taxonomy file (61 species):

| Guild level | Categories |
|---|---|
| `guild_major` (6) | Carnivore, Insectivore, Lagomorph, Primate, Rodent, Ungulate |
| `guild_minor_habitat` (11) | e.g., Large predator, Meso predator, Small predator, Forest ungulate, Open/mountain ungulate, Domestic ungulate, Aquatic rodent, ... |
| `guild_minor_diet` (12) | e.g., Large predator, Meso predator, Browser ungulate, Grazer ungulate, Mixed feeder ungulate, ... |

### 1.6 Environmental Covariates

Per-dataset environmental covariates are extracted from WorldClim v2.1 and MODIS:
- **Temperature seasonality** (`bio4_temp_seasonality`): SD of monthly temperature × 100. Best environmental predictor of deviation magnitude.
- **Annual precipitation** (`bio12_annual_precip`), **Precipitation seasonality** (`bio15_precip_seasonality`).
- **NDVI amplitude** (`ndvi_amplitude`): Annual max − min NDVI.

---

## 2. Data Preparation (`prep_sensitivity_data.R`)

The preparation script transforms the raw window metric outputs into model-ready data frames:

### 2.1 Species-Level Data (`sensitivity_species_data.rds`)
- **86,108 rows** (24 species × ~848 windows × 23 datasets, minus unavailable combinations)
- Extracts `day_start` from window identifiers; computes `day_center = (day_start + window_len/2) % 365`
- Computes absolute deviations: `abs_d_lambda`, `abs_d_rate`, `abs_d_matched_rate`
- Joins species traits (`guild_major`, `guild_minor_habitat`, `guild_minor_diet`)
- Joins environmental covariates (with special handling for sliced datasets and BE-Leuven)
- Standardises continuous covariates: `s_bio4`, `s_latitude`, `s_trap_array`
- Creates grouping factors: `species_f`, `dataset_f`, `ds_sp_f`

### 2.2 Richness-Level Data (`sensitivity_richness_data.rds`)
- **20,400 rows** (windows × 24 datasets including sliced datasets)
- Same covariate engineering as species-level

---

## 3. Statistical Models (`models_sensitivity_surface.R`)

All models are GAMMs fitted with `mgcv::bam()` using `method = "fREML"`, `discrete = TRUE`, and `nthreads = 4`. Outputs `sensitivity_gam_models.rds` and `sensitivity_models_env.RData` (cleaned environment loaded by downstream scripts). Total fitting time: ~2–3 minutes on 4 threads.

### 3.1 Model Comparison (abs_d_lambda only)

Eight model variants are compared to determine the optimal structure for the sensitivity surface:

| Model | Structure | AIC | Dev. expl. |
|-------|-----------|-----|------------|
| **M6** | Species-specific 2D surfaces (24 species) | **−810,440** | **89.8%** |
| M_diet | Diet-guild-specific surfaces (10 levels) | −801,378 | 88.7% |
| M5 | Major guild surfaces + bio4 × duration | −798,311 | 88.3% |
| M_hab | Habitat-guild-specific surfaces (9 levels) | −798,663 | 88.3% |
| M3 | Major guild-specific surfaces (5 levels) | −796,949 | 88.1% |
| M2 | Shared surface + guild-varying seasonality | −796,076 | 88.0% |
| M4 | Shared surface + bio4 × duration interaction | −792,754 | 87.5% |
| M1 | Single shared surface | −791,578 | 87.3% |

**Species-specific surfaces (M6) dominate**, with ΔAIC ≈ −9,063 over the best guild-level model (M_diet). The penalisation automatically regularises rare species based on data volume.

### 3.2 Best Model Formula (M6)

```r
abs_d_metric ~
  te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6), by = species_f) +
  species_f +
  s_bio4 + l_trapdays + l_nsites + s_latitude + s_trap_array +
  s(ds_sp_f, bs = "re")
```

**Key design choices:**
- **Cyclic spline** (`bs = "cc"`) for `day_start` with `knots = list(day_start = c(0, 365))`: ensures the seasonal pattern wraps smoothly from December to January. ~19% of windows cross the year boundary.
- **Tensor product** (`te()`): allows the timing × duration interaction to be non-separable — short windows are more sensitive to timing choice than long windows.
- **`by = species_f`**: separate surface per species. Penalisation controls complexity adaptively based on data volume.
- **`Gamma(link = "log")`**: appropriate for strictly positive, right-skewed absolute deviations.
- **Random effects** `s(ds_sp_f, bs = "re")`: equivalent to `(1 | dataset:species)` in mixed model notation.

### 3.3 Detection Models (M6 structure, all metrics)

| Response | Family | N | Dev. explained |
|----------|--------|---|----------------|
| `abs_d_lambda` (TTE daily detection rate) | Gamma(log) | 86,108 | 89.8% |
| `abs_d_matched_rate` (deconfounded spatial coverage rate) | Gamma(log) | 86,108 | 88.4% |
| `abs_d_rate` (encounter rate) | Gamma(log) | 86,108 | 67.4% |
| `d_rate` (signed encounter rate deviation) | Gaussian | 86,108 | 18.5% |

### 3.4 Community Models (shared surface)

| Response | Family | N | Dev. explained |
|----------|--------|---|----------------|
| `prop_sr_full` (proportion of species recovered) | Beta(logit) | 19,505 | 70.6% |
| `d_sr_raref` (rarefied richness deviation) | Gaussian | 20,353 | 27.2% |
| `rho_lambda` (rank correlation of detection rates) | Beta(logit) | 9,887 | 37.1% |

Community models use a single shared surface (no species dimension), with `s(dataset_f, bs = "re")` for dataset-level random effects. The `rho_lambda` model is restricted to windows with ≥5 shared species to avoid artefactual boundary inflation from Spearman's rho with very few species (see `methods_note_rho_filter.md`).

---

## 4. Research Questions (`sensitivity_results.R`)

| Q# | Question | Approach | Key output |
|----|----------|----------|------------|
| Q1 | How does deviation depend on window duration? | Marginal effect of `window_len`, aggregated from species surfaces | Duration-deviation curves (mean, IQR) per guild |
| Q2 | How does deviation depend on timing (season)? | Marginal effect of `day_start` at fixed durations | Seasonal profiles showing autumn peak |
| Q3 | How do duration and timing interact? | Full 2D tensor product surface predictions | Heatmap: timing × duration, per metric |
| Q4 | How do species traits modulate the surface? | Species-specific surfaces aggregated by guild | Season × species/guild deviation tables |
| Q5 | How does temperature seasonality modulate deviation? | Parametric `s_bio4` coefficient | Effect size per metric (all significantly negative) |
| Q6 | Where do named protocols sit on the surface? | Predict at CORE (day 244, 61d), BUFFER (day 230, 89d), EOW_EARLY (day 214, 60d), EOW_LATE (day 274, 60d) | Predicted deviation vs surface optimum |
| Q7 | What is the optimal window design? | Minimum-deviation timing per species/guild at key durations | Optimal start day per species at 30/60/90d |
| Q8 | How does richness recovery vary? | Community model predictions (prop_sr_full, d_sr_raref, rho) | Richness surface: spring-summer peak for recovery |

---

## 5. Figures (`sensitivity_figures.R`)

| Figure | Content |
|--------|---------|
| Fig 1 | Main sensitivity surface (all species pooled) with named protocol positions (CORE, BUFFER, EOW_EARLY, EOW_LATE) |
| Fig 2 | Guild-specific sensitivity surfaces (4 panels) |
| Fig 3 | Species-specific surfaces (6 focal species) |
| Fig 4 | Duration-deviation curves by season, with SE protocol durations marked |
| Fig 5 | Signed encounter rate deviation surface (bidirectional: red = overestimate, blue = underestimate) |
| Fig 6 | Community-level surfaces: species recovery, rank preservation, richness deviation |
| Fig 7 | Protocol evaluation: predicted deviation per guild × protocol (CORE, BUFFER, EOW_EARLY, EOW_LATE) × metric |
| Fig 8 | Model comparison AIC ladder |

---

## 6. Key Analytical Design Choices

1. **Sensitivity surface framing**: Instead of comparing two specific protocols, we model detection metric deviation as a continuous function of window timing and duration. This uses all ~86,000 species-level observations (vs ~170 in the old CORE-vs-BUFFER comparison) and generalises to any sampling window design.

2. **Benchmark**: The FULL (365-day) window is treated as a benchmark, not as ground truth. All deviations are relative to what a full year of sampling yields.

3. **Cyclic day-of-year**: The seasonal dimension uses cyclic cubic splines to ensure smooth wrapping from December to January.

4. **Species-specific surfaces**: Each species gets its own 2D sensitivity surface, with penalisation controlling complexity. This outperforms guild-level grouping by ΔAIC ≈ 9,000, reflecting genuine species-level differences in seasonal activity patterns.

5. **Structural confound resolution**: Raw spatial coverage (`spatial_cov`) has two confounds (denominator and accumulation effects). The matched-camera approach restricts computations to cameras active in both periods, and the `-log(1-scov)/L` transformation converts to a daily detection rate that is fully deconfounded. All models use `d_matched_rate` rather than the raw `d_spatial_cov`.

6. **Environmental predictor**: Temperature seasonality (BIO4) outperforms raw latitude as a predictor of deviation magnitude and has a clear mechanistic interpretation (sites with more pronounced seasons have more variable detection rates).

7. **Rarefaction for richness**: Incidence-based rarefaction (Colwell et al. 2012) standardises richness across windows with different sampling intensities. The sampling unit is the camera × occasion cell.

8. **Derived protocol evaluation**: All four named protocols (CORE, BUFFER, EOW_EARLY, EOW_LATE) are evaluated as specific coordinates on the fitted surface (Q6), allowing direct comparison with the surface optimum and quantification of the deviation "cost" of autumn timing. EOW_EARLY and EOW_LATE split the autumn at Oct 1 (the CORE midpoint), providing a pre-rut vs post-rut contrast at identical duration (60d).

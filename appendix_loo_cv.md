# Appendix: Leave-One-Dataset-Out Cross-Validation

The M6 model explains 87.1% of deviance in-sample, but this includes a dataset × species random intercept (SD = 0.94 on the log link) that absorbs site-level variation. To quantify how well the model's fixed-effect structure (species-specific surfaces + covariates) generalises to an unseen site, we performed leave-one-dataset-out cross-validation (LOO-CV) for all three absolute deviation metrics.

## Design

For each of the 27 base datasets, we re-fitted the M6 model on all remaining data, then predicted on the held-out observations with the random effect excluded. Base datasets rather than individual slices were held out, since slices from the same physical location share environmental covariates and would leak information.

Eleven species that appear in only one base dataset (*Castor fiber*, *Erinaceus europaeus*, *Genetta genetta*, *Gulo gulo*, *Herpestes ichneumon*, *Lepus granatensis*, *Lynx lynx*, *Martes foina*, *Myocastor coypus*, *Oryctolagus cuniculus*, *Ovis aries*) were necessarily dropped from the test set when their sole dataset was held out. Two additional species (*Procyon lotor*, *Ursus arctos*) appear in only two base datasets, so their training signal comes from a single site in each LOO iteration. In total, 207,047 of 218,663 observations received out-of-sample predictions (94.7%).

Each LOO iteration re-fits the full M6 model formula:

```
abs_d_metric ~ te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6), by = species_f) +
  species_f + l_trapdays + l_nsites + s_latitude + s_trap_array + s(ds_sp_f, bs = "re")
family = Gamma(link = "log")
```

Predictions are made on the response scale with the RE term excluded.

---

## Overall performance

| Metric | In-sample dev.expl | In-sample R² (full) | In-sample R² (fixed only) | LOO-CV R² | LOO-CV r | LOO-CV r (log) |
|--------|--------------------|---------------------|---------------------------|-----------|----------|----------------|
| \|d_lambda\| | 87.1% | 0.688 | 0.424 | 0.306 | 0.553 | 0.607 |
| \|d_matched_rate\| | 86.9% | 0.489 | 0.341 | 0.259 | 0.509 | 0.636 |
| \|d_rate\| | 71.4% | 0.615 | 0.203 | 0.147 | 0.383 | 0.496 |

The LOO-CV R² for lambda (0.306) attains 72% of the in-sample fixed-effects-only R² (0.424), indicating that the fixed-effect surface structure transfers well but with some loss from refitting on reduced data. The remaining gap between LOO-CV R² and the full in-sample R² (0.688) is the random effect's contribution — variance that is by definition site-specific and cannot transfer.

The pattern across metrics follows the in-sample ranking. Encounter rate has the lowest LOO performance (r = 0.38), consistent with its lower in-sample fit (71.4% dev.expl) and its substantially larger variance component (median variance fraction = 51% vs 9–10% for lambda and matched_rate).

---

## Per-dataset performance

| Dataset | r (lambda) | r (matched_rate) | r (rate) | Mean r |
|---------|-----------|-------------------|---------|--------|
| SE-grimso-high | 0.838 | 0.876 | 0.844 | 0.853 |
| SE-grimso-low | 0.842 | 0.875 | 0.790 | 0.836 |
| GE-Hunsrueck_Hochwald_NP | 0.829 | 0.868 | 0.725 | 0.807 |
| SE-Jarnashalvon | 0.868 | 0.891 | 0.613 | 0.791 |
| GE-BFNP_201920 | 0.836 | 0.832 | 0.701 | 0.790 |
| NO-nina | 0.690 | 0.818 | 0.765 | 0.758 |
| GE-Langenau | 0.805 | 0.792 | 0.581 | 0.726 |
| GE-BFNP_201819 | 0.693 | 0.728 | 0.668 | 0.697 |
| SI-serknica | 0.768 | 0.682 | 0.620 | 0.690 |
| GE-Kellerwald_Edersee_NP | 0.812 | 0.862 | 0.363 | 0.679 |
| BE-Leuven | 0.834 | 0.757 | 0.436 | 0.676 |
| SP-donana | 0.777 | 0.698 | 0.530 | 0.668 |
| SI-sinji | 0.826 | 0.507 | 0.669 | 0.668 |
| GE-Harz_NP | 0.693 | 0.736 | 0.540 | 0.656 |
| GE-Black_Forest_NP | 0.686 | 0.712 | 0.504 | 0.634 |
| GE-Hainich_NP | 0.628 | 0.674 | 0.589 | 0.630 |
| GE-Mueritz_NP | 0.650 | 0.736 | 0.495 | 0.627 |
| GE-Eifel_NP | 0.748 | 0.442 | 0.654 | 0.615 |
| SI-vrhe | 0.650 | 0.658 | 0.533 | 0.614 |
| NO-gravberget | 0.770 | 0.736 | 0.301 | 0.603 |
| SI-zabnik | 0.676 | 0.461 | 0.547 | 0.562 |
| IT-Alps | 0.791 | 0.405 | 0.218 | 0.471 |
| SI-rizana | 0.602 | 0.553 | 0.206 | 0.454 |
| NO-evenstadlia | 0.663 | 0.667 | 0.009 | 0.446 |
| PL-kampinos_NP | 0.579 | 0.623 | 0.094 | 0.432 |
| FR-montpellier | 0.359 | 0.289 | 0.150 | 0.266 |
| GE-Berchtesgaden_NP | 0.102 | 0.111 | 0.194 | 0.136 |

For lambda, 25 of 27 datasets achieve r > 0.5 and 17 exceed r > 0.65. The Swedish sites (SE-grimso-high, SE-grimso-low, SE-Jarnashalvon) and several German national parks consistently predict best across all three metrics. Per-dataset lambda and matched_rate correlations are themselves correlated (r = 0.76), while both correlate less with encounter rate (r ≈ 0.55).

Two datasets predict poorly across all metrics:

- **FR-montpellier** (r = 0.36/0.29/0.15): A Mediterranean site with two endemic species (*Genetta genetta*, *Myocastor coypus*) dropped during LOO and a distinct species community (including *Oryctolagus cuniculus*, *Lepus granatensis*, and *Herpestes ichneumon* — all single-dataset species). The model has limited training signal for Mediterranean ecology.

- **GE-Berchtesgaden_NP** (r = 0.10/0.11/0.19): An alpine site where per-species surface shapes are well-captured (within-species r = 0.60–0.94 for lambda) but deviation magnitudes are systematically wrong. Four of five species have deviations 3–11× smaller than the cross-site average, while *Rupicapra rupicapra* has deviations 9.5× larger. These opposing biases collapse the overall correlation. See § Case Study below.

---

## Per-species performance

| Species | n_obs | n_datasets | r (lambda) | r (matched) | r (rate) |
|---------|-------|-----------|-----------|------------|---------|
| Sciurus vulgaris | 11,922 | 12 | 0.77 | 0.82 | 0.41 |
| Meles meles | 14,995 | 17 | 0.64 | 0.78 | 0.57 |
| Ursus arctos | 2,338 | 2 | 0.76 | 0.78 | 0.40 |
| Lepus europaeus | 10,094 | 12 | 0.56 | 0.71 | 0.50 |
| Alces alces | 7,247 | 7 | 0.74 | 0.75 | 0.26 |
| Canis aureus | 1,057 | 3 | 0.51 | 0.40 | 0.78 |
| Vulpes vulpes | 34,656 | 25 | 0.66 | 0.55 | 0.36 |
| Bos taurus | 4,287 | 3 | 0.46 | 0.73 | 0.29 |
| Dama dama | 10,055 | 7 | 0.58 | 0.70 | 0.13 |
| Sus scrofa | 30,201 | 21 | 0.56 | 0.54 | 0.20 |
| Cervus elaphus | 29,003 | 21 | 0.47 | 0.43 | 0.32 |
| Capreolus capreolus | 36,821 | 26 | 0.50 | 0.45 | 0.27 |
| Lepus timidus | 3,835 | 5 | 0.39 | 0.52 | 0.00 |
| Procyon lotor | 2,496 | 2 | 0.55 | 0.31 | 0.04 |
| Martes martes | 856 | 5 | 0.37 | −0.11 | 0.53 |
| Rupicapra rupicapra | 3,414 | 3 | 0.08 | 0.18 | −0.33 |
| Canis lupus | 2,972 | 4 | −0.05 | 0.04 | 0.46 |
| Felis silvestris | 798 | 3 | −0.24 | −0.22 | 0.23 |

Lambda and matched_rate agree on which species predict well (*Sciurus vulgaris*, *Meles meles*, *Alces alces*, *Ursus arctos*) and which do not (*Felis silvestris*, *Canis lupus*, *Rupicapra rupicapra*). Encounter rate shows a partly different species ranking — it predicts *Canis aureus* (r = 0.78) and *Canis lupus* (r = 0.46) better than the TTE metrics, but fails for *Lepus timidus* (r = 0.00), *Dama dama* (r = 0.13), and *Procyon lotor* (r = 0.04) where the TTE metrics succeed.

Species with poor predictive performance fall into two categories:

1. **Data-limited species** (*Felis silvestris*: 798 obs, 3 datasets; *Martes martes*: 856 obs, 5 datasets): Too few training sites for the model to learn a generalisable surface.
2. **Ecologically idiosyncratic species** (*Rupicapra rupicapra*: 3,414 obs, 3 datasets; *Canis lupus*: 2,972 obs, 4 datasets): Sufficient data exists, but site-level detection ecology varies qualitatively across the species' range (see § Case Study).

---

## Performance by window duration

| Duration bin | r (lambda) | r (matched_rate) | r (rate) |
|-------------|-----------|-------------------|---------|
| 15–30 d | 0.408 | 0.350 | 0.285 |
| 31–60 d | 0.452 | 0.402 | 0.331 |
| 61–90 d | 0.479 | 0.463 | 0.350 |
| 91–120 d | 0.469 | 0.544 | 0.366 |
| 121–183 d | 0.427 | 0.575 | 0.385 |

Predictive performance for lambda peaks at 61–90 days and declines slightly at longer and shorter windows. For matched_rate, performance increases monotonically with duration. The declining lambda performance at very short windows reflects that these windows are most timing-sensitive and most variable across sites; at long windows, the signal (deviation from benchmark) shrinks toward the noise floor, reducing predictable variance.

---

## Systematic bias

| Metric | Mean bias | Median bias | % over-predicted |
|--------|-----------|-------------|------------------|
| \|d_lambda\| | −0.0029 | −0.0003 | 43.7% |
| \|d_matched_rate\| | −0.0026 | −0.0003 | 42.6% |
| \|d_rate\| | −0.0066 | −0.0008 | 45.6% |

All three metrics show a slight tendency toward under-prediction (negative mean bias), meaning the model predicts smaller deviations than observed. On the log scale, the median log(predicted/observed) ratio is −0.19 for lambda and −0.21 for matched_rate — the model tends to predict deviations ~20% too low. This is expected: the RE, which absorbs site-level deviation amplitudes, is excluded from predictions. Sites with unusually high deviations (e.g., Slovenian sites with strong ungulate rut signals) pull the mean bias negative.

---

## Cross-metric consistency

Pairwise correlation of per-dataset LOO r values:

| Pair | r |
|------|---|
| lambda ↔ matched_rate | 0.758 |
| lambda ↔ rate | 0.553 |
| matched_rate ↔ rate | 0.538 |

Lambda and matched_rate predict consistently across sites. Encounter rate's weaker cross-metric correlation reflects its different variance structure (larger stochastic component, bidirectional bias).

Datasets with r > 0.5 across metrics:

| Metric | Datasets with r > 0.5 |
|--------|----------------------|
| \|d_lambda\| | 25 / 27 |
| \|d_matched_rate\| | 22 / 27 |
| \|d_rate\| | 17 / 27 |

---

## Case study: Alpine chamois at GE-Berchtesgaden_NP

GE-Berchtesgaden_NP is the worst-predicted dataset (r = 0.10 for lambda). The cause is not a structural model failure but a site-level ecological signal that the fixed effects cannot anticipate.

### Per-species diagnosis

| Species | r (within-sp) | Mean obs | Mean pred | Obs / cross-site avg |
|---------|--------------|----------|-----------|---------------------|
| Vulpes vulpes | 0.911 | 0.0006 | 0.0064 | 0.09× |
| Capreolus capreolus | 0.937 | 0.0027 | 0.018 | 0.17× |
| Lepus timidus | 0.603 | 0.0006 | 0.0066 | 0.26× |
| Cervus elaphus | 0.791 | 0.006 | 0.018 | 0.31× |
| **Rupicapra rupicapra** | **0.641** | **0.010** | **0.0017** | **9.5×** |

Within each species, the model captures the surface shape well (r = 0.60–0.94). The problem is magnitude: four species have deviations far *below* the cross-site average, while chamois deviations are far *above*. These opposing biases cancel, collapsing the overall correlation.

### Chamois altitudinal migration

Monthly lambda values for chamois at Berchtesgaden (from 15–30 day windows) reveal a pronounced **spring peak** — detection rate in April (λ = 0.057) is 3.5× higher than in October (λ = 0.016). This timing is opposite to the autumn-rut signal that dominates other ungulates in the model.

The spring peak reflects **altitudinal migration**: chamois descend to lower elevations in winter, concentrating in the valley and forest zones where cameras are placed during March–May before ascending to high alpine pastures in summer. A secondary November peak (λ = 0.033) is consistent with the chamois rut, which peaks in November.

The same spring signal appears at IT-Alps (the only other Alpine site), where spring-start windows have lambda 3.6× higher than autumn-start windows — but at an 11× lower baseline detection rate, so absolute deviations remain small. SI-serknica (Slovenian karst, low vertical relief) shows no seasonal pattern for chamois (lambda ≈ 0.005 year-round), consistent with the absence of altitudinal migration in low-relief terrain.

| Site | Context | Mean lambda | Seasonal swing | Mean \|d_lambda\| |
|------|---------|-------------|----------------|-------------------|
| GE-Berchtesgaden_NP | Bavarian Alps, high relief | 0.018 | 3.5× (spring peak) | 0.010 |
| IT-Alps | Italian Alps, high relief | 0.0016 | ~3.5× (spring signal) | 0.0009 |
| SI-serknica | Slovenian karst, low relief | 0.0018 | ~1.3× (flat) | 0.0011 |

The deviation magnitude depends on the interaction of **population density at camera elevation** (driving baseline lambda) and **altitudinal migration amplitude** (driving seasonal swing). Both factors are site-specific and unobservable from the model's fixed covariates.

### Why the other Berchtesgaden species have low deviations

The four well-predicted species (*Vulpes vulpes*, *Capreolus capreolus*, *Cervus elaphus*, *Lepus timidus*) all show unusually small seasonal detection fluctuations at this site — deviations 3–11× below the cross-site average. This may reflect the mountainous topography dampening seasonal variation through altitudinal niche partitioning: when conditions change at camera elevation, animals can redistribute vertically rather than showing the large behavioural shifts (rut-driven activity changes, winter range contraction) seen at lowland sites.

### Implications

This case study illustrates three points:

1. **The random effect captures ecologically meaningful variation.** The RE (SD = 0.94 on the log link, corresponding to ~2.6× per SD) absorbs site-specific detection ecology — here, whether a montane species undergoes altitudinal migration through the camera-trap elevation band.

2. **Out-of-sample magnitude prediction has a principled ceiling.** The LOO-CV R² of ~31% (vs in-sample fixed-effects R² of 42%) reflects the inherent limit of predicting absolute deviation levels at unseen sites. The model can predict *which* seasons and durations are more problematic (surface shape) but not *how much* deviation to expect without site-specific calibration.

3. **The same species may require different monitoring designs at different sites.** Chamois monitoring at Berchtesgaden should avoid spring windows (March–May), while at a low-relief Slovenian site, timing matters much less for this species.

---

## Case study: Grey wolf at PL-kampinos_NP

*Canis lupus* has an overall LOO-CV r of −0.05 for lambda despite appearing in 4 base datasets (2,972 observations). As with chamois, per-dataset surface shapes are well-captured when the dataset is held out (within-dataset r = 0.48–0.89), but magnitude prediction fails.

### Cross-site detection ecology

| Site | N obs | Mean lambda | CV(lambda) | Mean \|d_lambda\| |
|------|-------|-------------|-----------|-------------------|
| PL-kampinos_NP | 1,073 | 0.0058 | 0.64 | 0.0042 |
| SI-sinji | 312 | 0.0029 | 0.59 | 0.0020 |
| SI-serknica | 1,033 | 0.0021 | 0.44 | 0.0010 |
| GE-Mueritz_NP | 554 | 0.0015 | 0.21 | 0.0004 |

PL-kampinos_NP has 2.8× higher baseline detection than the other sites and 3–12× larger deviations. Its seasonal profile shows a pronounced **summer peak** (June–August lambda ≈ 0.012–0.016), likely reflecting increased wolf activity during pup-rearing season. The Slovenian sites have intermediate detection rates with weaker seasonality, and GE-Mueritz_NP is nearly flat (CV = 0.21).

### LOO prediction failure

When Kampinos is held out, the model trains on two Slovenian and one German dataset — all with low, relatively stable wolf detection. It then predicts Kampinos at λ ≈ 0.0002, a 22-fold under-prediction (observed mean = 0.0042). This single dataset (36% of all wolf observations) dominates the overall negative correlation.

Kampinos NP is a lowland forest near Warsaw with a dense, recovering wolf population — an ecological context absent from the training data. The high detection rate and strong summer seasonality are genuine population-level signals, not measurement artefacts.

### Contrast with chamois

Both case studies share the same mechanism (site-level population density and ecology drive magnitude mismatch), but differ in the ecological driver. For chamois, the issue is **topographic context** (altitudinal migration through the camera-trap zone). For wolves, it is **population density and recovery status** — Kampinos wolves are simply far more abundant and active than wolves at the other three sites.

---

## Case study: European wildcat (*Felis silvestris*) — data limitation

*Felis silvestris* has the worst per-species LOO-CV correlation (r = −0.24 for lambda), but unlike the chamois and wolf cases, this reflects **data sparsity** rather than ecological idiosyncrasy.

The species appears in only 3 datasets: GE-Hainich_NP (723 obs, 91% of total), SI-serknica (73 obs from one slice), and GE-Eifel_NP (2 obs). When Hainich is held out, the model trains on just 75 wildcat observations — far too few for a meaningful species-specific tensor surface. The resulting predictions for Hainich are weakly negatively correlated with observed values (r = −0.30), driving the overall poor performance.

More fundamentally, wildcat deviations are extremely small: mean |d_lambda| = 0.0005, only 17% of the all-species median. Detection rate varies minimally across seasons (lambda range: 0.0009–0.003 at Hainich). At this scale, the signal is close to the noise floor and the model has little predictable variance to work with.

This case does not represent a model failure — it represents the practical limit of LOO-CV for rare, low-detectability species with insufficient geographic replication. Among the 18 species evaluated, it is the only one where the negative correlation is attributable to data limitation rather than ecological mismatch.

---

## Interpretation

The LOO-CV results support two conclusions:

**The sensitivity surface shape generalises well.** For the majority of datasets and species, the model predicts out-of-sample with moderate to good accuracy (median per-dataset r = 0.68 for lambda). The species-specific tensor surfaces — which encode how timing sensitivity varies with duration and how this pattern differs across species — transfer to unseen sites. This means the qualitative insights from the model (which seasons are risky, how much duration helps, which species are most sensitive) are not overfit to the training data.

**Absolute deviation magnitude requires site-specific information.** The gap between in-sample R² (0.69 full / 0.42 fixed-effects) and LOO-CV R² (0.31) is almost entirely attributable to the dataset × species random intercept. This component captures real ecological variation (site-level population density, local phenology, topographic context) that is invisible to the fixed effects. The site-level variation captured by the RE reflects ecological processes (e.g., altitudinal migration, local population density) that cannot be predicted from broad-scale covariates. For practical applications, the model's predictions should be interpreted in relative terms (surface shape, ranking of protocols, species-level comparisons) rather than as absolute deviation forecasts for a specific new site.

---

## Output files

| File | Content |
|------|---------|
| `loo_cv_predictions.csv` | Observation-level predictions (lambda) |
| `loo_cv_predictions_matched_rate.csv` | Observation-level predictions (matched_rate) |
| `loo_cv_predictions_rate.csv` | Observation-level predictions (rate) |
| `loo_cv_overall_summary.csv` | Overall performance (3 metrics) |
| `loo_cv_by_dataset.csv` | Per-dataset summary (lambda) |
| `loo_cv_by_dataset_all_metrics.csv` | Per-dataset summary (3 metrics) |
| `loo_cv_by_species_all_metrics.csv` | Per-species summary (3 metrics) |
| `figures/FigS_loo_cv_diagnostics.pdf` | 4-panel diagnostic (lambda) |
| `figures/FigS_loo_cv_cross_metric.pdf` | Cross-metric comparison |

Script: `loo_cv_analysis.R` (~150 min runtime).

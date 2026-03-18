# Appendix: Species-Inclusion Threshold Robustness

Species × window combinations enter the analysis only if they pass three minimum viable sample thresholds: minimum independent detection events, minimum cameras with at least one detection, and minimum occasion periods with detections. These thresholds prevent unreliable rate estimates from entering the sensitivity surface models, but could bias results by selectively excluding certain species, datasets, or window positions. We tested whether varying each threshold changes the main findings.

For each threshold variant, we re-prepared the species-level data and re-fitted the primary model (M6: species-specific 2D tensor product surfaces for |Δλ|, Gamma family with log link). We compared deviance explained, parametric coefficient estimates, and the shape of predicted surfaces (Pearson correlation of the full 2D surface and 60-day seasonal profile, evaluated on species shared across all threshold levels).

---

## A. Minimum events threshold

**Baseline: 20 events.** Tested at 10 (lenient) and 30 (strict). The lenient threshold required re-running the data pipeline; the strict threshold was applied as a post-hoc filter.

### Data retention

| Threshold | Species | Rows | Change vs baseline |
|-----------|---------|------|--------------------|
| Lenient (10) | 28 | 149,537 | +25,635 rows, +2 species |
| **Baseline (20)** | **26** | **123,902** | — |
| Strict (30) | 26 | 106,902 | −17,000 rows, same species |

The lenient threshold gains two additional species. The strict threshold loses no species entirely; the ~17,000 dropped rows come predominantly from short windows and rare species where event counts are marginal.

### Model comparison

| Threshold | Dev. expl. | R² | Surface *r* | Seasonal *r* (60 d) |
|-----------|------------|----|-------------|----------------------|
| Lenient (10) | 87.0% | 0.751 | 0.987 | 0.997 |
| **Baseline (20)** | **86.8%** | **0.752** | **1.000** | **1.000** |
| Strict (30) | 86.3% | 0.752 | 0.970 | 0.983 |

Deviance explained varies by less than 1 percentage point. Predicted surface correlations are ≥ 0.970 in all comparisons. All five parametric covariates (temperature seasonality, log trap-days, log number of sites, latitude, trap array size) retain consistent signs and comparable magnitudes; no coefficient changes significance direction.

**Fig. S2** shows duration curves and 60-day seasonal profiles for all three thresholds.

---

## B. Minimum positive sites threshold

**Baseline: 5 cameras.** Tested at 10 and 15 via post-hoc filtering. A lenient threshold of 3 was evaluated diagnostically but not pursued: it would add only ~3,100 rows (2.5%) with no new species surviving the minimum-rows-per-species filter.

### Data retention

| Threshold | Species | Rows | Change vs baseline |
|-----------|---------|------|--------------------|
| **Baseline (5)** | **26** | **123,902** | — |
| Strict (10) | 24 | 110,342 | −11%, −2 species |
| Very strict (15) | 23 | 89,534 | −28%, −3 species |

Species lost entirely at threshold 10: *Myocastor coypus* and *Canis aureus*. Species additionally lost at 15: *Martes martes*. Losses are concentrated among spatially sparse species — *Genetta genetta* loses 65% of observations, *Oryctolagus cuniculus* 53%, *Lepus europaeus* 47% — while core ungulate and generalist species lose < 7%.

### Model comparison

| Threshold | Dev. expl. | R² | Surface *r* | Seasonal *r* (60 d) |
|-----------|------------|----|-------------|----------------------|
| **Baseline (5)** | **86.8%** | **0.752** | **1.000** | **1.000** |
| Strict (10) | 86.8% | 0.794 | 0.974 | 0.977 |
| Very strict (15) | 85.8% | 0.792 | 0.925 | 0.994 |

Deviance explained is virtually unchanged at threshold 10 and declines by only 1 percentage point at threshold 15 despite a 28% data reduction. R² increases at stricter thresholds, consistent with removal of noisier marginal observations. The modest surface correlation reduction at threshold 15 (0.925) reflects a vertical level shift (higher mean |Δλ| after removing low-detection observations) rather than a change in surface shape, as confirmed by the near-perfect seasonal profile correlation (0.994).

### Parametric coefficients

| Covariate | Baseline (5) | Strict (10) | Very strict (15) |
|-----------|-------------|------------|-------------------|
| Temperature seasonality | −0.160 (*p* = .069) | −0.157 (*p* = .075) | −0.163 (*p* = .061) |
| Log trap-days | 0.728 | 0.782 | 0.796 |
| Log n sites | −1.17 | −1.13 | −1.15 |
| Latitude | −0.301 (*p* < .001) | −0.273 (*p* = .003) | −0.191 (*p* = .042) |
| Trap array size | −0.133 (*p* = .024) | −0.161 (*p* = .006) | −0.161 (*p* = .005) |

All coefficients retain the same sign. The latitude effect attenuates from −0.30 to −0.19 at the strictest threshold — expected because stricter filtering disproportionately removes spatially sparse species at northern sites — but remains significant.

**Fig. S3** shows duration curves and 60-day seasonal profiles for all three thresholds.

---

## C. Minimum positive occasions threshold

**Baseline: 5 occasions.** This threshold is non-binding: no species × window combination in the dropped data failed exclusively on it. Whenever a species-window has too few positive occasions, it also fails on the events or sites threshold. The minimum value in the retained data is 7 (median: 41), well above the threshold of 5. Even tripling the threshold to 15 would lose only 3.1% of rows with no species dropped. No model re-fitting is warranted.

---

## Conclusion

The sensitivity surface results are robust to species-inclusion thresholds. Deviance explained varies by ≤ 1 percentage point, parametric coefficient signs and significance are preserved, and predicted surface shapes are highly correlated with baseline (*r* ≥ 0.925 for the full surface, *r* ≥ 0.977 for seasonal profiles). The occasions threshold has no influence on data inclusion. Differences between threshold variants are driven by which rare, spatially sparse species enter the analysis, not by changes in the modelled relationship between window design and detection deviation.

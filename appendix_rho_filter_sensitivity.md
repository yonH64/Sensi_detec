# Appendix: Rank Correlation Filter Sensitivity

The rank correlation model (rho_lambda) uses beta regression on Spearman's ρ between sub-window and full-year detection rate rankings. With few shared species, Spearman's ρ has a small discrete support: 3 species yield only 7 possible values, 4 species yield 15, and 5 species yield 16. This coarse discretization produces substantial boundary mass — a large fraction of observations effectively at ρ = 1 — which violates beta regression assumptions and inflates deviance explained via artefactual perfect correlations rather than genuine surface structure.

The baseline analysis filters to windows with ≥ 5 shared species. Here we test cutoffs of 3, 4, 5, 6, and 8 to assess how the filter choice affects data retention, boundary mass, model fit, and the predicted surface shape.

---

## Data retention and boundary mass

| Min shared species | Observations | Datasets | % at boundary (ρ ≈ 1) | Dev. explained |
|-------------------|-------------|----------|------------------------|----------------|
| ≥ 3 | 39,862 | 35 | 44.5% | 19.0% |
| ≥ 4 | 32,750 | 34 | 40.0% | 30.6% |
| **≥ 5 (baseline)** | **23,943** | **30** | **33.8%** | **36.3%** |
| ≥ 6 | 17,530 | 27 | 29.3% | 43.0% |
| ≥ 8 | 5,424 | 12 | 15.6% | 55.3% |

Boundary mass and deviance explained move in opposite directions: each increment in the cutoff removes more boundary-inflated observations and allows the model to explain a larger share of the remaining variation. This tradeoff is the central tension — stricter filtering improves model quality but reduces data volume and geographic coverage.

### Note on boundary rates

Boundary mass (ρ ≈ 1) is not simply a function of small N — it alternates between even and odd numbers of shared species. With even N, rho = 1 occurs exactly when rankings are identical. With odd N, floating-point arithmetic produces values within 10⁻¹⁶ of 1 for identical rankings, which are functionally equivalent. The boundary rates by exact N are:

| N shared | Observations | % at ρ ≈ 1 |
|----------|-------------|------------|
| 3 | 7,112 | 65% |
| 4 | 8,807 | 57% |
| 5 | 6,413 | 46% |
| 6 | 6,988 | 41% |
| 7 | 5,118 | 28% |
| 8 | 2,385 | 19% |
| 9 | 2,135 | 14% |
| 10 | 677 | 6% |

---

## Surface comparison

| Cutoff | Surface *r* vs ≥ 5 baseline |
|--------|----------------------------|
| ≥ 3 | 0.911 |
| ≥ 4 | 0.937 |
| ≥ 6 | 0.953 |
| ≥ 8 | 0.916 |

Surface correlations range from 0.91 to 0.95, lower than for any other sensitivity check in this study. The surface shape is qualitatively preserved — all cutoffs show higher predicted ρ for longer windows and moderate seasonal variation — but the absolute predicted values shift. At ≥ 3, boundary mass compresses the predicted surface toward ρ ≈ 0.87 with a range of only 0.095. At ≥ 8, the surface spans ρ ≈ 0.72–1.00 (range 0.28), revealing timing × duration structure that is masked by boundary inflation at lower cutoffs.

The drop in correlation at ≥ 8 (r = 0.916, lower than ≥ 6) is driven by the severe data reduction: only 12 datasets and 5,424 observations survive, and the model becomes sensitive to which few datasets dominate.

---

## Justification for ≥ 5 shared species

The ≥ 5 cutoff balances three considerations:

1. **Boundary mass.** At ≥ 3, nearly half of observations are at the boundary — far too much for beta regression. The ≥ 5 cutoff reduces this to 34%, which is still high but manageable with the 0.001/0.999 squeezing applied before fitting. Each additional shared species yields diminishing returns in boundary reduction.

2. **Geographic coverage.** The ≥ 5 cutoff retains 30 of 35 datasets. At ≥ 6 this drops to 27, and at ≥ 8 to just 12 — fewer than half the study's datasets, severely limiting generalizability.

3. **Model stability.** Surface correlations with the ≥ 5 baseline peak at ≥ 6 (r = 0.953) and then decline at ≥ 8 (r = 0.916), indicating that stricter filtering introduces instability from small-sample effects that offsets the benefit of reduced boundary mass.

The ≥ 5 cutoff sits at the elbow where deviance explained begins accelerating (19% → 31% → 36%) while data retention and geographic coverage remain adequate. It is not formally optimal, but the sensitivity analysis confirms that no cutoff between 3 and 8 qualitatively changes the predicted surface shape or the conclusion that rank preservation is high (median predicted ρ ≈ 0.93 at the ≥ 5 cutoff) and increases monotonically with window duration.

---

## Residual concern

Even at ≥ 5, one-third of observations are at the upper boundary. The beta regression model handles this through the 0.001/0.999 squeezing, but a zero-one-inflated beta (ZOIB) model or a hurdle approach would be more principled. Such models are not currently available in `mgcv::bam()` with the tensor product smooth structure required here. The deviance explained for the rho model (36.3%) is the lowest of all seven models in the study, and some of this weakness is likely attributable to residual boundary effects rather than true lack of surface structure.

**Figure:** `figures/rho_filter_sensitivity.pdf`
**Output data:** `rho_filter_sensitivity_summary.csv`

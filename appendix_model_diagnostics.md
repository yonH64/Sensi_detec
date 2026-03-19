# Appendix: Model Specification Diagnostics

The primary model (M6) uses a Gamma family with log link, species-specific 2D tensor product surfaces with basis dimensions k = (8, 6), linear parametric covariates, and a single random intercept per dataset × species combination. Here we test whether alternative model specifications change the conclusions.

For each variant, we re-fitted the model on the same baseline data (222,749 observations, 29 species) and compared deviance explained, AIC, and predicted surface shape.

---

## A. Basis dimension adequacy (k-doubling)

**Concern:** The tensor product basis dimensions k = (8, 6) — yielding up to 48 basis functions per species surface — may be too restrictive to capture fine-grained timing × duration patterns.

**Test:** Doubled the basis to k = (16, 12), yielding up to 192 basis functions per species.

| Specification | Dev. expl. | AIC | Surface *r* |
|---------------|------------|-----|-------------|
| **Baseline k = (8, 6)** | **87.1%** | **−2,258,933** | **1.000** |
| Doubled k = (16, 12) | 87.3% | −2,260,000 | 0.996 |

Doubling the basis gains 0.2 percentage points of deviance explained and an AIC improvement of ~1,067. The predicted surface correlation is 0.996 — the additional flexibility captures minor fine-structure but does not alter the main surface features. This is consistent with the prior k-check (all k-index values ≥ 0.99, all p > 0.69) which indicated the baseline basis was adequate.

The marginal AIC improvement is offset by substantially increased computation time (17 minutes vs 2.5 minutes) and a larger model object, with no qualitative change in interpretation. The baseline k = (8, 6) is retained.

---

## B. Response distribution (Gamma vs Tweedie)

**Concern:** The Gamma family assumes a specific mean–variance relationship (variance ∝ mean²). The data might be better served by the more flexible Tweedie family, which estimates the power parameter *p* from data (Gamma corresponds to *p* = 2; compound Poisson–gamma to 1 < *p* < 2).

**Test:** Fitted the same M6 structure with `family = tw()`, which estimates *p* via maximum likelihood.

| Specification | Family | Estimated *p* | Surface *r* |
|---------------|--------|---------------|-------------|
| **Baseline** | **Gamma (log)** | **2 (fixed)** | **1.000** |
| Tweedie | Tweedie | 1.990 | 0.997 |

The Tweedie family estimated *p* = 1.990 — essentially Gamma. The predicted surface correlation is 0.997. AIC values are not directly comparable across families (different likelihoods), but the convergence to *p* ≈ 2 indicates the Gamma is the correct choice. There is no evidence of a compound Poisson component (which would push *p* below ~1.8).

---

## C. Random effects structure (nested vs flat)

**Concern:** Four multi-year base datasets (BE-Leuven, NO-evenstadlia, SI-serknica, SP-donana) contribute multiple annual slices — 13 slices total, accounting for 44.7% of observations. The current model uses a single random intercept per dataset-slice × species (254 levels). This treats slices from the same site as independent, ignoring potential site-level correlation.

**Test:** Added a nested random intercept for base-dataset × species (183 levels) alongside the existing slice-level term.

| Specification | RE levels | ΔAIC | Surface *r* |
|---------------|-----------|------|-------------|
| **Baseline (slice-level only)** | **254** | **0** | **1.000** |
| Nested (base + slice) | 183 + 254 | −5 | 0.985 |

The base-dataset random effect has substantial variance (SD = 0.845) compared to the slice-level term (SD = 0.350), confirming that slices from the same site share a common deviation level. However, this additional variance component produces a negligible AIC improvement (ΔAIC = −5) because the slice-level term already captures most of the site effect: each slice inherits its site's characteristics through the data, and the penalisation of the 254 slice-level random intercepts effectively pools toward the site mean.

The surface correlation of 0.985 — the lowest of any variant tested — reflects the redistribution of variance between the two RE terms rather than a change in the fixed-effect surface. Fixed-effect coefficients are virtually unchanged.

---

## D. Nonlinear environmental effect (smooth BIO4)

**Concern:** Temperature seasonality (BIO4) enters the model as a linear parametric term. The earlier observation that Slovenian sites drive a mid-seasonality deviation hump (see AGENTS.md § "Slovenian sites drive mid-seasonality band deviations") raises the question of whether a smooth (nonlinear) BIO4 effect better fits the data.

**Test:** Replaced the linear `s_bio4` term with a thin plate regression spline `s(s_bio4, k = 5)`.

| Specification | BIO4 edf | ΔAIC | BIO4 *p*-value |
|---------------|----------|------|----------------|
| **Baseline (linear)** | **1** | **0** | **0.125** |
| Smooth (k = 5) | 1.00 | 0 | 0.125 |

The smooth collapsed to an effective degrees of freedom of exactly 1.00 — a straight line. The AIC is identical to the baseline. The penalisation correctly determined that no curvature is supported by the data. The mid-seasonality hump observed in the raw data is a site-level effect (Slovenian ungulate rut) absorbed by the species-specific surfaces and random effects, not a systematic nonlinear BIO4 relationship.

---

## Summary

| Diagnostic | Change from baseline | Conclusion |
|------------|---------------------|------------|
| k-doubling | +0.2% dev.expl, *r* = 0.996 | Basis adequate |
| Tweedie family | *p* = 1.99 ≈ Gamma, *r* = 0.997 | Gamma validated |
| Nested RE | ΔAIC = −5, *r* = 0.985 | Slice-level RE sufficient |
| Smooth BIO4 | edf = 1.00, ΔAIC = 0 | Linear term adequate |

No model specification variant produces a qualitative change in the sensitivity surface or its interpretation. The baseline M6 specification is well-supported.

**Output data:** `model_diagnostics_phase3.csv`

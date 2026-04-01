# Appendix: Inverse-Slice Weighting for Multi-Year Site Overrepresentation

Four camera-trap sites contributed multiple annual slices to the dataset: BE-Leuven (5 slices), SP-donana (4), SI-serknica (2), and NO-evenstadlia (2). Together, these 13 slices account for 44.7% of all observations (99,467 / 222,748 rows). Because slices from the same physical location share identical environmental covariates (latitude, trap array size, etc.), they inflate the effective sample size for parametric coefficients without adding independent geographic information. This is particularly relevant given that the LOO influence analysis identified SP-donana (lowest latitude, 16.7% of rows) and BE-Leuven (15.6% of rows) as high-leverage datasets for the parametric coefficients.

We tested whether downweighting multi-year sites stabilises the parametric coefficients and preserves the species-specific sensitivity surfaces.

---

## Weighting scheme

Each observation receives weight w = 1 / n_slices, where n_slices is the number of temporal slices from its base dataset. Single-year datasets retain weight 1.

| Base dataset | Slices | Weight | Rows | Effective rows |
|---|---|---|---|---|
| BE-Leuven | 5 | 0.20 | 34,694 | 6,939 |
| SP-donana | 4 | 0.25 | 37,110 | 9,278 |
| SI-serknica | 2 | 0.50 | 19,534 | 9,767 |
| NO-evenstadlia | 2 | 0.50 | 8,129 | 4,065 |
| All others (22 datasets) | 1 | 1.00 | 123,281 | 123,281 |
| **Total** | | | **222,748** | **153,329** |

Only auto-sliced datasets (identified by `_slice\d+$` suffix) are grouped. GE-BFNP_201819/201920 (same park, consecutive years, slightly different camera centroids) and SE-grimso-high/low (same location, different elevation strata) are not grouped — they were entered as separate datasets, not produced by the temporal slicing algorithm.

Weights were passed to `mgcv::bam()` via the `weights` argument. All other model settings (formula, family, knots, method, discrete fitting) were identical to the baseline M6 model.

---

## Parametric coefficient comparison (abs_d_lambda)

| Term | β (baseline) | β (weighted) | % change | p (baseline) | p (weighted) |
|------|---|---|---|---|---|
| s_latitude | −0.456 | −0.483 | −6.0% | 6.1e-11 | 6.5e-11 |
| s_trap_array | −0.134 | −0.136 | −1.8% | 3.5e-02 | 4.3e-02 |
| l_nsites | −1.21 | −1.37 | −13.0% | <2e-16 | <2e-16 |
| l_trapdays | +0.82 | +0.98 | +19.6% | <2e-16 | <2e-16 |

All four coefficients retain the same sign and significance level. Latitude **strengthens slightly** (−0.456 → −0.483), confirming the effect is not inflated by SP-donana's overrepresentation. The effort controls (`l_nsites`, `l_trapdays`) shift by 13–20% — larger than the environmental covariates but consistent with those terms absorbing some of the variance previously explained by the downweighted multi-year site observations.

---

## Consistency across detection metrics

The same weighting was applied to all three absolute detection deviation models.

### Latitude coefficient

| Metric | β (baseline) | β (weighted) | % change |
|--------|---|---|---|
| abs_d_lambda | −0.456 | −0.483 | −6.0% |
| abs_d_rate | −0.269 | −0.279 | −3.6% |
| abs_d_matched_rate | −0.398 | −0.411 | −3.3% |

The latitude coefficient is stable across all three metrics (3–6% change, same sign, same significance).

### Trap array coefficient

| Metric | β (baseline) | β (weighted) | % change | p (baseline) | p (weighted) |
|--------|---|---|---|---|---|
| abs_d_lambda | −0.134 | −0.136 | −1.8% | 0.035 | 0.043 |
| abs_d_rate | +0.016 | −0.008 | sign flip | 0.80 | 0.89 |
| abs_d_matched_rate | +0.004 | +0.030 | +590% | 0.94 | 0.61 |

The `s_trap_array` coefficient is non-significant for abs_d_rate (p ≈ 0.80–0.89) and abs_d_matched_rate (p ≈ 0.94–0.61) in both models. The large percentage changes reflect noise around zero, not meaningful instability. Only for abs_d_lambda is the coefficient marginally significant (p ≈ 0.04), where it is stable (−1.8% change).

---

## Surface correlations

| Metric | r (overall) | r (min species) | r (median species) | Dev.expl (base → wt) |
|--------|---|---|---|---|
| abs_d_lambda | 0.990 | 0.975 | 0.997 | 87.1% → 84.2% |
| abs_d_rate | 0.992* | −0.012* | 0.977 | 71.4% → 72.6% |
| abs_d_matched_rate | 0.992 | 0.979 | 0.998 | 86.9% → 83.4% |

*For abs_d_rate, two marginal species distort the raw overall correlation (r = 0.490): Canis aureus (r = −0.01; 312 observations, almost exclusively in SI-eow, a single-slice dataset — weighting does not change its data, only the global intercept shifts its predicted surface) and Castor fiber (r = 0.67; 21 observations, only in BE-Leuven_slice5 — effective sample drops to ~4 with weight 0.2). Excluding these two species: r = 0.992. Spearman rank correlation across all 29 species: 0.953. Of 29 species, 23 have per-species r > 0.90 and 18 have r > 0.95.

Deviance explained drops by 2–3 percentage points for lambda and matched_rate, consistent with the reduced effective sample size (31% fewer effective observations). The rate model is essentially unchanged (71.4% → 72.6%).

---

## Community models

The same weighting was applied to the three community-level models (d_sr_raref, prop_sr_full, rho_lambda), which use `sens_richness` data at the dataset × window level. Effective sample sizes: d_sr_raref and prop_sr_full: 34,450 (from 46,375); rho_lambda: 14,978 (from 23,943). Community models lack a species dimension, so there is no `s_trap_array` term and no per-species surface breakdown.

### Latitude coefficient

| Metric | β (baseline) | β (weighted) | % change | p (baseline) | p (weighted) |
|--------|---|---|---|---|---|
| d_sr_raref | −0.037 | −0.039 | −5.5% | 5.3e-01 | 5.4e-01 |
| prop_sr_full | −0.203 | −0.171 | +16% | 3.5e-03 | 5.5e-03 |
| rho_lambda | −0.095 | −0.092 | +3.4% | 3.6e-01 | 5.1e-01 |

Latitude retains the same sign in all three models. For prop_sr_full (the only community model where latitude is significant), the coefficient weakens by 16% but remains significant (p = 0.006). For d_sr_raref and rho_lambda, latitude is non-significant in both baseline and weighted models — the changes are noise.

### Effort covariates

| Metric | Term | β (baseline) | β (weighted) | % change |
|--------|------|---|---|---|
| d_sr_raref | l_nsites | +0.398 | +0.374 | −6% |
| d_sr_raref | l_trapdays | +0.017 | +0.086 | +404%* |
| prop_sr_full | l_nsites | +0.500 | +0.770 | +54% |
| prop_sr_full | l_trapdays | +1.53 | +1.70 | +11% |
| rho_lambda | l_nsites | +0.227 | +0.438 | +93% |
| rho_lambda | l_trapdays | +0.917 | +1.21 | +32% |

*The 404% change for l_trapdays in d_sr_raref reflects a near-zero baseline (0.017), not a meaningful shift. All effort coefficients retain the same sign.

### Surface correlations

| Metric | r (overall) | Dev.expl (base → wt) |
|--------|---|---|
| d_sr_raref | 0.986 | 24.7% → 23.7% |
| prop_sr_full | 0.982 | 77.3% → 78.7% |
| rho_lambda | 0.991 | 36.3% → 40.5% |

All three community surfaces are well preserved (r ≥ 0.982). The rho_lambda model shows a modest improvement in deviance explained (36.3% → 40.5%), likely because downweighting BE-Leuven's 5 slices (which had high boundary mass in rho) reduces artefactual noise.

---

## Conclusion

Inverse-slice weighting confirms that the environmental coefficients — particularly latitude — are not artefacts of multi-year site overrepresentation. Across all six models (3 detection + 3 community), latitude retains the same sign and either strengthens or changes by <16%. The species-specific detection surfaces are preserved (median per-species r ≥ 0.977), and the community surfaces are similarly robust (r ≥ 0.982). The downweighting reduces effective sample size by 26–37%, which explains modest shifts in deviance explained, but qualitative conclusions are unchanged.

This check complements the LOO influence analysis (which identified *which* datasets have high leverage) by demonstrating that equalising their influence does not alter the results. The combination of LOO stability (25/26 iterations retain significant latitude, same sign) and inverse-slice weighting stability (3–6% coefficient change) provides strong evidence that the latitude effect is a genuine feature of the data, not a sampling artefact.

**Output files:** `inverse_slice_weighting_summary.csv`, `inverse_slice_weighting_coefficients.csv`

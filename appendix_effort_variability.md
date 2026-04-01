# Appendix: Within-Slice Effort Variability

Camera effort is not uniformly distributed across the 12-month benchmark period in all datasets. Deployments may start or end mid-year, cameras may fail and not be replaced, or new cameras may be added partway through. If the resulting effort imbalance is severe, the FULL benchmark could be weighted toward the season with more active cameras, potentially biasing the deviation metrics.

We quantified effort variability within each dataset-slice as the coefficient of variation (CV) of monthly active camera counts, estimated from the number of active sites (`n_sites`) across 29-day sliding windows.

---

## Effort CV distribution

| Statistic | Value |
|-----------|-------|
| Mean | 0.215 |
| SD | 0.147 |
| Range | 0.005 – 0.516 |
| Datasets with CV > 0.3 | 11 of 35 |

The 11 high-CV datasets include SI-eow (0.52), SP-donana slices 2–4 (0.32–0.50), SI-abnik (0.43), GE-Langenau (0.42), IT-Alps (0.34), SI-serknica_slice2 (0.32), NO-evenstadlia_slice2 (0.32), SE-grimso-high (0.32), and SE-grimso-low (0.31). Effort variability arises from different mechanisms: monotonic camera attrition (SI-eow drops from 12 to 4 cameras over 12 months), mid-year array expansion (SP-donana_slice2 doubles from 35 to 60 cameras in June), and seasonal camera failures.

---

## Residual diagnostic

We extracted deviance residuals from the baseline M6 model (abs_d_lambda, Gamma family, 222,748 observations) and compared high-CV (> 0.3) versus low-CV (≤ 0.3) datasets.

**Overall bias.** Mean deviance residuals were −0.080 (high-CV, N = 69,232) versus −0.092 (low-CV, N = 153,516). The Pearson correlation between effort CV and per-dataset mean residual was r = 0.14 (p = 0.43). Residual standard deviations were slightly smaller for high-CV datasets (0.50 vs 0.55), ruling out heteroscedasticity.

**Seasonal asymmetry.** If an effort-biased benchmark inflated or deflated deviations in specific seasons, residuals should show a consistent seasonal sign pattern in high-CV datasets. No such pattern was found. Mean residuals by season differed by less than 0.03 between the two groups, with no consistent direction. Per-dataset seasonal residuals were idiosyncratic: for example, SP-donana_slice2 showed more negative residuals in H1 (Jan–May, ~35 cameras) than H2 (Jun–Dec, ~58 cameras), but slices 3 and 4 — which have even effort — showed the same asymmetry in the opposite direction, indicating year-specific ecological variation rather than effort bias.

**Deviance contribution.** Datasets were ranked by leverage ratio (percentage of total squared deviance residuals divided by percentage of observations). SI-eow had the highest leverage ratio (2.49), but the next-highest were GE-BFNP_201819 (1.99, effort CV = 0.005) and NO-nina (1.95, CV = 0.21). The correlation between leverage ratio and effort CV was r = −0.18 (p = 0.31). Model misfit is not concentrated in high-CV datasets.

---

## Formal covariate test

Effort CV was standardised (centred at 0.215, scaled by SD = 0.147) and added as a linear parametric term to the M6 model for each of the three detection metrics.

| Metric | ΔAIC | Δ Dev. expl. | β(effort CV) | SE | p |
|--------|------|--------------|--------------|-----|---|
| abs_d_lambda | −2.3 | +0.000% | +0.271 | 0.066 | 3.6 × 10⁻⁵ |
| abs_d_rate | −0.2 | +0.000% | −0.042 | 0.066 | 0.52 |
| abs_d_matched_rate | −1.6 | +0.000% | +0.125 | 0.059 | 0.034 |

Deviance explained was identical to the second decimal place for all three metrics. The ΔAIC values (−0.2 to −2.3) are negligible given the sample size. The effort CV coefficient was nominally significant for lambda and matched_rate but non-significant and oppositely signed for encounter rate, indicating the signal is not consistent across metrics measuring the same underlying quantity. Latitude coefficients were stable: β(latitude) shifted from −0.456 to −0.401 for lambda, with comparable stability for rate (−0.269 to −0.278) and matched_rate (−0.398 to −0.372).

---

## SI-eow: highest effort CV and largest mean deviation

SI-eow (CV = 0.52, mean |d_lambda| = 0.032) was examined as a worst case. This site has 52 camera locations but only 4–12 were active within its 12-month slice, declining monotonically from August 2023 to July 2024. Despite the extreme effort imbalance, several diagnostics indicate the large deviations are driven by ecology and small sample size rather than benchmark bias:

1. **Metric agreement.** d_lambda and d_matched_rate correlated at r = 0.76–0.97 across species. The matched metric, which controls for camera composition, produces comparable or larger deviations than the unmatched metric.
2. **Species specificity.** Non-ungulate species at SI-eow had typical deviations (mean |d_lambda| = 0.009), matching the study average. Only the three ungulates (Capreolus capreolus, Cervus elaphus, Sus scrofa) were extreme, with Sus scrofa showing a 28.5-fold seasonal swing in lambda — far exceeding any other site.
3. **Residual direction.** The model overpredicted deviations at SI-eow (mean residual = −0.23), the opposite of what effort-driven benchmark bias would produce.
4. **Spatial coverage saturation.** With so few cameras, Capreolus hit spatial_cov_m = 1 (all cameras detecting) in 37% of windows, Cervus in 34%. This ceiling effect is a small-sample-size artefact, not an effort-distribution problem.

---

## Conclusion

Within-slice effort variability does not introduce systematic bias into the sensitivity surface. The three-metric covariate test, the residual diagnostics, and the SI-eow case study all converge on the same conclusion: the species-specific tensor product surfaces and dataset × species random effects in the M6 model adequately absorb any detection signal that effort variability might proxy for. Effort variability need not be included as a covariate or used as a dataset exclusion criterion.

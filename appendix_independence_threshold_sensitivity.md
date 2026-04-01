# Appendix: Independence Threshold Sensitivity

Raw camera-trap detections are collapsed into independent events by requiring a minimum time gap between consecutive detections of the same species at the same camera. Events separated by less than this gap are treated as a single detection. The baseline uses a 30-minute independence threshold. This is a standard choice in camera-trap studies, but values of 15 and 60 minutes are also common.

We tested whether the independence gap affects the sensitivity surface by re-running the full data pipeline for all 28 datasets at 15-minute and 60-minute gaps, then re-preparing the species-level data and re-fitting the primary model (M6).

---

## Data retention

| Independence gap | Species | Rows | Datasets |
|------------------|---------|------|----------|
| 15-min | 29 | 223,571 | 35 |
| **30-min (baseline)** | **29** | **222,748** | **35** |
| 60-min | 28 | 221,572 | 35 |

Shorter gaps produce slightly more rows (more events → more species passing thresholds in marginal windows). The 60-minute gap loses one species (*Ovis aries*) that becomes too sparse. The differences are small: ±0.4% in row count relative to baseline.

---

## Model comparison

| Gap | Dev. expl. | AIC | Surface *r* vs baseline |
|-----|------------|-----|-------------------------|
| 15-min | 87.1% | −2,269,125 | 0.9999 |
| **30-min (baseline)** | **87.1%** | **−2,258,933** | **1.000** |
| 60-min | 87.1% | −2,244,052 | 0.9998 |

Deviance explained is identical (87.1%) across all three thresholds. Surface correlations exceed **0.999** — the highest of any sensitivity check in this study. The three predicted surfaces are functionally indistinguishable.

---

## Parametric coefficients

| Covariate | 15-min | 30-min | 60-min |
|-----------|--------|--------|--------|
| Log trap-days | 0.818 | 0.819 | 0.819 |
| Log n sites | −1.21 | −1.21 | −1.21 |
| Latitude | −0.358 | −0.350 | −0.334 |
| Trap array size | −0.146 | −0.150 | −0.155 |

All coefficients retain the same sign and nearly identical magnitudes across all three independence gaps. The effect is negligible.

---

## Conclusion

The independence threshold has essentially no impact on the sensitivity surface. This is expected: the deviation metric (Δλ) compares a sub-window and its full-year benchmark using the *same* independence gap, so any systematic effect of the gap on event counts cancels in the difference. The choice of 30 minutes is neither fragile nor consequential.

**Figure:** `figures/independence_threshold_sensitivity.pdf`

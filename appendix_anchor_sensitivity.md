# Appendix: Anchor Detection Parameter Sensitivity

Multi-year datasets are split into annual slices using the `find_anchors()` algorithm, which detects 12-month deployment "seasons" based on effort continuity. The algorithm depends on five parameters: minimum total effort (trap-days), seasonal balance (minimum proportion of effort in the smaller half-year), minimum detection days, maximum zero-detection gap, and maximum fraction of zero-detection days. These parameters determine which portions of a multi-year deployment qualify as usable annual slices.

We tested whether relaxing or tightening these parameters changes the main findings by re-running the pipeline with three configurations.

---

## Parameter configurations

| Parameter | Relaxed | Current (baseline) | Strict |
|-----------|---------|-------------------|--------|
| `min_effort` (trap-days) | 1,000 | 2,000 | 3,000 |
| `balance` (half-year proportion) | 0.10 | 0.20 | 0.30 |
| `detect_days` (minimum) | 180 | 250 | 300 |
| `zero_frac` (max fraction) | 0.60 | 0.50 | 0.40 |
| `zero_gap` (max consecutive days) | 90 | 60 | 45 |

---

## Anchor detection results

| Config | Slices detected | Datasets with slices |
|--------|-----------------|----------------------|
| Relaxed | 34 | 26 |
| **Current (baseline)** | **34** | **26** |
| Strict | 31 | 24 |

The relaxed configuration produced **identical** slice boundaries to the current baseline — none of the five parameters were binding in the relaxed direction. The strict configuration dropped 3 slices: NO-gravberget (insufficient total effort), SE-grimso-low (insufficient effort), and one NO-evenstadlia slice (seasonal balance violation). The strict pipeline was re-run from scratch (27.4 minutes).

---

## Model comparison

| Config | Species | Rows | Dev. expl. | Surface *r* vs baseline |
|--------|---------|------|------------|-------------------------|
| Relaxed | — | — | — | 1.000 (identical to baseline) |
| **Current (baseline)** | **29** | **222,748** | **87.1%** | **1.000** |
| Strict | 26 | 208,455 | 87.0% | 0.980 |

The strict configuration loses 3 species (those that only appeared in the dropped slices) and ~14,000 rows. Deviance explained is virtually unchanged (87.0% vs 87.1%). The predicted surface correlation of 0.980 confirms that the three dropped slices do not drive the overall surface pattern.

---

## Conclusion

The anchor detection algorithm is robust to parameter perturbation. The relaxed configuration is functionally identical to the baseline, and even strict parameters — which drop 3 of 34 slices — leave the surface shape and model performance essentially unchanged. The current default parameters are not a fragile choice point; they identify the same set of usable annual slices across a wide parameter range.

**Figure:** `figures/anchor_sensitivity.pdf`

# Appendix: Study Design Covariates

**Script:** `appendix_study_design.R`

**Status:** Preliminary — to be re-run once `dataset_meta.xlsx` is complete for all datasets and the final dataset is incorporated.

## Summary

Camera-trap study design parameters (placement strategy, grid spacing, camera height, trail placement) were tested as covariates in the sensitivity surface model. **None represent genuine causal effects on detection deviation.** All apparent signals trace to single-dataset leverage points or geographic confounding, and are properly absorbed by the dataset × species random effect. The sensitivity surface predictions are unaffected.

## Design covariates available

| Covariate | Levels / range | Contrast | Notes |
|-----------|---------------|----------|-------|
| Placement | random (63%), systematic (37%) | Adequate | "mixed" (FR-montpellier only) merged with random |
| Grid spacing | 0.07–3 km (median 1 km) | Adequate but clustered | 15 of 27 datasets at 1 km; BE-Leuven at 0.07 km is an extreme outlier |
| Camera height | 50–120 cm | Low | 77% at 50 cm |
| Trail use | off (93%), on (7%) | Minimal | Only 2 datasets with "on" (FR-montpellier, PL-kampinos) |
| Target species | all "multiple" | None | |
| Bait | all "no" | None | |

Target and bait have zero contrast and are excluded from all models.

**Imputed datasets (7):** BE-Leuven, GE-Mueritz_NP, NL-MICA_C1, NL-MICA_C2, NO-nina, PL-kampinos_NP, SP-donana. Imputed values are flagged and sensitivity to their removal is tested in § D.

## Results

### Design covariates add zero predictive power to the full model (§ A)

Adding all four design covariates to the M6 model (with dataset × species RE) yields ΔAIC ≈ 0 and identical deviance explained (87.0%). Design covariates are constant within datasets and therefore compete directly with the dataset component of the RE, which absorbs whatever signal they carry.

### Without the RE, design explains 1.6 pp of deviance (§ B)

| Model | Dev. expl. | ΔAIC vs baseline |
|-------|-----------|-----------------|
| Baseline (no design, no RE) | 69.6% | 0 |
| + spacing | 70.8% | −8,347 |
| + height | 69.9% | −2,257 |
| + placement | 69.9% | −1,644 |
| + trail | 69.6% | −86 |
| + all 4 | 71.2% | −12,376 |

Spacing dominates (~75% of total design signal). Trail contributes nothing.

### Design absorbs ~32% of dataset-level RE variance (§ C)

With a dataset-only RE (35 levels, no species component), design covariates reduce the RE standard deviation from 0.727 to 0.598 (32% variance reduction). Placement and spacing each account for ~25%. However, deviance explained and AIC are identical with or without design — the RE compensates perfectly.

### All effects are driven by single-dataset leverage (§ D)

**BE-Leuven drives spacing and height.** Its rotating single-camera design produces an extreme log(spacing) = −2.66 (next closest: −0.69) with 34,694 observations (16% of data).

| Covariate | β (full data) | β (excl. BE-Leuven) | β (verified only) |
|-----------|-------------|--------------------|--------------------|
| Placement (systematic) | +0.182 | +0.125 | +0.074 |
| Spacing (log, std) | **−0.221** | **+0.370** | **+0.122** |
| Height (std) | **−0.110** | **+0.022** | **+0.020** |
| Trail (on) | +0.131 | +0.105 | +0.026 |
| Latitude (std) | −0.327 | −0.425 | −0.327 |

Spacing and height **reverse sign** when BE-Leuven is removed.

**SI-serknica drives placement.** The placement coefficient flips from +0.074 to −0.067 (sign reversal) when SI-serknica is excluded from the verified-only model. No other individual dataset produces a sign change.

### Root cause: species-nonspecific l_nsites coefficient (§ E)

The M6 model applies a uniform `l_nsites` coefficient (β ≈ −1.35) to all species, meaning each additional camera reduces predicted deviation by the same factor. But the actual l_nsites–deviation relationship varies by species ecology:

| Species type | Example | r(l_nsites, |d_lambda|) at 57d |
|-------------|---------|-------------------------------|
| Flat-season | Meles meles | −0.78 |
| Flat-season | Sciurus vulgaris | −0.78 |
| Moderate-season | Vulpes vulpes | −0.62 |
| Strong-season | Sus scrofa | −0.46 |
| Strong-season | Cervus elaphus | −0.33 |

For species with strong seasonal signals (ungulate rut, hibernation), deviations are predominantly temporal — driven by the mismatch between a sub-window's timing and the annual pattern. More cameras do not reduce this temporal signal. The uniform l_nsites coefficient over-credits camera count for these species.

Systematic designs tend to deploy more cameras (by design for spatial coverage). At SI-serknica (50–111 cameras on a systematic grid), the model expects very low Cervus deviations based on camera count, but the actual autumn rut signal persists regardless. The model underpredicts, producing a positive residual that is incorrectly attributed to "systematic placement."

The correlation between per-dataset residual gradient (resid_183d − resid_15d) and l_nsites scaling (how much camera count increases from 15d to 183d windows) is **r = 0.66** across all 27 datasets. This confirms the mechanism is general: datasets where camera count scales more steeply with window length show a systematic residual pattern where the model over-predicts at short windows and under-predicts at long windows.

### Sensitivity surface is unaffected (§ F)

An `l_nsites × window_len` interaction term improves AIC by −643 (significant) but produces surface predictions with r > 0.9999 for all 29 species compared to the baseline M6 model. The tensor product surfaces already capture the duration effect, so the l_nsites scaling misspecification does not propagate to the primary output.

## Conclusions

1. **No study design parameter has a genuine causal effect** on detection deviation in this dataset. The spacing effect is a BE-Leuven leverage artifact (sign-flip on removal); the placement effect is a SI-serknica artifact (sign-flip on removal); height and trail lack adequate contrast.

2. **The underlying mechanism** is a species-nonspecific `l_nsites` coefficient that over-credits camera count for species with strong seasonal signals. This produces site-specific residual patterns that correlate with study design by proxy (systematic designs → more cameras), but the design itself is not the cause.

3. **The sensitivity surface predictions are robust.** Surface correlations between the baseline and corrected models exceed 0.9999 for all species. No correction is needed for the published results.

4. **This analysis should be re-run** once `dataset_meta.xlsx` is complete and the final dataset is added, to confirm that the null result holds with full metadata coverage and maximum sample size.

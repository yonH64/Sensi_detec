# Appendix: Structural Confounds Identified and Resolved

Three structural confounds were discovered during metric development. Each involved a detection metric whose deviation from the 12-month benchmark was partially or wholly driven by mathematical properties of the metric formula rather than genuine differences in detection. All three were resolved by replacing the confounded metric with a deconfounded alternative.

---

## 1. Time-to-event probability — mathematical scaling with window length

The time-to-event detection probability *p* = 1 − exp(−λ*L*) is mathematically guaranteed to be smaller for shorter windows, even if the underlying daily detection rate λ is identical between the sub-window and the full year. The deviation Δ*p* = *p*_window − *p*_full was 75–79% negative, driven entirely by this exponential scaling.

**Resolution.** We extracted the daily detection rate λ = −log(1 − *p*) / *L*, which inverts the exponential and recovers a window-length-independent quantity. The deviation Δλ = λ_window − λ_full is free of the mathematical confound and serves as the primary detection metric throughout the analysis.

---

## 2. Spatial coverage — two opposing confounds

Spatial coverage is the proportion of active cameras that detected a species. Its deviation from the full-year benchmark was affected by two confounds acting in opposite directions.

**Denominator confound.** The camera pool changes between windows. Sub-windows have fewer active cameras than the full year, so the full-year denominator is diluted by cameras that had no opportunity to detect the species during the sub-window period. This pushes the full-year spatial coverage down, making the sub-window appear relatively better.

**Accumulation confound.** More observation time means more cameras eventually detect the species. The full year benefits from 365 days of accumulation; a 60-day sub-window has less time for rare detection events. This pushes the full-year spatial coverage up relative to sub-windows.

These confounds partially cancel in the unmatched metric (76% of values were negative). Matching cameras — restricting to those active in both periods — fixes the denominator confound but exposes the accumulation confound fully (93% negative).

**Resolution.** Two-step deconfounding: (i) camera matching eliminates the denominator confound; (ii) rate extraction via −log(1 − spatial coverage) / *L*, analogous to the time-to-event transformation, removes the accumulation confound. The resulting matched detection rate deviation is free of both confounds.

---

## 3. Log encounter rate — floor amplification at low baseline rates

The log-transformed encounter rate deviation log(rate_window + ε) − log(rate_full + ε) amplified small absolute differences at low baseline rates. For the same absolute rate difference of 0.005 events/trap-day, the log deviation was ~5.5× larger at a baseline of 0.02 than at 0.10. Because northern sites have systematically lower encounter rates (median ~0.02 vs ~0.04 at southern sites), this created a spurious positive latitude coefficient — the opposite direction from the other two metrics, which both show smaller deviations at higher latitudes.

**Resolution.** Replaced by the raw encounter rate difference (rate_window − rate_full), which has no floor amplification and produces a latitude coefficient consistent with the other metrics.

---

## Summary

| Original metric | Confound | Replacement |
|-----------------|----------|-------------|
| Time-to-event probability deviation | Exponential scaling with *L* | Daily detection rate deviation (Δλ) |
| Spatial coverage deviation | Denominator + accumulation | Matched detection rate deviation |
| Log encounter rate deviation | Log amplification at low rates | Raw encounter rate deviation |

The analysis uses three deconfounded metrics: Δλ (TTE daily detection rate), Δ rate (raw encounter rate), and Δ matched rate (camera-matched, accumulation-corrected detection rate).

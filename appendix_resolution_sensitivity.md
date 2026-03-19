# Appendix: Window-Start Resolution Sensitivity

The sliding window grid places start positions at regular intervals across the year. The baseline uses a 7-day step (53 start positions × 25 durations = 1,325 windows per dataset-slice). This choice balances computational cost against granularity. We tested whether using a finer (3-day) or coarser (14-day) step changes the main findings.

The 3-day step required a full pipeline re-run (122 start positions × 25 durations = 3,055 windows per dataset-slice; 69.7 minutes). The 14-day step was obtained by subsetting every other start position from the 7-day baseline data. For each resolution, we re-prepared the species-level data and re-fitted the primary model (M6: species-specific 2D tensor product surfaces for |Δλ|, Gamma family with log link).

---

## Data retention

| Resolution | Start positions | Windows per slice | Species | Rows |
|------------|-----------------|-------------------|---------|------|
| 3-day step | 122 | 3,055 | 28 | 512,990 |
| **7-day step (baseline)** | **53** | **1,325** | **29** | **222,748** |
| 14-day step | 26 | 654 | 29 | 109,298 |

The 3-day resolution more than doubles the data volume. One species (*Ovis aries*) drops out because at the finer resolution the minimum-rows-per-species threshold is not met in all datasets. The 14-day resolution halves the data but retains all 29 species.

---

## Model comparison

| Resolution | Dev. expl. | Surface *r* vs baseline |
|------------|------------|-------------------------|
| 3-day step | 87.1% | 0.979 |
| **7-day step (baseline)** | **87.1%** | **1.000** |
| 14-day step | 87.1% | 0.997 |

Deviance explained is identical (87.1%) across all three resolutions. The predicted surface correlation with the 7-day baseline exceeds 0.97 in both directions: finer sampling adds no information, and coarser sampling loses negligible structure. The modest reduction at 3-day (r = 0.979 vs 0.997 for 14-day) likely reflects slight changes in the species pool rather than temporal resolution itself.

---

## Conclusion

The 7-day step provides more than sufficient temporal resolution. The sensitivity surface shape is virtually identical whether windows are placed every 3, 7, or 14 days. The main features — the autumn deviation hotspot, the duration-smoothing curve, and the species-specific surface shapes — are fully captured at 7-day resolution.

**Figure:** `figures/resolution_sensitivity.pdf`

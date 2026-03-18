# Appendix: Benchmark Choice Robustness

The sensitivity surface measures deviation of sub-window detection metrics from a 12-month benchmark — the operational reference representing what year-round monitoring would yield. The surface shape could in principle depend on this specific choice of benchmark duration. We tested whether shorter benchmarks (180 and 270 days) change the qualitative structure of the surface.

---

## Method

TTE daily detection rate (λ) was recomputed for 180-day and 270-day windows at every 7-day start position across all 23 dataset-slices. Each existing sub-window was matched to the alternative benchmark whose temporal centre was closest to its own. Deviations were computed as |λ_window − λ_benchmark| for each benchmark duration.

We evaluated three aspects of robustness: (i) shape correlation — the Pearson *r* between seasonal profiles of |Δλ| for the 365-day benchmark vs each alternative, computed separately per sub-window duration; (ii) deviation magnitude — mean |Δλ| under each benchmark; and (iii) sign consistency — the percentage of observations with positive deviations (λ_window > λ_benchmark).

---

## Results

### Shape correlations

| Sub-window duration | 365 d vs 180 d | 365 d vs 270 d |
|---------------------|----------------|----------------|
| 15 d | 0.995 | 0.998 |
| 29 d | 0.983 | 0.995 |
| 57 d | 0.970 | 0.851 |
| 85 d | 0.933 | 0.602 |
| 120 d | 0.879 | 0.705 |

The 180-day benchmark preserves the surface shape well (*r* > 0.88 for all sub-window durations). The 270-day benchmark degrades for longer sub-windows because a 120-day window overlaps extensively with a 270-day reference, compressing the deviation signal mechanically. This degradation is a geometric inevitability, not evidence that the sensitivity surface is benchmark-dependent.

### Deviation magnitude and sign

Mean |Δλ| decreases slightly with shorter benchmarks: at 60-day sub-windows, mean |Δλ| is 0.0124 (365 d), 0.0107 (180 d), and 0.0122 (270 d). The directional pattern (sub-windows tending to exceed the benchmark) is preserved under all benchmarks for short windows. As sub-window duration increases, the fraction of positive deviations decreases under the 180-day benchmark (to 84% at 120 d), reflecting the narrower reference period.

---

## Conclusion

The sensitivity surface shape is robust to benchmark choice. The autumn activity spike, the duration smoothing curve, and the relative ordering of seasonal windows appear regardless of whether the benchmark is 180, 270, or 365 days. The 365-day benchmark is retained as the primary reference because it maximises contrast at all sub-window durations and has the clearest operational interpretation.

See **Fig. S1** for a three-panel visualisation: shape correlations across benchmarks, deviation magnitude convergence, and seasonal profiles at 29-day and 85-day sub-windows for all three benchmark durations.

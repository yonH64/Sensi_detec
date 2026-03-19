# Appendix: Joint Species-Inclusion Threshold Sensitivity

The individual threshold tests (Appendix: Species-Inclusion Threshold Robustness) varied one threshold at a time while holding the others at baseline. Here we test all three thresholds simultaneously — a more conservative test that captures any interactions between threshold choices.

For each configuration, we re-ran the data preparation pipeline with all three thresholds set jointly, then re-fitted the primary model (M6).

---

## Configurations

| Config | min_events | min_sites_pos | min_occasions_pos |
|--------|-----------|---------------|-------------------|
| Lenient | 10 | 3 | 3 |
| **Baseline** | **20** | **5** | **5** |
| Strict | 30 | 10 | 5 |

---

## Data retention

| Config | Species | Rows | Change vs baseline |
|--------|---------|------|--------------------|
| Lenient (10/3/3) | 30 | 276,986 | +54,337 rows, +4 species |
| **Baseline (20/5/5)** | **26** | **222,649** | — |
| Strict (30/10/5) | 24 | 183,511 | −39,138 rows, −2 species |

Note: The baseline row count here (222,649) differs slightly from the pipeline baseline (222,748) because the joint test applies the threshold filter post-hoc to the same raw data, rather than letting the pipeline apply them during metric computation. The 99-row difference reflects edge cases in the filtering sequence and has no practical impact.

---

## Model comparison

| Config | Dev. expl. | Surface *r* vs baseline |
|--------|------------|-------------------------|
| Lenient (10/3/3) | 87.2% | 0.976 |
| **Baseline (20/5/5)** | **87.1%** | **1.000** |
| Strict (30/10/5) | 87.1% | 0.970 |

Deviance explained is stable across all three configurations (87.1–87.2%). Surface correlations exceed 0.97 in both directions. The joint variation of all three thresholds simultaneously produces no greater instability than varying them individually — consistent with the `min_occasions_pos` threshold being non-binding (see Appendix: Species-Inclusion Threshold Robustness, §C).

---

## Conclusion

Simultaneous variation of all three species-inclusion thresholds confirms the individual-threshold findings: the sensitivity surface is robust. The 4 additional species gained under lenient thresholds and the 2 species lost under strict thresholds have negligible impact on the predicted surface shape. The baseline thresholds (20/5/5) are a defensible middle ground.

**Figure:** `figures/threshold_sensitivity_joint.pdf`  
**Output data:** `threshold_sensitivity_joint_summary.csv`

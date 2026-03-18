# Appendix: Spearman's ρ Boundary Inflation and Species Count Filter

## Background

Spearman's rank correlation coefficient (ρ) was computed between species-level detection rates in each sub-window and the corresponding 12-month benchmark, using species detected in both periods. This metric quantifies rank preservation — whether a sub-window maintains the relative detectability ordering of species, even if absolute rates are biased.

## Problem

With very few shared species, Spearman's ρ becomes unreliable. With 3 species, ρ can take only five discrete values (−1, −0.5, 0, 0.5, 1), and perfect agreement (ρ = 1) occurs frequently even under moderate noise. With 4 species the resolution improves but remains coarse. In our dataset, 33% of all ρ values were exactly 1, and 73% of these perfect correlations came from windows sharing only 3–4 species with the full-year benchmark. This boundary mass is artefactual — it reflects the discreteness of the statistic at low species counts rather than genuinely perfect rank preservation.

## Solution and effect

We modelled ρ on the transformed scale (ρ + 1) / 2 using beta regression (logit link), which assumes a continuous response on (0, 1) and is not designed to accommodate point masses at the boundaries. Including all observations yielded very low deviance explained, with the boundary mass at ρ = 1 dominating the residuals.

Restricting to windows sharing **≥ 5 species** with the benchmark removed the observations predominantly driving the artefact. The threshold of 5 was chosen pragmatically: it is the minimum at which Spearman's ρ has enough resolution to be meaningfully continuous. A threshold of 6 or 7 gives similar results but excludes more data; a threshold of 4 retains too much boundary inflation.

After filtering, the dataset retained 12,028 observations with a deviance explained of 40.4%.

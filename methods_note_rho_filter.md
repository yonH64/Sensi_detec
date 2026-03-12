# Methods note: Filtering Spearman's rho to windows with ≥ 5 shared species

## Rationale

Spearman's rank correlation coefficient (rho) was computed between species-level
detection rates in each sub-window and the corresponding 12-month benchmark, using
species detected in both periods (the "shared species"). In windows with very few
shared species, rho becomes unreliable as a measure of rank preservation:

- With 3 species, Spearman's rho can take only a small number of discrete values
  (−1, −0.5, 0, 0.5, 1), and perfect rank agreement (rho = 1) occurs frequently
  even under moderate noise.
- With 4 species, the number of possible values increases but remains coarse.

In our dataset, 33% of all rho values were exactly 1. Of these, 73% came from
windows sharing only 3–4 species with the full-year benchmark. This boundary mass
is artefactual — it reflects the coarseness of the rank-correlation statistic at
low species counts rather than genuinely perfect rank preservation.

## Effect on model fit

We modelled rho on the transformed scale ((rho + 1) / 2) using beta regression
(betar family with logit link). Beta regression assumes a continuous response on
(0, 1) and is not designed to accommodate point masses at the boundaries.

Including all observations yielded 17.3% deviance explained. Restricting to
windows with ≥ 5 shared species (removing 38% of observations, predominantly
those driving the boundary artefact) more than doubled deviance explained to 37.1%.

## Suggested methods text

> For the rank-preservation model (Spearman's rho between sub-window and full-year
> detection rates), we restricted the analysis to dataset × window combinations
> sharing at least five species with the 12-month benchmark. With fewer than five
> shared species, Spearman's rho is limited to a small number of discrete values
> and frequently produces artefactual perfect correlations (rho = 1); in our
> dataset, 73% of perfect correlations arose from windows with only 3–4 shared
> species. This filter retained 9,887 of the original 15,891 observations (62%)
> and substantially improved model fit (deviance explained: 37.1% vs 17.3%
> without filtering).

# Appendix: Model Comparison (M1–M6)

Eight model variants were compared for the primary detection metric (|Δλ|, Gamma family with log link). All models share the same fixed covariates (log trap-days, log number of sites, latitude, trap array size) and random intercepts (dataset × species). They differ in how the 2D tensor product surface (day-of-year × duration) varies across species or ecological guilds.

---

## Model structures

| Model | Surface structure | Rationale |
|-------|-------------------|-----------|
| M1 | Single shared surface | Baseline: same surface for all species |
| M2 | Shared surface + guild-varying seasonal spline | Guilds differ in seasonal amplitude |
| M3 | Full guild × surface (5 major guilds) | Each guild gets its own 2D surface |
| M4 | Shared surface + environment × duration interaction | Environmental modulation of the duration effect |
| M5 | Guild surface + environment × duration | Combines M3 and M4 |
| M_hab | Habitat guild surfaces (9 levels) | Finer grouping by habitat affinity |
| M_diet | Diet guild surfaces (10 levels) | Finer grouping by foraging ecology |
| **M6** | **Species-specific surfaces (26 species)** | **Each species gets its own 2D surface** |

All models use cyclic splines for day-of-year and tensor product smooths for the timing × duration interaction. Penalisation automatically regularises species or guilds with fewer observations.

---

## Results

| Model | AIC | Dev. expl. | edf | ΔAIC vs M6 |
|-------|-----|------------|-----|------------|
| **M6 (species)** | **−1,123,832** | **86.8%** | **839** | **0** |
| M_diet | −1,113,208 | 85.6% | 509 | +10,625 |
| M_hab | −1,110,209 | 85.2% | 469 | +13,624 |
| M5 (guild + env) | −1,109,765 | 85.2% | 486 | +14,067 |
| M3 (guild surface) | −1,108,433 | 85.0% | 483 | +15,399 |
| M2 (guild season) | −1,107,520 | 84.9% | 329 | +16,312 |
| M4 (env × duration) | −1,100,257 | 84.0% | 319 | +23,575 |
| M1 (base) | −1,099,062 | 83.9% | 315 | +24,770 |

---

## Interpretation

Species-specific surfaces (M6) are clearly preferred (ΔAIC = −10,625 vs the next-best model). Individual species have sufficiently distinct timing × duration profiles that neither major guild (5 levels) nor minor guild (9–10 levels) captures the variation adequately.

The progression from M1 through M3 shows incremental improvement from adding guild-level surface variation. The jump from guilds to species (M3 → M6) is the single largest improvement, consistent with species identity being the dominant source of surface heterogeneity. Among guild-level models, diet guilds (M_diet) outperform habitat guilds (M_hab), suggesting that foraging ecology — with its implications for activity patterns and detection seasonality — is a more relevant grouping than habitat affinity.

Despite the large AIC differences, deviance explained spans only 3 percentage points (83.9–86.8%). The shared covariates and random effects explain the majority of variance; the surface structure refinement captures the remaining species-level heterogeneity.

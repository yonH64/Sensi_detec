# Vignette: Chamois altitudinal migration as a case study in site-level detection idiosyncrasy

## Summary

Alpine chamois (*Rupicapra rupicapra*) at GE-Berchtesgaden National Park exhibits the largest species-level prediction error in our leave-one-dataset-out cross-validation — not because the model fails structurally, but because Berchtesgaden chamois have a seasonal detection profile fundamentally different from the same species at other sites. This difference is driven by **altitudinal migration**, an ecological process that creates site-level variation the fixed-effect structure cannot anticipate. The case illustrates why the dataset × species random effect is a necessary model component, and why out-of-sample prediction of deviation *magnitude* for a new site is inherently limited.

## Background

The LOO-CV for GE-Berchtesgaden_NP produced an overall predictive correlation of r = 0.10, by far the lowest of 27 base datasets. Initial investigation showed that per-species surface *shapes* were well-captured (within-species r = 0.60–0.94), but deviation *magnitudes* were systematically wrong. Four of five species had deviations 3–11× smaller than the cross-site average (model over-predicted), while *Rupicapra rupicapra* had deviations 9.5× larger (model under-predicted). These opposing biases cancel out in the overall correlation.

## The chamois detection profile at Berchtesgaden

Short-window (15–30 day) lambda values for chamois at GE-Berchtesgaden_NP reveal a pronounced **spring peak** in detection:

| Month | Lambda (daily detection rate) |
|-------|------------------------------|
| Jan   | 0.020 |
| Feb   | 0.025 |
| Mar   | 0.037 |
| **Apr** | **0.057** |
| **May** | **0.049** |
| Jun   | 0.039 |
| Jul   | 0.032 |
| Aug   | 0.029 |
| Sep   | 0.018 |
| **Oct** | **0.016** |
| Nov   | 0.033 |
| Dec   | 0.029 |

The seasonal swing is 3.5× (April peak / October trough). The dominant signal is a **spring peak**, not an autumn rut peak — the opposite timing from the ungulate rut that dominates other species in the model. A secondary November bump is consistent with the actual chamois rut (which peaks in November), but it is smaller than the spring migration signal.

## Ecological interpretation: altitudinal migration

Alpine chamois are well-documented altitudinal migrants. In late autumn and winter, deep snow at high elevations drives chamois to descend to forested slopes and valley floors. By spring (March–May), animals are concentrated at lower elevations — exactly where camera traps are typically deployed in mountainous terrain. As snow melts and alpine pastures become accessible in summer, chamois disperse upward, reducing encounter rates at camera stations.

Berchtesgaden National Park (47.6°N, Bavarian Alps) has substantial vertical relief, with camera traps placed in valley and forest zones. The spring detection peak reflects chamois moving through and concentrating in the camera-trap elevation band during the transition from winter range to summer alpine pastures.

## Contrast with other sites

*Rupicapra rupicapra* appears in three base datasets:

| Site | Context | Mean lambda | Seasonal swing | Mean |d_lambda| |
|------|---------|-------------|----------------|-------------------|
| GE-Berchtesgaden_NP | Bavarian Alps, high relief | 0.018 | 3.5× (spring peak) | 0.0103 |
| IT-Alps | Italian Alps, high relief | 0.0016 | ~3.5× (spring signal in long windows) | 0.0009 |
| SI-serknica | Slovenian karst, low relief | 0.0018 | ~1.3× (near-flat) | 0.0011 |

**IT-Alps** shows the same spring-biased seasonal pattern — windows starting in February–March have lambda 3.6× higher than those starting in August–September (0.0039 vs 0.0010), even though only long windows (113–141 days) are available for this site. However, the baseline detection rate is 11× lower than Berchtesgaden, so absolute deviations remain small.

**SI-serknica** (Slovenian karst plateau, lower elevation, less vertical relief) shows a nearly flat detection profile (lambda ≈ 0.005 year-round). Without significant altitudinal range, there is no migration-driven concentration of animals at camera elevation. Chamois are present but uncommon and temporally stable.

The key insight is that **the same ecological process (altitudinal migration) produces the seasonal pattern at both Alpine sites**, but the magnitude of its effect on deviations depends on local population density and camera placement relative to the altitudinal range. Both factors are site-specific and cannot be predicted from the fixed effects (species identity, latitude, effort).

## Why the model fails at Berchtesgaden

The M6 model learns *Rupicapra rupicapra*'s surface primarily from the Slovenian and Italian data (3 dataset-slices). The Slovenian surface is near-flat; the Italian data is sparse and only available at long window durations. The resulting species-specific tensor product surface reflects a species with modest, autumn-centred deviations — not the spring-dominated, high-amplitude signal at Berchtesgaden.

Without the random effect, the model predicts Berchtesgaden chamois deviations at approximately the Slovenian level (mean predicted: 0.0017 vs observed: 0.0103). The RE, with SD = 0.94 on the log link, is designed to absorb exactly this kind of site-level shift (exp(0.94) ≈ 2.6× per SD, so a 6× discrepancy is within ~2 SD).

Conversely, the four other Berchtesgaden species (*Capreolus capreolus*, *Cervus elaphus*, *Lepus timidus*, *Vulpes vulpes*) all have deviations 3–11× *lower* than the cross-site average. This is consistent with compressed seasonal detection variation at this high-elevation site — the same mountainous topography that concentrates chamois at camera elevation may disperse or stabilise detection of other species that are less strongly migratory.

## Broader implications

This case study illustrates three points relevant to the sensitivity surface analysis:

1. **The random effect is not nuisance variance.** The dataset × species RE (SD = 0.94) captures ecologically meaningful site-level variation — in this case, whether a montane species undergoes altitudinal migration through the camera-trap zone. This variation is real but unpredictable from the fixed effects.

2. **Out-of-sample magnitude prediction has a hard ceiling.** The LOO-CV R² of ~31% (vs 87% in-sample) reflects the RE's share of explained variance. For a genuinely new site, the model can predict the *shape* of the sensitivity surface (which seasons and durations are more or less problematic) but not the *absolute level* of deviation without site-specific calibration data.

3. **The same species can require different monitoring designs at different sites.** A chamois monitoring protocol at Berchtesgaden would need to account for the spring migration peak (avoid March–May windows if unbiased lambda is the goal), while at a Slovenian karst site, timing matters much less for this species. Species-specific recommendations that ignore topographic context may mislead.

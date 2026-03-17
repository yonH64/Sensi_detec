# Research Question Framing: What Are We Actually Measuring?

**Discussion document — March 2026**

---

## The core issue

The current study compares detection metrics from sub-annual sampling windows to a 12-month benchmark (the "FULL" window). The deviation surface tells us: *how much does a short window's estimate differ from what year-round monitoring would give?*

This is a valid question, but it is not the only valid question — and the choice of benchmark implicitly defines what counts as "error." The concern is that the study conflates two distinct things:

1. **Seasonal variation in detection** — a biological reality (e.g., the ungulate rut genuinely increases encounter rates in autumn)
2. **Estimation error from temporal subsampling** — a methodological problem

The annual-average benchmark treats *all* seasonal variation as if it were error. But for many applications, seasonal detection variation is signal, not noise.

Evidence of this conflation: ~99% of d_lambda values are positive. Sub-windows almost always yield higher daily detection rates than the annual average. This isn't because sub-windows are "biased" — it's because the annual mean is diluted by low-activity months. The sensitivity surface is essentially a map of seasonal amplitude, reframed as deviation.

---

## What question does a practitioner actually need answered?

Different research goals imply different benchmarks and different notions of what makes a "good" sampling window:

### Q1: "How well does this window approximate year-round monitoring?"

- **This is the current framing.** Benchmark = 12-month mean.
- Relevant when: the goal is to estimate an annual-average detection rate (e.g., for standardised monitoring indices).
- Limitation: penalises windows that capture high-activity periods, even though those windows may be operationally superior for other purposes.

### Q2: "Does this window give enough detection probability to estimate occupancy?"

- Benchmark = not a rate, but a precision threshold (SE or CI width of the occupancy estimate).
- A 60-day autumn window that "overestimates" the annual rate is actually *better* for occupancy — higher detection probability means fewer false absences and tighter occupancy estimates.
- The current framing penalises exactly the windows that are best for occupancy modelling.

### Q3: "Can I detect population trends using this window?"

- Benchmark = between-year consistency of the same window, not comparison to an annual mean.
- A window that systematically overestimates by 30% every year is fine for trend detection — the bias cancels in the year-to-year contrast. What matters is low inter-annual variance *within* the window.
- This reframes the question from "how biased is this window?" to "how repeatable is this window?"

### Q4: "Does this window preserve the relative detectability of species?"

- Already addressed by rho_lambda (Spearman rank correlation). Result: ~0.96 at 64 days across the surface (model restricted to windows with ≥5 shared species).
- This is arguably the most practically important finding. For multi-species studies, if rank order is preserved, the window is adequate for comparative analyses even if absolute rates are biased.
- Currently treated as supplementary.

### Q5: "Does this window capture the full community?"

- Addressed by richness models (prop_sr_full, d_sr_raref).
- Relevant for biodiversity inventories. Benchmark = full-year species list.
- Already in the study; framing is appropriate.

### Q6: "What is the detection rate during this season specifically?"

- No benchmark needed. The rate *is* the quantity of interest.
- Relevant when: studying seasonal ecology (e.g., autumn ungulate behaviour, winter predator activity).
- The annual-average comparison answers a question nobody asked in this context.

---

## What this means for the paper

### Option A: Keep the current framing, but make it explicit

State clearly that the 12-month benchmark represents "what year-round monitoring would yield" — an operational reference, not ecological truth. Acknowledge that the deviation surface primarily captures seasonal detection amplitude. Add a paragraph in the Discussion noting that high-deviation windows are not necessarily "bad" windows — they are windows where detection rates differ most from the annual average, which may be desirable for occupancy estimation and undesirable for standardised monitoring indices.

**Pros:** Minimal rework. Honest. The surface is still informative.
**Cons:** The "sensitivity" framing still implies deviation = bad, which readers may internalise.

### Option B: Reframe as "seasonal detection variation" with application-specific guidance

Present the surface as a characterisation of *how much detection rates vary* as a function of window timing, duration, species, and environment. Then provide application-specific interpretation:

- For standardised monitoring (annual index estimation): deviation from annual mean matters → use long or "neutral" windows.
- For occupancy estimation: high detection probability matters → autumn windows may be preferable despite (or because of) high seasonal rates.
- For trend detection: inter-annual consistency matters → show that bias is stable across years (you have multi-year slices to test this).
- For community comparisons: rank preservation matters → rho_lambda is high everywhere, so most windows are adequate.

**Pros:** More nuanced. More useful to practitioners. Differentiates the paper from a simple "deviation = bad" message.
**Cons:** More complex framing. Requires some additional analyses (e.g., occupancy precision simulation, inter-annual consistency check).

### Option C: Dual-target approach

Model both (a) deviation from the annual mean (current approach) and (b) estimation precision (via MSE or SE directly). The first tells you about seasonal amplitude; the second tells you about statistical quality. Short high-activity windows may score poorly on (a) but well on (b).

**Pros:** Directly addresses the conflation. MSE models are straightforward to add (same M6 structure, data already computed).
**Cons:** More models, more results to present.

---

## The benchmark validity question

Separate from framing, the 12-month benchmark has some technical vulnerabilities:

- **Inter-annual variability:** Among sliced datasets, the 12-month lambda_full has a median CV of ~13% across years (up to 87% for rare species at individual sites). The signal-to-noise ratio is decent for short windows (median deviation/inter-annual SD ≈ 14:1), but for long windows (90–120 days) where deviations shrink, benchmark noise becomes a non-trivial fraction of the signal.
- **The 365-day cutoff is arbitrary:** Some datasets have >2 years of data, but the anchor algorithm enforces exactly 365-day windows. There's no principled reason 365 is the "right" reference duration.
- **Benchmark uncertainty is not propagated:** Each lambda_full is treated as known. In reality, it has its own SE (lambda_se_full), which could be incorporated into the MSE decomposition or as measurement error in the model.

A robustness check — does the surface shape change if you use, say, a 6-month or 9-month benchmark? — would help establish whether the findings are benchmark-dependent or reflect genuine structure.

---

## Missing: MSE models

The MSE columns (bias² + variance) are computed in the data pipeline but never modeled. This is a gap:

- MSE is the proper loss function for estimation quality.
- Currently only |bias| is modeled.
- Variance contributes 6–17% of MSE depending on window length (rising for longer windows as bias shrinks).
- For practical window design recommendations, a practitioner needs to know the total estimation error, not just bias.
- Adding MSE models with the same M6 structure would be straightforward and would address the "what about precision?" question directly.

---

## Suggested discussion points

1. Which framing (A, B, or C) best serves the paper's contribution?
2. Is the 12-month benchmark sufficiently justified, or do we need a robustness check?
3. Should MSE models be added as primary or supplementary results?
4. Should the rank-preservation result (rho_lambda ≈ 0.96 at typical durations) be elevated to a main finding?
5. Do we need application-specific guidance (occupancy, trends, inventories) in the Discussion, or is that out of scope?

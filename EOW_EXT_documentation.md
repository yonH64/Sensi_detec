# EOW_EXT: Extended Observation Window

> **⚠️ SUPERSEDED** — EOW_EXT has been replaced in the pipeline by **EOW_EARLY** (Aug 2–Sep 30, 60d) and **EOW_LATE** (Oct 1–Nov 29, 60d), which split the autumn window at the CORE midpoint (Oct 1) to provide a pre-rut vs post-rut contrast at matching duration. EOW_EXT is no longer computed by `Full1.R`. This document is retained for historical reference only.

## Overview

**EOW_EXT** (Extended Observation Window) was a 6-week (42-day) sampling window centered on the **Snapshot Europe Core** protocol period. It provided a focused alternative to the full SNAP_EU_CORE (61 days) and SNAP_EU_BUFFER (89 days) windows.

---

## Window Definitions

| Window | Start | End | Length | Description |
|--------|-------|-----|--------|-------------|
| **SNAP_EU_CORE** | Sep 1 (DOY 244) | Oct 31 (DOY 304) | 61 days | Snapshot Europe Core protocol |
| **SNAP_EU_BUFFER** | Aug 18 (DOY 230) | Nov 14 (DOY 318) | 89 days | Extended protocol with 2-week buffer |
| **EOW_EXT** | Sep 10 (DOY 253) | Oct 21 (DOY 294) | 42 days (6 weeks) | **New: Centered 6-week window** |

---

## Design Rationale

### Why 6 weeks?

1. **Practical duration**: 6 weeks (42 days) is:
   - Long enough to capture meaningful detection patterns
   - Short enough to be feasible for many study designs
   - A standard field season duration in ecology

2. **Centered on peak activity**: The window is centered on the midpoint of SNAP_EU_CORE:
   - **Midpoint**: DOY 274 (Oct 1)
   - Captures the peak autumn activity period for many species
   - Balances early September and mid-October conditions

3. **Comparison to existing protocols**:
   - **31% shorter** than SNAP_EU_CORE (42 vs 61 days)
   - **53% shorter** than SNAP_EU_BUFFER (42 vs 89 days)
   - Allows testing of shorter standardized protocols

---

## Window Overlap

```
Aug        Sep              Oct              Nov
|----------|----------------|----------------|----------|
           CORE (61d)
           ==================
    BUFFER (89d)
    ==========================================
                EOW_EXT (42d)
                ===============
           |        ↓        |
         Sep 1   Oct 1    Oct 31
                (midpoint)
```

- EOW_EXT starts **9 days after** SNAP_EU_CORE
- EOW_EXT ends **10 days before** SNAP_EU_CORE ends
- EOW_EXT is **fully contained** within both CORE and BUFFER

---

## Implementation in Code

The window is defined in `Full1.R` in the `add_snapshot_europe_windows()` function:

```r
# SNAP_EU_CORE reference
core_start <- doy(9, 1)    # Sep 1 = DOY 244
core_end   <- doy(10, 31)  # Oct 31 = DOY 304
core_len   <- 61

# Calculate EOW_EXT centered on SNAP_EU_CORE midpoint
core_midpoint <- core_start + floor(core_len / 2)  # DOY 274 (Oct 1)
eow_len       <- 42L                                # 6 weeks
eow_start     <- core_midpoint - floor(eow_len / 2) # DOY 253 (Sep 10)
eow_end       <- eow_start + eow_len - 1L           # DOY 294 (Oct 21)
```

The window is automatically added to `window_grid` and processed by `dataset_wrapper1()` like any other window.

---

## Usage in Analysis

### Access EOW_EXT results

```r
# Filter for EOW_EXT windows
eow_results <- all_window_species %>%
  filter(window_id == "EOW_EXT")

# Compare to other protocol windows
protocol_comparison <- all_window_species %>%
  filter(window_id %in% c("SNAP_EU_CORE", "SNAP_EU_BUFFER", "EOW_EXT"))
```

### Visualize EOW_EXT

```r
# Calendar heatmap showing EOW_EXT
plot_calendar_heatmaps(
  all_window_species,
  dropped_df = all_dropped_species,
  merge_datasets = TRUE,
  view = "offset"
)

# EOW_EXT will appear with a blue outline (protocol window marker)
```

### Compare metrics

```r
# Compare detection metrics across protocol windows
all_window_species %>%
  filter(
    window_id %in% c("SNAP_EU_CORE", "SNAP_EU_BUFFER", "EOW_EXT"),
    species == "Capreolus capreolus"
  ) %>%
  select(dataset, window_id, lambda, log_rate, matched_rate)
```

---

## Research Questions

EOW_EXT enables investigation of:

1. **Optimal protocol duration**:
   - Is 6 weeks sufficient for reliable species detection?
   - How do metrics (lambda, matched_rate) compare to 9-week CORE?

2. **Centered vs full protocol**:
   - Does centering on peak activity improve detection?
   - What species benefit most from the focused timing?

3. **Effort-efficiency tradeoff**:
   - How much detection probability is lost with a shorter window?
   - Does spatial coverage differ between EOW_EXT and CORE?

4. **Seasonal patterns**:
   - Is the centered 6-week period more stable across datasets?
   - How does EOW_EXT perform in early vs late autumn conditions?

---

## Expected Results

### Hypothesis: EOW_EXT performance

**Common species** (e.g., Capreolus capreolus, Sus scrofa):
- `matched_rate` ≈ SNAP_EU_CORE (detection rate stable)
- `lambda` slightly lower (less time to detect)
- `d_matched_rate` ≈ 0 (minimal loss of detection quality)

**Rare species** (low abundance):
- `matched_rate` < SNAP_EU_CORE (lower detection rate)
- `lambda` moderately lower
- `d_matched_rate` < 0 (reduced detection quality)

**Seasonal migrants** (late arrivals):
- Performance depends on arrival timing
- May perform poorly if peak activity is late October

---

## Comparison with Other Windows

| Metric | SNAP_EU_CORE | SNAP_EU_BUFFER | EOW_EXT | Expected Δ |
|--------|--------------|----------------|---------|------------|
| Duration | 61 days | 89 days | 42 days | -31% vs CORE |
| Trap-days | ~3500* | ~5000* | ~2400* | -31% vs CORE |
| Species richness | High | Highest | Moderate | -5 to -10% |
| Detection rate (lambda) | High | Highest | Moderate | slight decrease |
| Matched rate | High | High | Moderate | slight decrease |

*Assuming ~60 cameras

---

## Recommendations

### When to use EOW_EXT

✅ **Use EOW_EXT when:**
- You want a standardized protocol shorter than 2 months
- Budget/effort constraints limit field season duration
- Peak activity period is well-defined (autumn)
- Focus is on common species

❌ **Avoid EOW_EXT when:**
- Rare species are the primary target
- Year-round activity patterns are of interest
- Maximum spatial coverage is critical
- Arrival/departure timing is uncertain

---

## Next Steps

To validate EOW_EXT performance:

1. **Run analysis** with updated `window_grid` including EOW_EXT
2. **Compare MSE** between EOW_EXT and SNAP_EU_CORE
3. **Test across datasets** to assess generalizability
4. **Examine species-specific patterns** to identify which benefit most

---

## References

- Snapshot Europe protocol: https://www.mammalweb.org/snapshot-europe/
- This analysis: `Full1.R`, `dataset_wrapper1()`, `add_snapshot_europe_windows()`

---

## Changelog

- **2026-02-19**: Created EOW_EXT (6-week window centered on SNAP_EU_CORE)
- Initial definition: DOY 253-294 (Sep 10 - Oct 21)

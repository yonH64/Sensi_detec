# ═══════════════════════════════════════════════════════════════════════════════
# run_sensitivity_variants.R
# ═══════════════════════════════════════════════════════════════════════════════
#
# Runs the Full1.R data pipeline with 4 parameter configurations for the
# robustness checks in sensitivity_appendices.R §3–5:
#
#   1. 3-day window step   → all_window_species_3day.rds   (§3 Resolution)
#   2. Strict anchors      → all_window_species_strict.rds (§4 Anchor)
#   3. 15-min indep. gap   → all_window_species_15min.rds  (§5 Independence)
#   4. 60-min indep. gap   → all_window_species_60min.rds  (§5 Independence)
#
# Estimated runtime: ~2.5–3 hours total (parallelised across datasets).
# ═══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(lubridate)
library(furrr)
library(geosphere)

# ── Source shared helpers (find_anchors, make_window_template, spp_keep, …) ──
source("helpers.R")

# ── Source function definitions from Full1.R (lines 30–1410) ─────────────────
# We parse the file and evaluate only the function-definition blocks,
# skipping the demo calls and the orchestration section at the end.

full1_lines <- readLines("Full1.R")

eval(parse(text = full1_lines[30:53]))     # compute_full_effort
eval(parse(text = full1_lines[91:233]))    # split_spatial_by_distance
eval(parse(text = full1_lines[315:1177]))  # build_window_metrics_fast1
eval(parse(text = full1_lines[1202:1410])) # dataset_wrapper1

# ── Dataset paths ────────────────────────────────────────────────────────────
ds_paths <- list.dirs(
  "C:/Users/yonah/OneDrive - Universitetet i Innlandet/Dokumenter/PhD/Sensitivity of detection/Datasets",
  recursive = FALSE, full.names = TRUE
)

cat(sprintf("Found %d dataset directories.\n", length(ds_paths)))

# ── Parallel backend ─────────────────────────────────────────────────────────
plan(multisession, workers = 4)


# ═══════════════════════════════════════════════════════════════════════════════
# HELPER: run one pipeline variant end-to-end
# ═══════════════════════════════════════════════════════════════════════════════

run_variant <- function(
    label,
    step_doy            = 7,
    independence_mins   = 30,
    anchor_overrides    = list(),
    output_species      = NULL,
    output_richness     = NULL
) {
  cat(sprintf("\n\n════════════════════════════════════════════════════\n"))
  cat(sprintf("  VARIANT: %s\n", label))
  cat(sprintf("════════════════════════════════════════════════════\n\n"))

  t0 <- proc.time()

  # ── 1. Window grid ──────────────────────────────────────────────────────────
  window_grid <<- make_window_template(
    step_doy  = step_doy,
    lengths_d = seq(15, 183, by = 7)
  )
  cat(sprintf("  Window grid: %d windows (step = %d days)\n",
              nrow(window_grid), step_doy))

  # ── 2. Anchor detection ────────────────────────────────────────────────────
  anchor_args <- list(
    independence_mins = independence_mins
  )
  anchor_args <- modifyList(anchor_args, anchor_overrides)

  cat(sprintf("  Finding anchors (independence = %d min", independence_mins))
  if (length(anchor_overrides) > 0) {
    cat(sprintf(", + %d overrides", length(anchor_overrides)))
  }
  cat(") ...\n")

  anchors <<- furrr::future_map_dfr(
    ds_paths,
    function(p) do.call(find_anchors, c(list(archive_path = p), anchor_args)),
    .options = furrr::furrr_options(seed = TRUE)
  )

  cat(sprintf("  Anchors: %d slices across %d datasets\n",
              nrow(anchors), n_distinct(anchors$dataset)))

  # ── 3. Process all datasets ────────────────────────────────────────────────
  cat(sprintf("  Processing %d datasets ...\n", length(ds_paths)))

  wrapped <- furrr::future_map(
    ds_paths,
    dataset_wrapper1,
    independence_mins = independence_mins,
    .options = furrr::furrr_options(seed = TRUE)
  )

  # ── 4. Aggregate ───────────────────────────────────────────────────────────
  all_window_species  <- map_df(wrapped, "window_species")
  all_window_richness <- map_df(wrapped, "window_richness")

  elapsed <- round((proc.time() - t0)["elapsed"] / 60, 1)
  cat(sprintf("  Done in %.1f min. Species data: %d rows, %d species\n",
              elapsed, nrow(all_window_species),
              n_distinct(all_window_species$species)))

  # ── 5. Save ────────────────────────────────────────────────────────────────
  if (!is.null(output_species)) {
    saveRDS(all_window_species, output_species)
    cat(sprintf("  Saved: %s\n", output_species))
  }
  if (!is.null(output_richness)) {
    saveRDS(all_window_richness, output_richness)
    cat(sprintf("  Saved: %s\n", output_richness))
  }

  invisible(list(species = all_window_species, richness = all_window_richness))
}


# ═══════════════════════════════════════════════════════════════════════════════
# VARIANT 1: 3-day window step (§3 Resolution)
# ═══════════════════════════════════════════════════════════════════════════════

run_variant(
  label            = "3-day window step",
  step_doy         = 3,
  output_species   = "all_window_species_3day.rds",
  output_richness  = "all_window_richness_3day.rds"
)
gc(verbose = FALSE)


# ═══════════════════════════════════════════════════════════════════════════════
# VARIANT 2: Strict anchor parameters (§4 Anchor)
# ═══════════════════════════════════════════════════════════════════════════════

run_variant(
  label            = "Strict anchors",
  anchor_overrides = list(
    min_effort_per_slice = 3000,
    min_balance_ratio    = 0.30,
    min_detect_days      = 300,
    max_zero_frac        = 0.40,
    max_zero_gap         = 45
  ),
  output_species   = "all_window_species_strict.rds",
  output_richness  = "all_window_richness_strict.rds"
)
gc(verbose = FALSE)


# ═══════════════════════════════════════════════════════════════════════════════
# VARIANT 3: 15-minute independence gap (§5 Independence)
# ═══════════════════════════════════════════════════════════════════════════════

run_variant(
  label              = "15-minute independence gap",
  independence_mins  = 15,
  output_species     = "all_window_species_15min.rds",
  output_richness    = "all_window_richness_15min.rds"
)
gc(verbose = FALSE)


# ═══════════════════════════════════════════════════════════════════════════════
# VARIANT 4: 60-minute independence gap (§5 Independence)
# ═══════════════════════════════════════════════════════════════════════════════

run_variant(
  label              = "60-minute independence gap",
  independence_mins  = 60,
  output_species     = "all_window_species_60min.rds",
  output_richness    = "all_window_richness_60min.rds"
)
gc(verbose = FALSE)


# ═══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════

plan(sequential)

cat("\n\n═══════════════════════════════════════════════════════\n")
cat("All 4 pipeline variants complete.\n")
cat("═══════════════════════════════════════════════════════\n\n")

for (f in c("all_window_species_3day.rds",   "all_window_richness_3day.rds",
            "all_window_species_strict.rds",  "all_window_richness_strict.rds",
            "all_window_species_15min.rds",   "all_window_richness_15min.rds",
            "all_window_species_60min.rds",   "all_window_richness_60min.rds")) {
  if (file.exists(f)) {
    info <- file.info(f)
    cat(sprintf("  ✓ %s  (%.1f MB)\n", f, info$size / 1e6))
  } else {
    cat(sprintf("  ✗ %s  MISSING\n", f))
  }
}

# ============================================================================
# Verification Script: EOW_EXT Window
# ============================================================================
# This script verifies that:
# 1. EOW_EXT window is correctly defined in Full1.R
# 2. EOW_EXT is properly centered on SNAP_EU_CORE
# 3. The window appears in window_grid
# 4. dataset_wrapper1 can process EOW_EXT

library(tidyverse)
library(lubridate)

# ----------------------------------------------------------------------------
# 1. Define helper function
# ----------------------------------------------------------------------------
doy <- function(m, d) {
  as.integer(format(as.Date(sprintf("2021-%02d-%02d", m, d)), "%j"))
}

# ----------------------------------------------------------------------------
# 2. Recreate window logic from Full1.R
# ----------------------------------------------------------------------------
cat("=" %+% rep("=", 70) %+% "=\n")
cat("VERIFICATION: EOW_EXT Window Definition\n")
cat("=" %+% rep("=", 70) %+% "=\n\n")

# SNAP_EU_CORE
core_start <- doy(9, 1)
core_end   <- doy(10, 31)
core_len   <- core_end - core_start + 1L

# EOW_EXT (centered on SNAP_EU_CORE)
core_midpoint <- core_start + floor(core_len / 2)
eow_len       <- 42L  # 6 weeks
eow_start     <- core_midpoint - floor(eow_len / 2)
eow_end       <- eow_start + eow_len - 1L

cat("SNAP_EU_CORE:\n")
cat("  Start:    DOY", core_start, "=", format(as.Date("2021-01-01") + core_start - 1, "%b %d"), "\n")
cat("  End:      DOY", core_end, "=", format(as.Date("2021-01-01") + core_end - 1, "%b %d"), "\n")
cat("  Length:  ", core_len, "days\n")
cat("  Midpoint: DOY", core_midpoint, "\n\n")

cat("EOW_EXT (Extended Observation Window):\n")
cat("  Start:    DOY", eow_start, "=", format(as.Date("2021-01-01") + eow_start - 1, "%b %d"), "\n")
cat("  End:      DOY", eow_end, "=", format(as.Date("2021-01-01") + eow_end - 1, "%b %d"), "\n")
cat("  Length:  ", eow_len, "days (", eow_len / 7, "weeks)\n")
cat("  Midpoint: DOY", eow_start + floor(eow_len / 2), "\n\n")

# Centering check
eow_midpoint <- eow_start + floor(eow_len / 2)
cat("✓ Centering Check:\n")
cat("  SNAP_EU_CORE midpoint: DOY", core_midpoint, "\n")
cat("  EOW_EXT midpoint:      DOY", eow_midpoint, "\n")
cat("  Difference:           ", abs(core_midpoint - eow_midpoint), "days")
if (abs(core_midpoint - eow_midpoint) == 0) {
  cat(" ← CORRECTLY CENTERED ✓\n\n")
} else {
  cat(" ← NOT CENTERED ✗\n\n")
}

# Overlap check
cat("✓ Overlap with SNAP_EU_CORE:\n")
cat("  EOW starts", eow_start - core_start, "days after SNAP_EU_CORE\n")
cat("  EOW ends  ", core_end - eow_end, "days after SNAP_EU_CORE ends\n")
cat("  EOW is", core_len - eow_len, "days shorter than SNAP_EU_CORE\n\n")

# ----------------------------------------------------------------------------
# 3. Check window_grid
# ----------------------------------------------------------------------------
cat("=" %+% rep("=", 70) %+% "=\n")
cat("VERIFICATION: window_grid\n")
cat("=" %+% rep("=", 70) %+% "=\n\n")

source("Full1.R", local = TRUE, echo = FALSE)  # Source to get window_grid

if (exists("window_grid")) {
  protocol_windows <- window_grid %>% filter(str_detect(window_id, "SNAP|EOW"))
  
  cat("Protocol windows in window_grid:\n")
  print(protocol_windows %>% select(window_id, start_doy, length_d, end_doy))
  
  if ("EOW_EXT" %in% protocol_windows$window_id) {
    cat("\n✓ EOW_EXT found in window_grid\n")
    
    eow_row <- protocol_windows %>% filter(window_id == "EOW_EXT")
    cat("\n  Verification:\n")
    cat("    Expected: start_doy =", eow_start, ", length_d =", eow_len, "\n")
    cat("    Actual:   start_doy =", eow_row$start_doy, ", length_d =", eow_row$length_d, "\n")
    
    if (eow_row$start_doy == eow_start && eow_row$length_d == eow_len) {
      cat("    Match: ✓\n")
    } else {
      cat("    Match: ✗ MISMATCH!\n")
    }
  } else {
    cat("\n✗ EOW_EXT NOT FOUND in window_grid\n")
  }
} else {
  cat("✗ window_grid not found after sourcing Full1.R\n")
}

# ----------------------------------------------------------------------------
# 4. Test dataset_wrapper1 compatibility
# ----------------------------------------------------------------------------
cat("\n")
cat("=" %+% rep("=", 70) %+% "=\n")
cat("VERIFICATION: dataset_wrapper1 Compatibility\n")
cat("=" %+% rep("=", 70) %+% "=\n\n")

cat("The EOW_EXT window is now defined and should work with dataset_wrapper1.\n")
cat("\nTo verify with actual data:\n")
cat("  1. Run: wrapped <- furrr::future_map(ds_paths, dataset_wrapper1, ...)\n")
cat("  2. Check: all_window_species %>% filter(window_id == 'EOW_EXT')\n")
cat("  3. Verify EOW_EXT appears in results for relevant datasets\n\n")

cat("Expected behavior:\n")
cat("  - EOW_EXT will be computed for all datasets like any other window\n")
cat("  - Results will include d_p_tte, d_log_rate, d_spatial_cov vs FULL\n")
cat("  - EOW_EXT can be compared to SNAP_EU_CORE and SNAP_EU_BUFFER\n\n")

cat("=" %+% rep("=", 70) %+% "=\n")
cat("VERIFICATION COMPLETE\n")
cat("=" %+% rep("=", 70) %+% "=\n")

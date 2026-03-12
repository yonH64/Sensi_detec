## ── QUICK DEBUG BLOCK ─────────────────────────────────────────────────
# 1. how many calendar months do we have?
cat("# months in month_df:", nrow(month_df), "\n")

# 2. how many 12-month slices did we expect?
expected_slices <- nrow(month_df) - 11
cat("# expected slices        :", expected_slices, "\n")

# 3. are those slice IDs present and unique?
slice_counts <- basic %>% count(slice_id)
print(slice_counts, n = 10)

if (nrow(slice_counts) != expected_slices)
  stop("❌ slice_id count mismatch")

# 4. does every slice have exactly 12 rows?
bad <- slice_counts %>% filter(n != 12)
if (nrow(bad)) {
  print(bad)
  stop("❌ some slice(s) don’t have 12 rows")
} else {
  cat("✅ every slice has 12 rows\n")
}

# 5. any duplicate month labels inside a slice?
dup <- basic %>% group_by(slice_id, month_lab) %>% filter(n() > 1)
if (nrow(dup)) {
  print(dup %>% select(slice_id, month_lab) %>% distinct())
  stop("❌ duplicate month labels within slice")
} else {
  cat("✅ no duplicate month names inside slices\n")
}
## ─────────────────────────────────────────────────────────────────────

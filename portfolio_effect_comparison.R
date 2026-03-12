############################################################
# PORTFOLIO EFFECT COMPARISON
# MSE vs. Absolute Deviation Models
#
# PURPOSE: Compare how species richness (portfolio effect)
# influences MSE vs. absolute deviation to understand
# whether richness primarily reduces:
#   - Bias (central tendency of error)
#   - Variance (precision/spread of estimates)
#   - Both equally
#
# Date: 2026-02-20
############################################################

library(dplyr)
library(tibble)
library(posterior)
library(brms)
library(knitr)

############################################################
# HELPER FUNCTION: Extract portfolio effect from fitted model
############################################################

extract_portfolio_effect <- function(fit, model_name, metric, model_type) {
  
  # Get posterior draws
  draws <- as_draws_df(fit)
  
  # Extract l_n_species coefficient
  if (!"b_l_n_species" %in% names(draws)) {
    warning("Model ", model_name, " does not have l_n_species coefficient")
    return(NULL)
  }
  
  x <- draws$b_l_n_species
  
  tibble(
    metric = metric,
    model_type = model_type,
    model_name = model_name,
    median = median(x),
    mean = mean(x),
    sd = sd(x),
    q2.5 = quantile(x, 0.025),
    q97.5 = quantile(x, 0.975),
    p_negative = mean(x < 0),
    p_positive = mean(x > 0),
    interpretation = case_when(
      mean(x < 0) > 0.95 ~ "Strong negative (reduces error)",
      mean(x > 0) > 0.95 ~ "Strong positive (increases error)",
      TRUE ~ "Uncertain"
    )
  )
}


############################################################
# EXTRACT EFFECTS FROM ALL MODELS
############################################################

# NOTE: This assumes you have fitted the UPDATED models
# If models are not yet fitted, this will produce errors
# Run this AFTER fitting models_abs_UPDATED.R and models_mse_UPDATED.R

portfolio_effects <- bind_rows(
  
  # ════════════════════════════════════════════════════════
  # ABSOLUTE DEVIATION MODELS (|Δ|)
  # ════════════════════════════════════════════════════════
  
  # TTE - Time to Event detection probability
  extract_portfolio_effect(
    fit_abs_tte_base,  # You'll need to name your fitted models this way
    "abs_tte_base", 
    "TTE (p_tte)", 
    "Absolute Deviation"
  ),
  
  # LRR - Log event rate
  extract_portfolio_effect(
    fit_abs_lrr_base,
    "abs_lrr_base",
    "Event Rate (log scale)",
    "Absolute Deviation"
  ),
  
  # SCOV - Spatial coverage
  extract_portfolio_effect(
    fit_abs_scov_base,
    "abs_scov_base",
    "Spatial Coverage",
    "Absolute Deviation"
  ),
  
  # ════════════════════════════════════════════════════════
  # MSE MODELS (Bias² + Variance)
  # ════════════════════════════════════════════════════════
  
  # TTE MSE
  extract_portfolio_effect(
    fit_mse_tte_base,
    "mse_tte_base",
    "TTE (p_tte)",
    "MSE"
  ),
  
  # LRR MSE
  extract_portfolio_effect(
    fit_mse_lrr_base,
    "mse_lrr_base",
    "Event Rate (log scale)",
    "MSE"
  ),
  
  # SCOV MSE
  extract_portfolio_effect(
    fit_mse_scov_base,
    "mse_scov_base",
    "Spatial Coverage",
    "MSE"
  )
)


############################################################
# CREATE COMPARISON TABLE
############################################################

portfolio_comparison <- portfolio_effects %>%
  select(metric, model_type, median, q2.5, q97.5, p_negative) %>%
  mutate(
    # Format credible interval
    ci_95 = sprintf("[%.3f, %.3f]", q2.5, q97.5),
    # Format median with 3 decimals
    median_fmt = sprintf("%.3f", median),
    # Format probability
    p_neg_fmt = sprintf("%.3f", p_negative)
  ) %>%
  select(metric, model_type, median_fmt, ci_95, p_neg_fmt) %>%
  tidyr::pivot_wider(
    names_from = model_type,
    values_from = c(median_fmt, ci_95, p_neg_fmt),
    names_sep = "_"
  )

# Print formatted table
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("PORTFOLIO EFFECT COMPARISON: MSE vs. Absolute Deviation\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
cat("Coefficient: log1p(n_species)\n")
cat("Interpretation: Negative = richness reduces error\n\n")

print(portfolio_comparison, n = Inf)


############################################################
# CALCULATE EFFECT SIZE DIFFERENCES
############################################################

effect_differences <- portfolio_effects %>%
  select(metric, model_type, median, mean) %>%
  tidyr::pivot_wider(
    names_from = model_type,
    values_from = c(median, mean)
  ) %>%
  mutate(
    median_diff = `median_MSE` - `median_Absolute Deviation`,
    mean_diff = `mean_MSE` - `mean_Absolute Deviation`,
    interpretation = case_when(
      median_diff < -0.05 ~ "MSE effect STRONGER (richness reduces variance more)",
      median_diff > 0.05  ~ "Abs effect STRONGER (richness reduces bias more)",
      TRUE ~ "Similar effects (balanced bias-variance reduction)"
    )
  )

cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("EFFECT SIZE DIFFERENCES (MSE - Absolute Deviation)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

print(effect_differences %>% select(metric, median_diff, interpretation), n = Inf)


############################################################
# VISUALIZATION: FOREST PLOT
############################################################

library(ggplot2)

portfolio_plot <- portfolio_effects %>%
  mutate(
    metric = factor(metric, levels = c("TTE (p_tte)", 
                                       "Event Rate (log scale)", 
                                       "Spatial Coverage")),
    model_type = factor(model_type, levels = c("Absolute Deviation", "MSE"))
  ) %>%
  ggplot(aes(x = median, y = metric, color = model_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(xmin = q2.5, xmax = q97.5),
                  position = position_dodge(width = 0.4),
                  size = 0.8, linewidth = 1) +
  scale_color_manual(
    values = c("Absolute Deviation" = "#E69F00", "MSE" = "#56B4E9"),
    name = "Model Type"
  ) +
  labs(
    title = "Portfolio Effect Across Detection Metrics",
    subtitle = "Species richness coefficient: log1p(n_species)",
    x = "Coefficient estimate (95% CrI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_face(face = "bold"),
    plot.subtitle = element_face(color = "gray40")
  )

print(portfolio_plot)

# Save plot
ggsave("portfolio_effect_comparison.png", 
       portfolio_plot, 
       width = 10, height = 6, dpi = 300)


############################################################
# STATISTICAL COMPARISON
############################################################

# Compare posterior distributions directly
# Are MSE and |Δ| effects significantly different?

compare_posteriors <- function(fit_abs, fit_mse, metric_name) {
  
  draws_abs <- as_draws_df(fit_abs)$b_l_n_species
  draws_mse <- as_draws_df(fit_mse)$b_l_n_species
  
  # Difference in effects
  diff <- draws_mse - draws_abs
  
  tibble(
    metric = metric_name,
    diff_median = median(diff),
    diff_q2.5 = quantile(diff, 0.025),
    diff_q97.5 = quantile(diff, 0.975),
    p_mse_stronger = mean(diff < 0),  # MSE more negative
    p_abs_stronger = mean(diff > 0),  # |Δ| more negative
    credibly_different = !(diff_q2.5 < 0 & diff_q97.5 > 0)
  )
}

posterior_comparisons <- bind_rows(
  compare_posteriors(fit_abs_tte_base, fit_mse_tte_base, "TTE (p_tte)"),
  compare_posteriors(fit_abs_lrr_base, fit_mse_lrr_base, "Event Rate (log scale)"),
  compare_posteriors(fit_abs_scov_base, fit_mse_scov_base, "Spatial Coverage")
)

cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("POSTERIOR COMPARISON: Are MSE and |Δ| effects different?\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

print(posterior_comparisons, n = Inf)


############################################################
# INTERPRETATION GUIDE
############################################################

cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("INTERPRETATION GUIDE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Portfolio Effect (l_n_species coefficient):\n")
cat("  • Negative: Higher richness → Lower error\n")
cat("  • Mechanism: Asynchronous activity patterns create temporal averaging\n\n")

cat("Comparing MSE vs. |Δ|:\n")
cat("  • Both negative, similar size:\n")
cat("    → Richness reduces bias AND variance proportionally\n\n")
cat("  • MSE more negative:\n")
cat("    → Richness primarily reduces VARIANCE (precision gains)\n")
cat("    → Estimates more stable but not necessarily more accurate\n\n")
cat("  • |Δ| more negative:\n")
cat("    → Richness primarily reduces BIAS (accuracy gains)\n")
cat("    → Estimates more accurate but not necessarily more precise\n\n")

cat("Expected pattern:\n")
cat("  • Similar effects across MSE and |Δ|\n")
cat("  • Reason: Portfolio effect smooths temporal variation,\n")
cat("    affecting both central tendency (bias) and spread (variance)\n\n")


############################################################
# EXPORT TABLE FOR MANUSCRIPT
############################################################

# Create formatted table for publication
manuscript_table <- portfolio_effects %>%
  mutate(
    estimate = sprintf("%.3f [%.3f, %.3f]", median, q2.5, q97.5),
    p_value = ifelse(p_negative > 0.95, "< 0.05", 
                     ifelse(p_negative > 0.90, "< 0.10", "> 0.10"))
  ) %>%
  select(metric, model_type, estimate, p_value) %>%
  tidyr::pivot_wider(
    names_from = model_type,
    values_from = c(estimate, p_value)
  )

# Export as CSV
write.csv(manuscript_table, 
          "portfolio_effect_manuscript_table.csv", 
          row.names = FALSE)

cat("\n\nTable exported to: portfolio_effect_manuscript_table.csv\n")
cat("Figure saved to: portfolio_effect_comparison.png\n\n")

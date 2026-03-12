############################################################
# VISUALIZING THE LATITUDE SUPPRESSION EFFECT
#
# Creates publication-quality figures demonstrating the
# statistical suppression effect of latitude on protocol
# approximation quality.
############################################################

library(tidyverse)
library(patchwork)
library(ggpubr)
library(scales)

# ══════════════════════════════════════════════════════════
# FIGURE 1: SUPPRESSION EFFECT DEMONSTRATION
# ══════════════════════════════════════════════════════════

# Panel A: Bivariate relationship (latitude alone)
plot_bivariate <- dataset_stats %>%
  left_join(dataset_latitude, by = "dataset") %>%
  ggplot(aes(x = latitude, y = mean_spatial_cov)) +
  geom_point(size = 3, alpha = 0.7, color = "#2C3E50") +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", fill = "#E74C3C", alpha = 0.2) +
  annotate("text", x = -0.5, y = 0.25, 
           label = "β = 0.007, p = 0.62\nR² = 0.015", 
           hjust = 0, size = 3.5, color = "#E74C3C") +
  labs(
    title = "A. Latitude alone (bivariate)",
    x = "Latitude (scaled)",
    y = "Spatial coverage deviation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

# Panel B: Species richness alone
plot_species <- dataset_stats %>%
  ggplot(aes(x = n_species, y = mean_spatial_cov)) +
  geom_point(size = 3, alpha = 0.7, color = "#2C3E50") +
  geom_smooth(method = "lm", se = TRUE, color = "#27AE60", fill = "#27AE60", alpha = 0.2) +
  annotate("text", x = 2.5, y = 0.25, 
           label = "β = -0.030, p < 0.001\nR² = 0.556", 
           hjust = 0, size = 3.5, color = "#27AE60") +
  labs(
    title = "B. Species richness alone (bivariate)",
    x = "Number of species",
    y = "Spatial coverage deviation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

# Panel C: Partial residual plot (latitude controlling for species)
# Residualize both variables
lm_dev_species <- lm(mean_spatial_cov ~ n_species, data = test_data)
lm_lat_species <- lm(latitude ~ n_species, data = test_data)

partial_data <- test_data %>%
  mutate(
    dev_resid = residuals(lm_dev_species),
    lat_resid = residuals(lm_lat_species)
  )

plot_partial <- partial_data %>%
  ggplot(aes(x = lat_resid, y = dev_resid)) +
  geom_point(size = 3, alpha = 0.7, color = "#2C3E50") +
  geom_smooth(method = "lm", se = TRUE, color = "#9B59B6", fill = "#9B59B6", alpha = 0.2) +
  annotate("text", x = -1.5, y = 0.08, 
           label = "β = -0.023, p = 0.04\nR² = 0.661 (full model)", 
           hjust = 0, size = 3.5, color = "#9B59B6") +
  labs(
    title = "C. Latitude effect (controlling for species richness)",
    x = "Latitude residual (after removing species effect)",
    y = "Deviation residual (after removing species effect)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

# Panel D: Path diagram (conceptual)
# Create with ggplot arrows and text
plot_path <- ggplot() +
  # Boxes
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 2.5, ymax = 3.5, 
           fill = NA, color = "black", size = 0.8) +
  annotate("text", x = 1, y = 3, label = "Latitude", size = 4, fontface = "bold") +
  
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 0.5, ymax = 1.5, 
           fill = NA, color = "black", size = 0.8) +
  annotate("text", x = 1, y = 1, label = "Species\nRichness", size = 4, fontface = "bold") +
  
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = 1.5, ymax = 2.5, 
           fill = NA, color = "black", size = 0.8) +
  annotate("text", x = 4, y = 2, label = "Protocol\nDeviation", size = 4, fontface = "bold") +
  
  # Arrows
  # Latitude → Species (negative)
  geom_segment(aes(x = 1, y = 2.5, xend = 1, yend = 1.5), 
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               size = 1, color = "#E74C3C") +
  annotate("text", x = 1.3, y = 2, label = "r = -0.48", size = 3, color = "#E74C3C") +
  
  # Species → Deviation (negative, strong)
  geom_segment(aes(x = 1.5, y = 1, xend = 3.5, yend = 2), 
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               size = 1.5, color = "#27AE60") +
  annotate("text", x = 2.5, y = 1.3, label = "β = -0.038\n(portfolio effect)", 
           size = 3, color = "#27AE60", fontface = "bold") +
  
  # Latitude → Deviation direct (negative, suppressed)
  geom_segment(aes(x = 1.5, y = 3.2, xend = 3.5, yend = 2.3), 
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               size = 1, color = "#9B59B6", linetype = "dashed") +
  annotate("text", x = 2.5, y = 3, label = "β = -0.023\n(suppressed)", 
           size = 3, color = "#9B59B6") +
  
  # Annotations
  annotate("text", x = 2.5, y = 0.2, 
           label = "Indirect path: +0.007 (via richness gradient)\nDirect path: -0.023 (phenological compression)\nNet bivariate effect: ≈ 0 (paths cancel)", 
           size = 3, hjust = 0.5, color = "#34495E") +
  
  labs(title = "D. Path diagram showing suppression mechanism") +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 4)) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

# Combine panels
fig1_suppression <- (plot_bivariate | plot_species) / (plot_partial | plot_path) +
  plot_annotation(
    title = "Figure 1. Latitude Suppression Effect on Protocol Approximation Quality",
    subtitle = "Latitude's effect is masked in bivariate analysis due to confounding with species richness gradient",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "#7F8C8D")
    )
  )

# Save
ggsave("figures/Fig_latitude_suppression.pdf", fig1_suppression, 
       width = 12, height = 10, dpi = 300)
ggsave("figures/Fig_latitude_suppression.png", fig1_suppression, 
       width = 12, height = 10, dpi = 300)


# ══════════════════════════════════════════════════════════
# FIGURE 2: GEOGRAPHIC VISUALIZATION
# ══════════════════════════════════════════════════════════

# Join geographic data
map_data <- dataset_stats %>%
  left_join(dataset_latitude, by = "dataset") %>%
  left_join(
    all_window_species %>%
      group_by(dataset) %>%
      summarise(
        lat_deg = first(latitude),
        lon_deg = first(longitude),
        .groups = "drop"
      ),
    by = "dataset"
  )

# Create map
fig2_map <- ggplot() +
  # Europe basemap
  geom_sf(data = world_eu, fill = "#ECF0F1", color = "#BDC3C7", size = 0.2) +
  # Study sites
  geom_point(data = map_data, 
             aes(x = lon_deg, y = lat_deg, 
                 size = n_species, 
                 fill = mean_spatial_cov),
             shape = 21, color = "black", stroke = 0.8, alpha = 0.8) +
  scale_fill_gradient2(
    low = "#27AE60", mid = "#F39C12", high = "#E74C3C",
    midpoint = median(map_data$mean_spatial_cov),
    name = "Spatial coverage\ndeviation",
    guide = guide_colorbar(barwidth = 1, barheight = 8)
  ) +
  scale_size_continuous(
    range = c(3, 10),
    name = "Species\nrichness",
    breaks = c(2, 4, 6, 8)
  ) +
  coord_sf(xlim = c(-5, 25), ylim = c(42, 63)) +
  labs(
    title = "Figure 2. Geographic Distribution of Protocol Performance",
    subtitle = "Northern sites show lower deviations despite lower species richness",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "#7F8C8D"),
    panel.grid = element_line(color = "#ECF0F1"),
    legend.position = "right"
  )

ggsave("figures/Fig_geographic_pattern.pdf", fig2_map, 
       width = 10, height = 8, dpi = 300)


# ══════════════════════════════════════════════════════════
# FIGURE 3: INTERACTION SURFACE (3D CONTOUR)
# ══════════════════════════════════════════════════════════

# Create prediction grid
pred_grid <- expand.grid(
  latitude_scaled = seq(-1.5, 1.5, length.out = 50),
  n_species = seq(2, 8, length.out = 50)
)

# Fit model to get predictions
lm_full_pred <- lm(mean_spatial_cov ~ latitude * n_species, data = test_data)
pred_grid$deviation <- predict(lm_full_pred, newdata = pred_grid)

# Create contour plot
fig3_contour <- ggplot(pred_grid, aes(x = n_species, y = latitude_scaled, z = deviation)) +
  geom_contour_filled(breaks = seq(0, 0.3, by = 0.03), alpha = 0.8) +
  geom_contour(color = "white", size = 0.3, breaks = seq(0, 0.3, by = 0.03)) +
  geom_point(data = test_data, 
             aes(x = n_species, y = latitude, z = NULL), 
             size = 3, shape = 21, fill = "white", color = "black", stroke = 1) +
  scale_fill_viridis_d(
    option = "plasma",
    name = "Deviation",
    guide = guide_legend(reverse = TRUE)
  ) +
  scale_y_continuous(
    breaks = c(-1.5, -0.75, 0, 0.75, 1.5),
    labels = c("43°N\n(Montpellier)", "48°N", "52°N", "57°N", "62°N\n(Norway)")
  ) +
  labs(
    title = "Figure 3. Joint Effects of Latitude and Species Richness",
    subtitle = "Contour plot showing predicted protocol deviation as function of both predictors",
    x = "Species richness (number of species)",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "#7F8C8D"),
    legend.position = "right"
  )

ggsave("figures/Fig_interaction_surface.pdf", fig3_contour, 
       width = 10, height = 6, dpi = 300)


# ══════════════════════════════════════════════════════════
# SUPPLEMENTARY: COEFFICIENT CHANGE VISUALIZATION
# ══════════════════════════════════════════════════════════

# Compare coefficients across models
coef_comparison <- tibble(
  model = c("Latitude only", "Species only", "Both (latitude)", "Both (species)"),
  estimate = c(
    coef(lm_latitude_only)[2],
    coef(lm_species_only)[2],
    coef(lm_full)[2],
    coef(lm_full)[3]
  ),
  se = c(
    summary(lm_latitude_only)$coef[2,2],
    summary(lm_species_only)$coef[2,2],
    summary(lm_full)$coef[2,2],
    summary(lm_full)$coef[3,2]
  ),
  p_value = c(
    summary(lm_latitude_only)$coef[2,4],
    summary(lm_species_only)$coef[2,4],
    summary(lm_full)$coef[2,4],
    summary(lm_full)$coef[3,4]
  )
) %>%
  mutate(
    ci_lower = estimate - 1.96 * se,
    ci_upper = estimate + 1.96 * se,
    significant = p_value < 0.05
  )

figS1_coef_change <- coef_comparison %>%
  ggplot(aes(x = model, y = estimate, color = significant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#95A5A6") +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper), 
                  size = 1, linewidth = 1.5) +
  geom_text(aes(label = sprintf("p = %.3f", p_value)), 
            nudge_x = 0.3, size = 3) +
  scale_color_manual(values = c("TRUE" = "#27AE60", "FALSE" = "#E74C3C"),
                     labels = c("Non-significant", "Significant (p<0.05)")) +
  coord_flip() +
  labs(
    title = "Supplementary Figure S1. Coefficient Stability Across Model Specifications",
    x = NULL,
    y = "Standardized coefficient estimate (β) ± 95% CI",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

ggsave("figures/FigS1_coefficient_change.pdf", figS1_coef_change, 
       width = 8, height = 5, dpi = 300)


# ══════════════════════════════════════════════════════════
# SUMMARY TABLE FOR MANUSCRIPT
# ══════════════════════════════════════════════════════════

table_suppression <- tibble(
  Model = c("Latitude only", "Species richness only", "Both predictors"),
  Formula = c("deviation ~ latitude", 
              "deviation ~ n_species",
              "deviation ~ latitude + n_species"),
  R_squared = c(
    summary(lm_latitude_only)$r.squared,
    summary(lm_species_only)$r.squared,
    summary(lm_full)$r.squared
  ),
  Beta_latitude = c(
    coef(lm_latitude_only)[2],
    NA,
    coef(lm_full)[2]
  ),
  P_latitude = c(
    summary(lm_latitude_only)$coef[2,4],
    NA,
    summary(lm_full)$coef[2,4]
  ),
  Beta_species = c(
    NA,
    coef(lm_species_only)[2],
    coef(lm_full)[3]
  ),
  P_species = c(
    NA,
    summary(lm_species_only)$coef[2,4],
    summary(lm_full)$coef[3,4]
  )
) %>%
  mutate(
    R_squared = sprintf("%.3f", R_squared),
    Beta_latitude = sprintf("%.4f", Beta_latitude),
    P_latitude = ifelse(is.na(P_latitude), "—", 
                        ifelse(P_latitude < 0.001, "<0.001", sprintf("%.3f", P_latitude))),
    Beta_species = sprintf("%.4f", Beta_species),
    P_species = ifelse(is.na(P_species), "—", 
                       ifelse(P_species < 0.001, "<0.001", sprintf("%.3f", P_species)))
  )

write_csv(table_suppression, "tables/Table_suppression_effect.csv")

cat("═══ VISUALIZATION COMPLETE ═══\n\n")
cat("Generated figures:\n")
cat("  • figures/Fig_latitude_suppression.pdf/png\n")
cat("  • figures/Fig_geographic_pattern.pdf\n")
cat("  • figures/Fig_interaction_surface.pdf\n")
cat("  • figures/FigS1_coefficient_change.pdf\n\n")
cat("Generated tables:\n")
cat("  • tables/Table_suppression_effect.csv\n\n")
cat("All figures ready for manuscript submission.\n")

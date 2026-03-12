# sensitivity_figures.R
# ============================================================================
# Publication figures for the sensitivity surface analysis.
#
# Input:  sensitivity_models_env.RData  (from models_sensitivity_surface.R)
#         Q*.csv files from sensitivity_results.R
#
# Output: Fig1_sensitivity_surface.pdf
#         Fig2_guild_surfaces.pdf
#         Fig3_species_surfaces.pdf
#         Fig4_duration_curves.pdf
#         Fig5_signed_rate_surface.pdf
#         Fig6_richness_surfaces.pdf
#         Fig7_snapshot_evaluation.pdf
#         Fig8_model_comparison.pdf
# ============================================================================

library(mgcv)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

theme_set(theme_minimal(base_size = 11))

cat("=== Sensitivity Surface Figures ===\n\n")

# ── Load ─────────────────────────────────────────────────────────────────────

load("sensitivity_models_env.RData")

q1 <- read.csv("Q1_duration_effect.csv")
q2 <- read.csv("Q2_seasonal_profiles.csv")
q3 <- read.csv("Q3_surface_predictions.csv")
q4 <- read.csv("Q4_species_guild_surfaces.csv")
q6 <- read.csv("Q6_protocol_evaluation.csv")
q8 <- read.csv("Q8_richness_surface.csv")

guild_info <- sens_species |>
  distinct(species_f, guild_major, guild_minor_habitat, guild_minor_diet)

# Month labels for x-axis
month_breaks <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
month_labels <- month.abb

# Protocol coordinates (from helpers.R::protocol_windows())
protocol_df <- tibble(
  protocol   = c("CORE", "BUFFER", "EOW_EARLY", "EOW_LATE"),
  day_start  = c(244, 230, 214, 274),
  window_len = c(61, 89, 60, 60)
)


# ============================================================================
# Fig 1: Main sensitivity surface (all species pooled)
# ============================================================================

cat("  Fig 1: Main surface...\n")

fig1_data <- q3 |>
  filter(metric == "abs_d_lambda") |>
  summarise(mean_pred = mean(mean_pred), .by = c(day_start, window_len))

fig1 <- ggplot(fig1_data, aes(day_start, window_len, fill = mean_pred)) +
  geom_tile(width = 7, height = 7) +
  scale_fill_viridis_c(option = "inferno", direction = -1,
                       name = expression("Predicted |" * Delta * lambda * "|")) +
  scale_x_continuous(breaks = month_breaks, labels = month_labels) +
  # Protocol positions
  geom_point(data = protocol_df, aes(fill = NULL),
             color = "cyan", size = 3, shape = 18) +
  annotate("text", x = 252, y = 55, label = "CORE", color = "cyan",
           size = 2.8, fontface = "bold", hjust = 0) +
  annotate("text", x = 238, y = 83, label = "BUFFER", color = "cyan",
           size = 2.8, fontface = "bold", hjust = 0) +
  annotate("text", x = 200, y = 54, label = "EOW\nEARLY", color = "cyan",
           size = 2.5, fontface = "bold", hjust = 1) +
  annotate("text", x = 282, y = 54, label = "EOW\nLATE", color = "cyan",
           size = 2.5, fontface = "bold", hjust = 0) +
  labs(x = "Window start (day of year)",
       y = "Window duration (days)") +
  theme(legend.position = "right")

ggsave("Fig1_sensitivity_surface.pdf", fig1, width = 8, height = 5)


# ============================================================================
# Fig 2: Guild-specific surfaces
# ============================================================================

cat("  Fig 2: Guild surfaces...\n")

fig2_data <- q3 |>
  filter(metric == "abs_d_lambda", guild_major != "Insectivore")

fig2 <- ggplot(fig2_data, aes(day_start, window_len, fill = mean_pred)) +
  geom_tile(width = 7, height = 7) +
  scale_fill_viridis_c(option = "inferno", direction = -1,
                       name = expression("|" * Delta * lambda * "|")) +
  scale_x_continuous(breaks = c(1, 91, 182, 274), labels = c("J", "A", "J", "O")) +
  facet_wrap(~guild_major, nrow = 1) +
  labs(x = "Window start", y = "Duration (days)") +
  theme(legend.position = "bottom")

ggsave("Fig2_guild_surfaces.pdf", fig2, width = 10, height = 4)


# ============================================================================
# Fig 3: Species-specific surfaces (focal species)
# ============================================================================

cat("  Fig 3: Species surfaces...\n")

focal_species <- c("Vulpes vulpes", "Cervus elaphus", "Capreolus capreolus",
                    "Sus scrofa", "Meles meles", "Lepus europaeus")

# Re-predict on the fine grid for focal species only
mod_lambda <- all_models$detection$abs_d_lambda
sp_levels <- levels(sens_species$species_f)

focal_grid <- expand.grid(
  day_start  = seq(1, 358, by = 7),
  window_len = c(15, 22, 29, 36, 43, 50, 57, 64, 71, 78, 85, 92, 99, 106, 113, 120),
  species_f  = factor(focal_species, levels = sp_levels),
  stringsAsFactors = FALSE
) |>
  as_tibble() |>
  mutate(
    s_bio4 = 0,
    l_trapdays = median(sens_species$l_trapdays),
    l_nsites = median(sens_species$l_nsites),
    s_latitude = 0,
    s_trap_array = 0,
    ds_sp_f = sens_species$ds_sp_f[1]
  )

focal_grid$pred <- predict(mod_lambda, newdata = focal_grid, type = "response",
                           exclude = "s(ds_sp_f)")

fig3 <- ggplot(focal_grid, aes(day_start, window_len, fill = pred)) +
  geom_tile(width = 7, height = 7) +
  scale_fill_viridis_c(option = "inferno", direction = -1,
                       name = expression("|" * Delta * lambda * "|")) +
  scale_x_continuous(breaks = c(1, 91, 182, 274), labels = c("J", "A", "J", "O")) +
  facet_wrap(~species_f, nrow = 2) +
  labs(x = "Window start", y = "Duration (days)") +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "italic"))

ggsave("Fig3_species_surfaces.pdf", fig3, width = 10, height = 6)


# ============================================================================
# Fig 4: Duration-deviation curves by season
# ============================================================================

cat("  Fig 4: Duration curves...\n")

fig4_data <- q2 |>
  filter(guild_major == "All") |>
  mutate(
    season = case_when(
      day_start >= 1   & day_start <= 78  ~ "Winter (Jan-Mar)",
      day_start >= 85  & day_start <= 162 ~ "Spring (Apr-Jun)",
      day_start >= 169 & day_start <= 246 ~ "Summer (Jul-Sep)",
      day_start >= 253 & day_start <= 337 ~ "Autumn (Sep-Dec)",
      TRUE ~ "Winter (Jan-Mar)"
    )
  ) |>
  summarise(mean_pred = mean(mean_pred), .by = c(window_len, season))

fig4 <- ggplot(fig4_data, aes(window_len, mean_pred, color = season)) +
  geom_line(linewidth = 1) +
  # Mark SE protocol durations
  geom_vline(xintercept = 61, linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = 89, linetype = "dotted", alpha = 0.4) +
  annotate("text", x = 63, y = max(fig4_data$mean_pred) * 0.95,
           label = "CORE\n(61d)", size = 2.5, hjust = 0) +
  annotate("text", x = 91, y = max(fig4_data$mean_pred) * 0.85,
           label = "BUFFER\n(89d)", size = 2.5, hjust = 0) +
  scale_color_brewer(type = "qual", palette = "Set1") +
  labs(x = "Window duration (days)",
       y = expression("Mean predicted |" * Delta * lambda * "|"),
       color = "Window\ncentered in") +
  theme(legend.position = c(0.8, 0.75))

ggsave("Fig4_duration_curves.pdf", fig4, width = 6, height = 4.5)


# ============================================================================
# Fig 5: Signed encounter rate surface
# ============================================================================

cat("  Fig 5: Signed rate surface...\n")

fig5_data <- q3 |>
  filter(metric == "d_rate_signed") |>
  summarise(mean_pred = mean(mean_pred), .by = c(day_start, window_len))

fig5 <- ggplot(fig5_data, aes(day_start, window_len, fill = mean_pred)) +
  geom_tile(width = 7, height = 7) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0,
                       name = expression(Delta * "rate (signed)")) +
  scale_x_continuous(breaks = month_breaks, labels = month_labels) +
  labs(x = "Window start (day of year)",
       y = "Window duration (days)") +
  theme(legend.position = "right")

ggsave("Fig5_signed_rate_surface.pdf", fig5, width = 8, height = 5)


# ============================================================================
# Fig 6: Richness & similarity surfaces
# ============================================================================

cat("  Fig 6: Richness surfaces...\n")

q8_wide <- q8 |>
  mutate(
    metric_label = case_when(
      metric == "d_sr_raref"   ~ "Rarefied richness deviation",
      metric == "prop_sr_full" ~ "Proportion of species recovered",
      metric == "rho_lambda"   ~ "Rank correlation (Spearman rho)"
    )
  )

fig6_prop <- q8_wide |>
  filter(metric == "prop_sr_full") |>
  ggplot(aes(day_start, window_len, fill = pred)) +
  geom_tile(width = 7, height = 7) +
  scale_fill_viridis_c(option = "viridis", name = "Proportion") +
  scale_x_continuous(breaks = c(1, 91, 182, 274), labels = c("J", "A", "J", "O")) +
  labs(x = NULL, y = "Duration (days)", title = "Species recovered")

fig6_rho <- q8_wide |>
  filter(metric == "rho_lambda") |>
  ggplot(aes(day_start, window_len, fill = pred)) +
  geom_tile(width = 7, height = 7) +
  scale_fill_viridis_c(option = "viridis", name = expression(rho)) +
  scale_x_continuous(breaks = c(1, 91, 182, 274), labels = c("J", "A", "J", "O")) +
  labs(x = NULL, y = "Duration (days)", title = "Rank preservation")

fig6_sr <- q8_wide |>
  filter(metric == "d_sr_raref") |>
  ggplot(aes(day_start, window_len, fill = pred)) +
  geom_tile(width = 7, height = 7) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                       midpoint = 0, name = expression(Delta * "S")) +
  scale_x_continuous(breaks = c(1, 91, 182, 274), labels = c("J", "A", "J", "O")) +
  labs(x = "Window start", y = "Duration (days)", title = "Richness deviation")

fig6 <- fig6_prop + fig6_rho + fig6_sr + plot_layout(nrow = 1)

ggsave("Fig6_richness_surfaces.pdf", fig6, width = 12, height = 4)


# ============================================================================
# Fig 7: Snapshot Europe evaluation — guild × metric comparison
# ============================================================================

cat("  Fig 7: Protocol evaluation...\n")

q6_long <- q6 |>
  pivot_longer(c(CORE, BUFFER, EOW_EARLY, EOW_LATE),
               names_to = "protocol", values_to = "pred") |>
  filter(guild_major != "Insectivore") |>
  mutate(
    protocol = factor(protocol, levels = c("EOW_EARLY", "CORE", "EOW_LATE", "BUFFER")),
    metric_label = case_when(
      metric == "abs_d_lambda"       ~ "TTE detection rate",
      metric == "abs_d_rate"         ~ "Encounter rate",
      metric == "abs_d_matched_rate" ~ "Spatial coverage rate"
    )
  )

fig7 <- ggplot(q6_long, aes(guild_major, pred, fill = protocol)) +
  geom_col(position = "dodge") +
  facet_wrap(~metric_label, scales = "free_y") +
  scale_fill_manual(
    values = c(EOW_EARLY = "#ff7f0e", CORE = "#1f77b4",
               EOW_LATE = "#d62728", BUFFER = "#2ca02c"),
    labels = c(EOW_EARLY = "EOW Early (60d, Aug 2)",
               CORE = "CORE (61d, Sep 1)",
               EOW_LATE = "EOW Late (60d, Oct 1)",
               BUFFER = "BUFFER (89d, Aug 18)")
  ) +
  labs(x = NULL,
       y = "Predicted absolute deviation",
       fill = "Protocol") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top")

ggsave("Fig7_protocol_evaluation.pdf", fig7, width = 10, height = 4.5)


# ============================================================================
# Fig 8: Model comparison (AIC ladder)
# ============================================================================

cat("  Fig 8: Model comparison...\n")

comp_tbl <- read.csv("model_comparison_table.csv") |>
  mutate(model = reorder(model, -AIC))

fig8 <- ggplot(comp_tbl, aes(model, delta_AIC)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0(round(dev_expl * 100, 1), "%")),
            hjust = -0.1, size = 3) +
  coord_flip() +
  labs(x = NULL, y = expression(Delta * "AIC (relative to best)"),
       title = "Model comparison for |d_lambda|") +
  theme(plot.title = element_text(size = 11))

ggsave("Fig8_model_comparison.pdf", fig8, width = 7, height = 4)


cat("\n=== All figures saved as PDF ===\n")

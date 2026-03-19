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
#         Fig9_benchmark_noise_floor.pdf
#         FigS1_robustness_benchmark.pdf  (supplementary)
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
  window_len = seq(15, 183, by = 7),
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



# ============================================================================
# Fig 9: Benchmark noise floor
# ============================================================================
# Panel A: Inter-annual CV of lambda_full by species × site
# Panel B: Overall SNR degradation with window duration
# Panel C: % observations below noise floor
# Panel D: Species-specific SNR curves for focal species

cat("  Fig 9: Benchmark noise floor...\n")

q9 <- read.csv("Q9_noise_floor.csv")
cv_data <- read.csv("interannual_cv_lambda_full.csv")

# --- Panel A: CV by species ---
fig9a <- cv_data |>
  mutate(species = reorder(species, cv)) |>
  ggplot(aes(cv, species, fill = base_dataset)) +
  geom_col() +
  labs(x = "Inter-annual CV (%)", y = NULL, fill = "Site",
       title = "(a) Benchmark variability") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(face = "italic", size = 7))

# --- Panel B: Overall SNR by window duration ---
snr_duration <- q9 |> filter(section == "duration_summary")

fig9b <- ggplot(snr_duration, aes(window_len, median_snr)) +
  geom_ribbon(aes(ymin = q25_snr, ymax = q75_snr), alpha = 0.2) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  annotate("text", x = 115, y = 1.5, label = "SNR = 1", color = "red",
           size = 3, hjust = 1) +
  scale_y_log10() +
  labs(x = "Window duration (days)", y = "SNR (log scale)",
       title = "(b) Overall SNR by duration")

# --- Panel C: % below noise floor by duration ---
fig9c <- ggplot(snr_duration, aes(window_len, pct_below_1)) +
  geom_col(width = 6) +
  labs(x = "Window duration (days)", y = "% obs. with SNR < 1",
       title = "(c) Below noise floor")

# --- Panel D: Species-specific SNR curves ---
# Recover noise_floor from all_window_species + cv_data
multi_year_bases <- unique(cv_data$base_dataset)
noise_floor_local <- all_window_species |>
  filter(window_id != "FULL",
         !window_id %in% c("SNAP_EU_CORE", "SNAP_EU_BUFFER",
                            "EOW_EARLY", "EOW_LATE")) |>
  mutate(base_dataset = gsub("_slice\\d+$", "", dataset)) |>
  filter(base_dataset %in% multi_year_bases) |>
  left_join(cv_data |> select(base_dataset, species, sd_lambda, cv),
            by = c("base_dataset", "species")) |>
  filter(!is.na(sd_lambda)) |>
  mutate(snr = abs(d_lambda) / sd_lambda)

focal_species <- c("Cervus elaphus", "Alces alces", "Sciurus vulgaris",
                    "Vulpes vulpes", "Ursus arctos", "Lepus europaeus")

focal_curves <- noise_floor_local |>
  filter(species %in% focal_species) |>
  summarise(
    median_snr = median(snr),
    q25_snr    = quantile(snr, 0.25),
    q75_snr    = quantile(snr, 0.75),
    .by = c(species, window_len)
  ) |>
  left_join(cv_data |> summarise(cv = mean(cv), .by = species),
            by = "species") |>
  mutate(label = paste0(species, " (CV=", round(cv), "%)"))

fig9d <- ggplot(focal_curves, aes(window_len, median_snr, color = label)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = q25_snr, ymax = q75_snr, fill = label),
              alpha = 0.08, color = NA) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_log10() +
  labs(x = "Window duration (days)", y = "SNR (log scale)",
       color = NULL, fill = NULL,
       title = "(d) Species-specific noise floor convergence") +
  theme(legend.position = "bottom",
        legend.text = element_text(face = "italic", size = 8)) +
  guides(color = guide_legend(ncol = 2), fill = "none")

# --- Assemble ---
fig9 <- fig9a + (fig9b / fig9c / fig9d) +
  plot_layout(widths = c(1.1, 1)) +
  plot_annotation(
    title = "Benchmark noise floor: inter-annual variability of 12-month detection rate",
    subtitle = "4 multi-year sites, 20 species"
  )

ggsave("figures/Fig9_benchmark_noise_floor.pdf", fig9, width = 14, height = 10)



# ============================================================================
# Fig S1: Benchmark robustness check (supplementary)
# ============================================================================
# Shows that the sensitivity surface shape is preserved when using 180-day
# or 270-day benchmarks instead of the default 365-day benchmark.

cat("  Fig S1: Benchmark robustness check...\n")

rob_summary <- read.csv("robustness_benchmark_summary.csv")
rob         <- readRDS("robustness_check_data.rds")

# --- Panel A: Shape correlations by window duration ---
corr_data <- rob_summary |>
  filter(!is.na(shape_r_365_180)) |>
  select(window_len, shape_r_365_180, shape_r_365_270) |>
  pivot_longer(-window_len, names_to = "comparison", values_to = "correlation") |>
  mutate(comparison = ifelse(comparison == "shape_r_365_180",
                              "365d vs 180d", "365d vs 270d"))

figS1a <- ggplot(corr_data, aes(window_len, correlation, color = comparison)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.9, linetype = "dashed", alpha = 0.4) +
  annotate("text", x = 118, y = 0.91, label = "r = 0.9", size = 2.8, alpha = 0.5) +
  scale_color_manual(values = c("365d vs 180d" = "#1b9e77", "365d vs 270d" = "#d95f02")) +
  labs(x = "Sub-window duration (days)",
       y = "Shape correlation (r)",
       color = "Benchmark\ncomparison",
       title = "(a) Seasonal shape preservation") +
  theme(legend.position = c(0.35, 0.25))

# --- Panel B: Mean |d_lambda| by benchmark ---
mean_data <- rob_summary |>
  filter(!is.na(shape_r_365_180)) |>
  select(window_len, mean_abs_d_365, mean_abs_d_270, mean_abs_d_180) |>
  pivot_longer(-window_len, names_to = "benchmark", values_to = "mean_abs_d") |>
  mutate(benchmark = case_when(
    benchmark == "mean_abs_d_365" ~ "365-day",
    benchmark == "mean_abs_d_270" ~ "270-day",
    benchmark == "mean_abs_d_180" ~ "180-day"
  ))

figS1b <- ggplot(mean_data, aes(window_len, mean_abs_d, color = benchmark)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_color_manual(values = c("365-day" = "grey30", "270-day" = "#d95f02",
                                "180-day" = "#1b9e77")) +
  labs(x = "Sub-window duration (days)",
       y = expression("Mean |" * Delta * lambda * "|"),
       color = "Benchmark\nduration",
       title = "(b) Deviation magnitude by benchmark") +
  theme(legend.position = c(0.75, 0.75))

# --- Panel C: Example seasonal profiles at 29d and 85d ---
profile_data <- rob |>
  filter(window_len %in% c(29, 85)) |>
  summarise(
    abs_d_365 = mean(abs_d_lambda_365, na.rm = TRUE),
    abs_d_270 = mean(abs_d_lambda_270, na.rm = TRUE),
    abs_d_180 = mean(abs_d_lambda_180, na.rm = TRUE),
    .by = c(day_start, window_len)
  ) |>
  pivot_longer(starts_with("abs_d_"), names_to = "benchmark", values_to = "abs_d") |>
  mutate(
    benchmark = case_when(
      benchmark == "abs_d_365" ~ "365-day",
      benchmark == "abs_d_270" ~ "270-day",
      benchmark == "abs_d_180" ~ "180-day"
    ),
    duration_label = paste0(window_len, "-day window")
  )

figS1c <- ggplot(profile_data, aes(day_start, abs_d, color = benchmark)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~duration_label, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("365-day" = "grey30", "270-day" = "#d95f02",
                                "180-day" = "#1b9e77")) +
  scale_x_continuous(breaks = month_breaks, labels = month_labels) +
  labs(x = "Window start (day of year)",
       y = expression("Mean |" * Delta * lambda * "|"),
       color = "Benchmark",
       title = "(c) Seasonal profiles by benchmark") +
  theme(legend.position = "bottom")

# --- Assemble ---
figS1 <- (figS1a | figS1b) / figS1c +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title = "Benchmark robustness check",
    subtitle = "Sensitivity surface shape is preserved across 180-day, 270-day, and 365-day benchmarks"
  )

ggsave("figures/FigS1_robustness_benchmark.pdf", figS1, width = 12, height = 9)


cat("\n=== All figures saved as PDF ===\n")

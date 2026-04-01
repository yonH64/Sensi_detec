# ───────────────────────────────────────────────────────────────────────────────
# Leave-one-dataset-out cross-validation for sensitivity surface models
# ───────────────────────────────────────────────────────────────────────────────
#
# Purpose: Evaluate out-of-sample predictive performance of the M6 model
#   (species-specific tensor product surfaces) by iteratively holding out each
#   base dataset, refitting on the remaining data, and predicting on the
#   held-out observations with the random effect excluded.
#
# Inputs:
#   - sensitivity_models_env.RData  (fitted models + prepared data)
#
# Outputs:
#   - loo_cv_predictions.csv             (obs-level predictions, lambda)
#   - loo_cv_predictions_matched_rate.csv
#   - loo_cv_predictions_rate.csv
#   - loo_cv_overall_summary.csv
#   - loo_cv_by_dataset.csv              (per-dataset lambda summary)
#   - loo_cv_by_dataset_all_metrics.csv  (per-dataset, all 3 metrics)
#   - loo_cv_by_species.csv              (per-species lambda summary)
#   - loo_cv_by_species_all_metrics.csv  (per-species, all 3 metrics)
#   - figures/FigS_loo_cv_diagnostics.pdf
#   - figures/FigS_loo_cv_cross_metric.pdf
#
# Runtime: ~150 minutes (27 datasets × 3 metrics × ~2 min/fit)
# ───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(mgcv)
library(patchwork)

# ── Load data ─────────────────────────────────────────────────────────────────

load("sensitivity_models_env.RData")

# Derive base_dataset (strip _sliceN suffix)
sens_species <- sens_species |>
  mutate(base_dataset = sub("_slice\\d+$", "", dataset))

base_datasets <- unique(sens_species$base_dataset)
cat("Base datasets:", length(base_datasets), "\n")
cat("Total observations:", nrow(sens_species), "\n")

# ── Metrics to evaluate ──────────────────────────────────────────────────────

metrics <- c("abs_d_lambda", "abs_d_matched_rate", "abs_d_rate")

# ── LOO-CV loop ──────────────────────────────────────────────────────────────

loo_results <- setNames(
  lapply(metrics, function(m) vector("list", length(base_datasets))),
  metrics
)

cat(sprintf(
  "\nStarting LOO-CV: %d datasets × %d metrics = %d fits\n",
  length(base_datasets), length(metrics), length(base_datasets) * length(metrics)
))
cat("Started:", format(Sys.time(), "%H:%M:%S"), "\n\n")

t_total <- proc.time()

for (i in seq_along(base_datasets)) {
  ds <- base_datasets[i]

  train_idx <- sens_species$base_dataset != ds
  test_idx  <- sens_species$base_dataset == ds

  train_data <- sens_species[train_idx, ]
  test_data  <- sens_species[test_idx, ]

  # Drop unused RE factor levels

  train_data$ds_sp_f <- droplevels(train_data$ds_sp_f)

  # Identify species present in test but absent from training
  missing_spp <- setdiff(
    unique(as.character(test_data$species_f)),
    unique(as.character(train_data$species_f))
  )
  if (length(missing_spp) > 0) {
    test_data <- test_data[!test_data$species_f %in% missing_spp, ]
    cat(sprintf("  [%s] Dropped %d spp: %s\n",
                ds, length(missing_spp), paste(missing_spp, collapse = ", ")))
  }

  if (nrow(test_data) == 0) {
    cat(sprintf("  [%2d/%d] %-30s  SKIPPED (no test rows)\n",
                i, length(base_datasets), ds))
    next
  }

  # Prepare test data for prediction (dummy RE level)
  test_data_pred <- test_data
  test_data_pred$ds_sp_f <- train_data$ds_sp_f[1]

  # Fit and predict for each metric
  r_vals <- c()
  for (metric in metrics) {
    t0 <- proc.time()

    fml <- as.formula(paste0(
      metric,
      " ~ te(day_start, window_len, bs = c('cc', 'tp'), k = c(8, 6), by = species_f) + ",
      "species_f + l_trapdays + l_nsites + s_latitude + s_trap_array + s(ds_sp_f, bs = 're')"
    ))

    fit_i <- tryCatch(
      bam(fml, data = train_data, family = Gamma(link = "log"),
          knots = list(day_start = c(0, 365)), discrete = TRUE, nthreads = 4),
      error = function(e) {
        cat(sprintf("    ERROR [%s]: %s\n", metric, conditionMessage(e)))
        NULL
      }
    )

    elapsed_i <- (proc.time() - t0)["elapsed"]

    if (!is.null(fit_i)) {
      preds <- predict(fit_i, newdata = test_data_pred, type = "response",
                       exclude = "s(ds_sp_f)")

      loo_results[[metric]][[i]] <- tibble(
        base_dataset = ds,
        dataset = test_data$dataset,
        species = as.character(test_data$species_f),
        day_start = test_data$day_start,
        window_len = test_data$window_len,
        observed = test_data[[metric]],
        predicted = as.numeric(preds),
        metric = metric,
        fit_time = elapsed_i
      )

      r_vals <- c(r_vals, round(cor(test_data[[metric]], preds), 3))
    } else {
      r_vals <- c(r_vals, NA)
    }

    rm(fit_i)
    gc(verbose = FALSE)
  }

  cat(sprintf("  [%2d/%d] %-30s  r: %s\n",
              i, length(base_datasets), ds,
              paste(sprintf("%s=%.3f", c("lam", "mrt", "rat"), r_vals), collapse = "  ")))
}

t_elapsed <- (proc.time() - t_total)["elapsed"]
cat(sprintf("\nDone. Total: %.1f minutes\n", t_elapsed / 60))

# ── Combine results ──────────────────────────────────────────────────────────

loo_lambda  <- bind_rows(loo_results[["abs_d_lambda"]])
loo_matched <- bind_rows(loo_results[["abs_d_matched_rate"]])
loo_rate    <- bind_rows(loo_results[["abs_d_rate"]])

# ── Overall summary ──────────────────────────────────────────────────────────

overall <- tibble(
  metric = metrics,
  n_obs = c(nrow(loo_lambda), nrow(loo_matched), nrow(loo_rate)),
  n_species_evaluated = c(n_distinct(loo_lambda$species),
                          n_distinct(loo_matched$species),
                          n_distinct(loo_rate$species)),
  r_response = c(
    cor(loo_lambda$observed, loo_lambda$predicted),
    cor(loo_matched$observed, loo_matched$predicted),
    cor(loo_rate$observed, loo_rate$predicted)
  ),
  r_log = c(
    cor(log(loo_lambda$observed), log(loo_lambda$predicted)),
    cor(log(loo_matched$observed), log(loo_matched$predicted)),
    cor(log(loo_rate$observed), log(loo_rate$predicted))
  ),
  rmse = c(
    sqrt(mean((loo_lambda$observed - loo_lambda$predicted)^2)),
    sqrt(mean((loo_matched$observed - loo_matched$predicted)^2)),
    sqrt(mean((loo_rate$observed - loo_rate$predicted)^2))
  ),
  mae = c(
    mean(abs(loo_lambda$observed - loo_lambda$predicted)),
    mean(abs(loo_matched$observed - loo_matched$predicted)),
    mean(abs(loo_rate$observed - loo_rate$predicted))
  )
) |>
  mutate(r2_response = r_response^2, r2_log = r_log^2)

cat("\n=== Overall LOO-CV summary ===\n")
print(overall)

# ── Per-dataset summary (lambda) ─────────────────────────────────────────────

loo_by_ds <- loo_lambda |>
  group_by(base_dataset) |>
  summarise(
    n_test = n(),
    n_species = n_distinct(species),
    r_response = cor(observed, predicted),
    r_log = cor(log(observed), log(predicted)),
    rmse = sqrt(mean((observed - predicted)^2)),
    mae = mean(abs(observed - predicted)),
    mean_obs = mean(observed),
    mean_pred = mean(predicted),
    bias = mean(predicted - observed),
    .groups = "drop"
  ) |>
  arrange(desc(r_response))

# ── Per-dataset summary (all 3 metrics) ─────────────────────────────────────

ds_all_metrics <- loo_by_ds |>
  select(base_dataset, r_lambda = r_response) |>
  left_join(
    loo_matched |>
      group_by(base_dataset) |>
      summarise(r_matched = cor(observed, predicted), .groups = "drop"),
    by = "base_dataset"
  ) |>
  left_join(
    loo_rate |>
      group_by(base_dataset) |>
      summarise(r_rate = cor(observed, predicted), .groups = "drop"),
    by = "base_dataset"
  ) |>
  mutate(mean_r = (r_lambda + r_matched + r_rate) / 3) |>
  arrange(desc(mean_r))

# ── Per-species summary (all 3 metrics) ─────────────────────────────────────

make_sp_summary <- function(df, metric_name) {
  df |>
    group_by(species) |>
    summarise(
      n_obs = n(),
      n_datasets = n_distinct(base_dataset),
      r = cor(observed, predicted),
      .groups = "drop"
    ) |>
    mutate(metric = metric_name)
}

loo_sp_all <- bind_rows(
  make_sp_summary(loo_lambda, "abs_d_lambda"),
  make_sp_summary(loo_matched, "abs_d_matched_rate"),
  make_sp_summary(loo_rate, "abs_d_rate")
)

# ── Save results ─────────────────────────────────────────────────────────────

write_csv(loo_lambda, "loo_cv_predictions.csv")
write_csv(loo_matched, "loo_cv_predictions_matched_rate.csv")
write_csv(loo_rate, "loo_cv_predictions_rate.csv")
write_csv(overall, "loo_cv_overall_summary.csv")
write_csv(loo_by_ds, "loo_cv_by_dataset.csv")
write_csv(ds_all_metrics, "loo_cv_by_dataset_all_metrics.csv")
write_csv(loo_sp_all, "loo_cv_by_species_all_metrics.csv")

cat("\nCSV files saved.\n")

# ── Figure: Lambda LOO-CV diagnostics (4-panel) ─────────────────────────────

r_overall   <- cor(loo_lambda$observed, loo_lambda$predicted)
r_log_overall <- cor(log(loo_lambda$observed), log(loo_lambda$predicted))

# Panel 1: Observed vs predicted (log scale)
set.seed(4821)
p1 <- ggplot(loo_lambda |> sample_n(min(20000, n())),
             aes(log(predicted), log(observed))) +
  geom_point(alpha = 0.1, size = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.6) +
  labs(x = "log(predicted)", y = "log(observed)",
       title = sprintf("LOO-CV: r = %.2f (log scale)", r_log_overall)) +
  theme_minimal(base_size = 10)

# Panel 2: Per-dataset correlation
ds_order <- loo_by_ds |> arrange(r_response) |> pull(base_dataset)
loo_by_ds$base_dataset_f <- factor(loo_by_ds$base_dataset, levels = ds_order)

p2 <- ggplot(loo_by_ds, aes(r_response, base_dataset_f)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = r_overall, colour = "red", linetype = "dashed") +
  labs(x = "Pearson r (response scale)", y = NULL,
       title = "Per-dataset predictive correlation") +
  theme_minimal(base_size = 9)

# Panel 3: Per-species correlation
sp_lambda <- loo_sp_all |> filter(metric == "abs_d_lambda")
sp_order <- sp_lambda |> arrange(r) |> pull(species)
sp_lambda$species_f <- factor(sp_lambda$species, levels = sp_order)

p3 <- ggplot(sp_lambda, aes(r, species_f)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dotted") +
  labs(x = "Pearson r (response scale)", y = NULL,
       title = "Per-species predictive correlation") +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(face = "italic"))

# Panel 4: By duration bin
loo_lambda <- loo_lambda |>
  mutate(duration_bin = cut(window_len, breaks = c(0, 30, 60, 90, 120, 183),
                            labels = c("15-30", "31-60", "61-90", "91-120", "121-183")))

loo_by_dur <- loo_lambda |>
  group_by(duration_bin) |>
  summarise(r = cor(observed, predicted), .groups = "drop")

p4 <- ggplot(loo_by_dur, aes(duration_bin, r)) +
  geom_col(width = 0.6) +
  labs(x = "Window duration (days)", y = "Pearson r",
       title = "Predictive correlation by duration") +
  theme_minimal(base_size = 10)

fig_diag <- (p1 + p4) / (p2 + p3) +
  plot_annotation(title = "Leave-one-dataset-out cross-validation (|d_lambda|, M6)")

ggsave("figures/FigS_loo_cv_diagnostics.pdf", fig_diag, width = 12, height = 10)
cat("Saved figures/FigS_loo_cv_diagnostics.pdf\n")

# ── Figure: Cross-metric comparison (3-panel) ───────────────────────────────

# Panel A: Overall r barplot
p_overall <- ggplot(overall, aes(reorder(metric, r_response), r_response)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", r_response)), hjust = -0.1, size = 3.5) +
  coord_flip(ylim = c(0, 0.7)) +
  labs(x = NULL, y = "Overall Pearson r (response scale)",
       title = "LOO-CV overall performance") +
  scale_x_discrete(labels = c("abs_d_rate" = "|d_rate|",
                               "abs_d_matched_rate" = "|d_matched_rate|",
                               "abs_d_lambda" = "|d_lambda|")) +
  theme_minimal(base_size = 10)

# Panel B: Per-dataset scatter (lambda vs matched_rate)
p_ds <- ggplot(ds_all_metrics, aes(r_lambda, r_matched)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  ggrepel::geom_text_repel(
    data = ds_all_metrics |> filter(mean_r < 0.45 | abs(r_lambda - r_matched) > 0.25),
    aes(label = base_dataset), size = 2.5, max.overlaps = 15
  ) +
  labs(x = "LOO r (|d_lambda|)", y = "LOO r (|d_matched_rate|)",
       title = "Per-dataset LOO r") +
  coord_equal(xlim = c(-0.1, 1), ylim = c(-0.1, 1)) +
  theme_minimal(base_size = 10)

# Panel C: Per-species heatmap
sp_heat <- loo_sp_all |>
  mutate(metric = recode(metric,
    "abs_d_lambda" = "|d_lambda|",
    "abs_d_matched_rate" = "|d_matched_rate|",
    "abs_d_rate" = "|d_rate|"
  )) |>
  mutate(species = reorder(species, r, FUN = function(x) mean(x, na.rm = TRUE)))

p_sp <- ggplot(sp_heat, aes(metric, species, fill = r)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.5) +
  scale_fill_gradient2(low = "firebrick", mid = "white", high = "steelblue",
                       midpoint = 0.4, limits = c(-0.3, 0.9)) +
  labs(x = NULL, y = NULL, title = "Per-species LOO r by metric", fill = "r") +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(face = "italic"))

fig_cross <- (p_overall + p_ds) / p_sp + plot_layout(heights = c(1, 1.5))

ggsave("figures/FigS_loo_cv_cross_metric.pdf", fig_cross, width = 11, height = 10)
cat("Saved figures/FigS_loo_cv_cross_metric.pdf\n")

cat("\nLOO-CV analysis complete.\n")

# resolution_sensitivity.R
# ============================================================================
# Resolution sensitivity: does the 7-day step size for window start positions
# affect the fitted sensitivity surface?
#
# Tests three step sizes:
#   3-day  : ~122 start positions × 25 durations (full pipeline re-run)
#   7-day  : ~53 start positions × 25 durations (current baseline)
#   14-day : ~27 start positions × 25 durations (subset of 7-day data)
#
# Output:
#   resolution_sensitivity_summary.csv
#   figures/resolution_sensitivity.pdf
# ============================================================================

library(tidyverse)
library(mgcv)
library(furrr)
library(patchwork)

N_THREADS <- 4
cc_knots  <- list(day_start = c(0, 365))

cat("=== Resolution Sensitivity Analysis ===\n\n")

# ── Shared utilities (same as threshold_sensitivity_joint.R) ─────────────────

species_meta  <- read.csv("species_meta.csv")
dataset_meta  <- read.csv("dataset_metadata.csv")

be_leuven_env <- tibble(
  bio4_temp_seasonality   = 559.2,
  bio12_annual_precip     = 774,
  bio15_precip_seasonality = 11.2,
  ndvi_amplitude          = 0.449
)

join_env_covariates <- function(df, meta, be_fix) {
  env_cols  <- c("bio4_temp_seasonality", "bio12_annual_precip",
                 "bio15_precip_seasonality", "ndvi_amplitude")
  join_cols <- c(env_cols, "centroid_lat")
  df <- df |> mutate(dataset_base = str_remove(dataset, "_slice\\d+$"))
  df <- df |>
    left_join(meta |> select(dataset_name, all_of(join_cols)) |>
                rename(dataset = dataset_name), by = "dataset")
  lookup <- meta |> select(dataset_name, all_of(join_cols)) |>
    rename(dataset_base = dataset_name)
  for (col in join_cols) {
    missing_idx <- which(is.na(df[[col]]))
    if (length(missing_idx) > 0) {
      fill_vals <- lookup[[col]][match(df$dataset_base[missing_idx], lookup$dataset_base)]
      df[[col]][missing_idx] <- fill_vals
    }
  }
  for (col in env_cols) {
    be_rows <- grepl("^BE-Leuven", df$dataset) & is.na(df[[col]])
    if (any(be_rows)) df[[col]][be_rows] <- be_fix[[col]]
  }
  df
}

prep_species_data <- function(ws_data, min_rows_per_species = 50) {
  out <- ws_data |>
    filter(!window_id %in% c("SNAP_EU_CORE", "SNAP_EU_BUFFER",
                              "EOW_EARLY", "EOW_LATE")) |>
    mutate(
      day_start          = as.integer(str_extract(window_id, "(?<=d)\\d+")),
      abs_d_lambda       = pmax(abs(d_lambda), 1e-10),
      abs_d_rate         = pmax(abs(d_rate), 1e-10),
      abs_d_matched_rate = abs(d_matched_rate),
      l_trapdays         = log1p(trap_days_window),
      l_nsites           = log1p(n_sites),
      ds_sp              = interaction(dataset, species, drop = TRUE)
    ) |>
    left_join(
      species_meta |> select(species, guild_major, guild_minor_habitat, guild_minor_diet),
      by = "species"
    ) |>
    join_env_covariates(dataset_meta, be_leuven_env)

  spp_counts <- out |> count(species)
  keep_spp   <- spp_counts |> filter(n >= min_rows_per_species) |> pull(species)

  out |>
    filter(!is.na(guild_major), species %in% keep_spp) |>
    mutate(
      s_bio4       = as.numeric(scale(bio4_temp_seasonality)),
      s_latitude   = as.numeric(scale(latitude)),
      s_trap_array = as.numeric(scale(log1p(trap_array))),
      guild_major  = factor(guild_major),
      species_f    = factor(species),
      dataset_f    = factor(dataset),
      ds_sp_f      = factor(ds_sp)
    )
}

fit_m6_lambda <- function(data, label = "") {
  cat(sprintf("  Fitting M6 %-40s ... ", label))
  t0 <- proc.time()
  fit <- bam(
    abs_d_lambda ~
      te(day_start, window_len, bs = c("cc", "tp"), k = c(8, 6),
         by = species_f) +
      species_f +
      s_bio4 + l_trapdays + l_nsites + s_latitude + s_trap_array +
      s(ds_sp_f, bs = "re"),
    data     = data,
    family   = Gamma(link = "log"),
    method   = "fREML",
    discrete = TRUE,
    nthreads = N_THREADS,
    knots    = cc_knots
  )
  elapsed <- round((proc.time() - t0)[3], 1)
  cat(sprintf("done (%ss, dev.expl = %.1f%%)\n", elapsed,
              100 * summary(fit)$dev.expl))
  fit
}


# ── 1. 3-day step: full pipeline re-run ─────────────────────────────────────

cat("── Step 1: 3-day resolution (pipeline re-run) ────────────────────\n\n")

window_grid_3d <- make_window_template(step_doy = 3, lengths_d = seq(15, 183, by = 7))
cat("  3-day grid:", nrow(window_grid_3d), "rows\n")

plan(multisession, workers = 4)

# Temporarily override window_grid for the pipeline
window_grid_orig <- window_grid
window_grid <<- window_grid_3d

t0 <- proc.time()
wrapped_3d <- furrr::future_map(
  ds_paths,
  dataset_wrapper1,
  .options = furrr::furrr_options(seed = TRUE)
)
t_pipe <- round((proc.time() - t0)[3] / 60, 1)
cat(sprintf("  Pipeline done in %.1f minutes.\n", t_pipe))

# Restore original
window_grid <<- window_grid_orig

aws_3d <- purrr::map_df(wrapped_3d, "window_species")
sens_3d <- aws_3d |> prep_species_data()
cat(sprintf("  3-day: %d rows, %d species\n", nrow(sens_3d), n_distinct(sens_3d$species)))


# ── 2. 7-day: load baseline ─────────────────────────────────────────────────

cat("\n── Step 2: 7-day baseline (from disk) ──────────────────────────\n\n")

sens_7d <- readRDS("sensitivity_species_data.rds")
cat(sprintf("  7-day: %d rows, %d species\n", nrow(sens_7d), n_distinct(sens_7d$species)))


# ── 3. 14-day: subset from 7-day data ───────────────────────────────────────

cat("\n── Step 3: 14-day (subset of 7-day) ────────────────────────────\n\n")

# 7-day starts are at 1, 8, 15, 22, 29, ...
# 14-day starts: keep every other one → 1, 15, 29, 43, ...
keep_starts_14d <- seq(1, 358, by = 14)
sens_14d <- sens_7d |>
  filter(day_start %in% keep_starts_14d)
cat(sprintf("  14-day: %d rows, %d species, %d start positions\n",
            nrow(sens_14d), n_distinct(sens_14d$species),
            n_distinct(sens_14d$day_start)))


# ── 4. Fit M6 for each resolution ───────────────────────────────────────────

cat("\n── Step 4: Model fitting ──────────────────────────────────────\n\n")

mod_3d  <- fit_m6_lambda(sens_3d,  "3-day step")
mod_7d  <- readRDS("sensitivity_gam_models.rds")$detection$abs_d_lambda
cat(sprintf("  7-day: loaded (dev.expl = %.1f%%)\n", 100 * summary(mod_7d)$dev.expl))
mod_14d <- fit_m6_lambda(sens_14d, "14-day step")


# ── 5. Compare on a common prediction grid ──────────────────────────────────

cat("\n── Step 5: Comparison ─────────────────────────────────────────\n\n")

shared_spp <- Reduce(intersect, list(
  unique(sens_3d$species), unique(sens_7d$species), unique(sens_14d$species)
))
cat("  Shared species:", length(shared_spp), "\n")

pred_grid <- expand.grid(
  day_start  = seq(1, 358, by = 7),
  window_len = seq(15, 183, by = 1),
  species_f  = factor(shared_spp, levels = levels(sens_7d$species_f)),
  stringsAsFactors = FALSE
) |>
  as_tibble() |>
  mutate(s_bio4 = 0, l_trapdays = median(sens_7d$l_trapdays),
         l_nsites = median(sens_7d$l_nsites), s_latitude = 0,
         s_trap_array = 0, ds_sp_f = sens_7d$ds_sp_f[1])

pred_grid$pred_3d  <- predict(mod_3d,  newdata = pred_grid, type = "response", exclude = "s(ds_sp_f)")
pred_grid$pred_7d  <- predict(mod_7d,  newdata = pred_grid, type = "response", exclude = "s(ds_sp_f)")
pred_grid$pred_14d <- predict(mod_14d, newdata = pred_grid, type = "response", exclude = "s(ds_sp_f)")

cat("Surface correlations vs 7-day baseline:\n")
cat("  3d vs 7d: ", round(cor(pred_grid$pred_3d,  pred_grid$pred_7d),  4), "\n")
cat("  14d vs 7d:", round(cor(pred_grid$pred_14d, pred_grid$pred_7d), 4), "\n")

# Duration curves
dur_comp <- pred_grid |>
  summarise(
    pred_3d  = mean(pred_3d),
    pred_7d  = mean(pred_7d),
    pred_14d = mean(pred_14d),
    .by = window_len
  ) |>
  pivot_longer(-window_len, names_to = "resolution", values_to = "mean_pred") |>
  mutate(resolution = case_when(
    resolution == "pred_3d"  ~ "3-day step",
    resolution == "pred_7d"  ~ "7-day step",
    resolution == "pred_14d" ~ "14-day step"
  ))

cat("\nDuration curve values:\n")
dur_comp |>
  pivot_wider(names_from = resolution, values_from = mean_pred) |>
  filter(window_len %in% c(15, 60, 120, 183)) |>
  mutate(across(-window_len, ~ round(.x, 5))) |>
  print()

# Seasonal profiles at 60d
seas_comp <- pred_grid |>
  filter(window_len == 60) |>
  summarise(
    pred_3d  = mean(pred_3d),
    pred_7d  = mean(pred_7d),
    pred_14d = mean(pred_14d),
    .by = day_start
  ) |>
  pivot_longer(-day_start, names_to = "resolution", values_to = "mean_pred") |>
  mutate(resolution = case_when(
    resolution == "pred_3d"  ~ "3-day step",
    resolution == "pred_7d"  ~ "7-day step",
    resolution == "pred_14d" ~ "14-day step"
  ))

# Parametric coefficients
extract_p <- function(mod, lab) {
  s <- summary(mod)
  as.data.frame(s$p.table) |>
    mutate(term = rownames(s$p.table), resolution = lab) |>
    as_tibble()
}

params <- bind_rows(extract_p(mod_3d, "3d"), extract_p(mod_7d, "7d"), extract_p(mod_14d, "14d")) |>
  filter(term %in% c("s_bio4", "l_trapdays", "l_nsites", "s_latitude", "s_trap_array")) |>
  select(resolution, term, Estimate, `Std. Error`, `Pr(>|t|)`) |>
  arrange(term, resolution)

cat("\nParametric coefficients:\n")
params |>
  mutate(across(c(Estimate, `Std. Error`), ~round(.x, 4)),
         `Pr(>|t|)` = formatC(`Pr(>|t|)`, format = "g", digits = 3)) |>
  print(n = 20)

# Summary table
res_summary <- tibble(
  resolution = c("3-day", "7-day", "14-day"),
  n_start_positions = c(n_distinct(sens_3d$day_start),
                        n_distinct(sens_7d$day_start),
                        n_distinct(sens_14d$day_start)),
  n_rows     = c(nrow(sens_3d), nrow(sens_7d), nrow(sens_14d)),
  n_species  = c(n_distinct(sens_3d$species), n_distinct(sens_7d$species),
                 n_distinct(sens_14d$species)),
  dev_expl   = round(c(summary(mod_3d)$dev.expl, summary(mod_7d)$dev.expl,
                        summary(mod_14d)$dev.expl), 3),
  surface_cor_vs_7d = c(
    round(cor(pred_grid$pred_3d, pred_grid$pred_7d), 4),
    1.0,
    round(cor(pred_grid$pred_14d, pred_grid$pred_7d), 4)
  )
)
print(res_summary)
write.csv(res_summary, "resolution_sensitivity_summary.csv", row.names = FALSE)


# ── 6. Figures ──────────────────────────────────────────────────────────────

p_dur <- ggplot(dur_comp, aes(window_len, mean_pred, color = resolution)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("3-day step" = "#2166AC", "7-day step" = "#1B1B1B",
                                "14-day step" = "#B2182B")) +
  labs(x = "Window duration (days)", y = "Mean |d_lambda|",
       color = "Start-day\nresolution", title = "Duration curves") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

p_seas <- ggplot(seas_comp, aes(day_start, mean_pred, color = resolution)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("3-day step" = "#2166AC", "7-day step" = "#1B1B1B",
                                "14-day step" = "#B2182B")) +
  scale_x_continuous(breaks = c(1, 91, 182, 274), labels = c("Jan", "Apr", "Jul", "Oct")) +
  labs(x = "Window start", y = "Mean |d_lambda|",
       color = "Start-day\nresolution", title = "Seasonal profiles (60-day windows)") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

fig <- p_dur + p_seas + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("figures/resolution_sensitivity.pdf", fig, width = 9, height = 4.5)
cat("\nSaved figures/resolution_sensitivity.pdf\n")
cat("Saved resolution_sensitivity_summary.csv\n")
cat("\nDone.\n")

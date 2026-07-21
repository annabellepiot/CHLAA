# =========================================================================
# Collect and Compare Fitted Parameters Across Health Zones
# =========================================================================
#
# This script collects posterior estimates from all successful fits and
# generates comparative summaries and visualizations. Useful for identifying
# parameter variation across geographies and detecting outliers.
#
# Reads fit artifacts from: figures/.rds files/*_fit.rds
# Outputs tables to:        figures/tables/
# Outputs figures to:       figures/
#
# Usage: Rscript 01_05_collect_fitted_parameters.R
#
# =========================================================================

library(tidyverse)
library(chlaa)

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- Setup ----

rds_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/figures/.rds files"
fig_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/figures"
tab_dir <- file.path(fig_dir, "tables")

dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# R0 scaling constant: K_R0 = POP_REF * 5.1446e-3 / H_REF (must match pipeline)
K_R0 <- 2654.6

# ---- Get list of successful fits ----

fit_files <- list.files(rds_dir, pattern = "_fit\\.rds$", full.names = TRUE)

if (length(fit_files) == 0) {
  stop("No successful fits found in: ", rds_dir)
}

cat("\n", rep("=", 70), "\n", sep = "")
cat("Collecting Fitted Parameters\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Found", length(fit_files), "fit artifacts.\n\n")

# ---- Extract data from each artifact ----

all_posteriors  <- list()
all_traces      <- list()
all_diagnostics <- list()
all_r0          <- list()
all_budgets     <- list()

for (fit_file in fit_files) {
  art <- readRDS(fit_file)
  hz  <- art$hz_name

  # --- Posterior summary (log/logit scale from report) ---
  post <- art$report$posterior_summary
  post$hz           <- hz
  post$total_cases  <- art$total_cases
  post$n_weeks      <- art$n_weeks
  post$acceptance_rate <- art$report$acceptance_rate

  # --- Natural-scale posterior from trace (+ R0 per-draw) ---
  tr <- chlaa_fit_trace(art$fit, burnin = 0.25, scale = "natural")
  tr_wide <- tr |>
    pivot_wider(names_from = parameter, values_from = value)
  tr_wide$R0 <- tr_wide$frac_neff * K_R0 * tr_wide$trans_prob

  # Pivot back to long for summaries and plotting
  draws_long <- tr_wide |>
    select(trans_prob, frac_neff, obs_size, E0, R0) |>
    pivot_longer(everything(), names_to = "parameter", values_to = "value")

  nat_summary <- draws_long |>
    group_by(parameter) |>
    summarise(
      median = median(value),
      q025   = quantile(value, 0.025),
      q975   = quantile(value, 0.975),
      mean   = mean(value),
      sd     = sd(value),
      .groups = "drop"
    )
  nat_summary$hz           <- hz
  nat_summary$total_cases  <- art$total_cases
  nat_summary$n_weeks      <- art$n_weeks
  nat_summary$acceptance_rate <- art$report$acceptance_rate

  all_posteriors[[hz]] <- nat_summary

  # --- Raw trace draws (for density plots) ---
  draws_long$hz <- hz
  all_traces[[hz]] <- draws_long

  # --- Pre-computed tables ---
  if (!is.null(art$diagnostics)) all_diagnostics[[hz]] <- art$diagnostics
  if (!is.null(art$r0_table))    all_r0[[hz]]          <- art$r0_table
  if (!is.null(art$budget))      all_budgets[[hz]]     <- art$budget

  cat(sprintf("  %-20s  cases=%5d  weeks=%3d  accept=%.3f\n",
              hz, art$total_cases, art$n_weeks, art$report$acceptance_rate))
}

# Combine
all_params <- bind_rows(all_posteriors) |>
  relocate(hz, total_cases, n_weeks, acceptance_rate)

all_draws <- bind_rows(all_traces)

diag_table <- bind_rows(all_diagnostics)
r0_table   <- bind_rows(all_r0)
budget_table <- bind_rows(all_budgets)

cat("\nSuccessfully extracted parameters from", length(all_posteriors), "health zones.\n\n")

# ---- Summary statistics by parameter ----

cat(rep("=", 70), "\n", sep = "")
cat("Parameter Summary Statistics (natural scale)\n")
cat(rep("=", 70), "\n\n", sep = "")

param_stats <- all_params |>
  group_by(parameter) |>
  summarise(
    n_hzs      = n(),
    mean_med   = mean(median),
    sd_med     = sd(median),
    min_med    = min(median),
    max_med    = max(median),
    mean_q025  = mean(q025),
    mean_q975  = mean(q975),
    .groups = "drop"
  )

print(param_stats, n = Inf)
cat("\n")

# ---- Parameter estimates by HZ (wide format) ----

cat(rep("=", 70), "\n", sep = "")
cat("Parameter Estimates by Health Zone\n")
cat(rep("=", 70), "\n\n", sep = "")

param_wide <- all_params |>
  select(hz, parameter, median, q025, q975) |>
  pivot_wider(
    names_from  = parameter,
    values_from = c(median, q025, q975),
    names_glue  = "{parameter}_{.value}"
  )

print(param_wide, n = Inf)
cat("\n")

# ---- Identify outliers ----

cat(rep("=", 70), "\n", sep = "")
cat("Potential Outliers (>2 SD from mean)\n")
cat(rep("=", 70), "\n\n", sep = "")

outliers_found <- FALSE

for (param in unique(all_params$parameter)) {
  param_data <- all_params |> filter(parameter == !!param)

  mean_val <- mean(param_data$median)
  sd_val   <- sd(param_data$median)

  outliers <- param_data |> filter(abs(median - mean_val) > 2 * sd_val)

  if (nrow(outliers) > 0) {
    outliers_found <- TRUE
    cat(sprintf("%s:\n", param))
    for (i in 1:nrow(outliers)) {
      cat(sprintf("  %s: %.4f (mean=%.4f, sd=%.4f)\n",
                  outliers$hz[i], outliers$median[i], mean_val, sd_val))
    }
    cat("\n")
  }
}

if (!outliers_found) cat("No outliers detected.\n\n")

# ---- Correlation between outbreak size and parameters ----

cat(rep("=", 70), "\n", sep = "")
cat("Correlation Between Outbreak Size and Parameters\n")
cat(rep("=", 70), "\n\n", sep = "")

correlations <- all_params |>
  group_by(parameter) |>
  summarise(
    cor_with_cases = cor(median, total_cases, use = "complete.obs"),
    cor_with_weeks = cor(median, n_weeks,     use = "complete.obs"),
    .groups = "drop"
  )

print(correlations, n = Inf)
cat("\n")

# ---- R0 summary (from pre-computed r0_table) ----

if (nrow(r0_table) > 0) {
  cat(rep("=", 70), "\n", sep = "")
  cat("R0 Estimates Across Health Zones\n")
  cat(rep("=", 70), "\n\n", sep = "")
  print(r0_table |> select(hz, pop, R0_med, R0_lo, R0_hi, N_eff_med), n = Inf)
  cat("\n")
}

# ---- Export tables ----

write.csv(all_params,
  file.path(tab_dir, "fitted_parameters_summary.csv"),
  row.names = FALSE)

write.csv(param_wide,
  file.path(tab_dir, "fitted_parameters_wide.csv"),
  row.names = FALSE)

if (nrow(r0_table) > 0) {
  write.csv(r0_table,
    file.path(tab_dir, "r0_estimates_all_hzs.csv"),
    row.names = FALSE)
}

if (nrow(diag_table) > 0) {
  write.csv(diag_table,
    file.path(tab_dir, "diagnostics_all_hzs.csv"),
    row.names = FALSE)
}

cat("Exported tables to:", tab_dir, "\n\n")

# =========================================================================
# VISUALIZATIONS
# =========================================================================

# ---- 1. Posterior distributions across HZs (5 params, density curves) ----
#
# ggh4x::facetted_pos_scales is incompatible with ggdist::stat_halfeye
# (produces a blank plot), so we split into two sub-plots:
#   - linear-scale:  R0, frac_neff, obs_size   (3 panels)
#   - log10-scale:   trans_prob, E0             (2 panels)
# Combined with patchwork.

library(ggdist)
library(patchwork)

param_labels <- c(
  R0         = "Basic reproduction number (R0)",
  trans_prob = "Transmission probability",
  frac_neff  = "Effective population fraction",
  obs_size   = "Observation size",
  E0         = "Initial seed (E0)"
)

# Order HZs by median R0
hz_order <- all_draws |>
  filter(parameter == "R0") |>
  group_by(hz) |>
  summarise(med = median(value), .groups = "drop") |>
  arrange(med) |>
  pull(hz)

draws_all <- all_draws |>
  filter(parameter %in% names(param_labels)) |>
  mutate(hz = factor(hz, levels = hz_order))

# Per-facet summary labels (median [95% CrI])
summ <- draws_all |>
  group_by(hz, parameter) |>
  ggdist::median_qi(value, .width = 0.95) |>
  ungroup()

# Shared theme
theme_halfeye <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    axis.text.y      = element_text(size = 9),
    plot.title       = element_text(face = "bold")
  )

# --- Standalone R0 plot ---
draws_r0 <- draws_all |> filter(parameter == "R0")
summ_r0  <- summ |> filter(parameter == "R0")

p_r0_density <- ggplot(draws_r0, aes(x = value, y = hz)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey40") +
  stat_halfeye(
    .width = c(0.8, 0.95), fill = "#6baed6",
    normalize = "panels", slab_alpha = 0.8,
    point_size = 1.5, interval_size_range = c(0.6, 1.2)
  ) +
  geom_text(
    data = summ_r0,
    aes(x = Inf, y = hz,
        label = sprintf("%.2f [%.2f, %.2f]", value, .lower, .upper)),
    hjust = "inward", size = 3.0, colour = "grey30"
  ) +
  labs(
    x        = expression(R[0]),
    y        = NULL,
    title    = "R0 estimates by health zone",
    subtitle = "Densities with median, 80% and 95% credible intervals"
  ) +
  theme_halfeye

ggsave(file.path(fig_dir, "fitted_r0_density.png"),
  p_r0_density, width = 8, height = 7, dpi = 300)

# --- 4-parameter comparison (2x2): linear pair + log pair via patchwork ---
lin_params  <- c("frac_neff", "obs_size")
draws_lin   <- draws_all |>
  filter(parameter %in% lin_params) |>
  mutate(parameter = factor(parameter, levels = lin_params))
summ_lin    <- summ |>
  filter(parameter %in% lin_params) |>
  mutate(parameter = factor(parameter, levels = lin_params))

p_lin <- ggplot(draws_lin, aes(x = value, y = hz)) +
  stat_halfeye(
    .width = c(0.8, 0.95), fill = "#6baed6",
    normalize = "panels", slab_alpha = 0.8,
    point_size = 1.5, interval_size_range = c(0.6, 1.2)
  ) +
  geom_text(
    data = summ_lin,
    aes(x = Inf, y = hz,
        label = sprintf("%.3g [%.3g, %.3g]", value, .lower, .upper)),
    hjust = "inward", size = 2.0, colour = "grey30"
  ) +
  facet_wrap(~ parameter, scales = "free_x", nrow = 1,
             labeller = labeller(parameter = param_labels)) +
  labs(x = "Posterior value", y = NULL) +
  theme_halfeye

log_params  <- c("trans_prob", "E0")
draws_log   <- draws_all |>
  filter(parameter %in% log_params) |>
  mutate(parameter = factor(parameter, levels = log_params))
summ_log    <- summ |>
  filter(parameter %in% log_params) |>
  mutate(parameter = factor(parameter, levels = log_params))

p_log <- ggplot(draws_log, aes(x = value, y = hz)) +
  stat_halfeye(
    .width = c(0.8, 0.95), fill = "#6baed6",
    normalize = "panels", slab_alpha = 0.8,
    point_size = 1.5, interval_size_range = c(0.6, 1.2)
  ) +
  geom_text(
    data = summ_log,
    aes(x = Inf, y = hz,
        label = sprintf("%.3g [%.3g, %.3g]", value, .lower, .upper)),
    hjust = "inward", size = 2.0, colour = "grey30"
  ) +
  facet_wrap(~ parameter, scales = "free_x", nrow = 1,
             labeller = labeller(parameter = param_labels)) +
  scale_x_log10() +
  labs(x = "Posterior value (log scale)", y = NULL) +
  theme_halfeye

# Stack linear (top) and log (bottom) into a 2x2 grid
p_dist <- (p_lin / p_log) +
  plot_annotation(
    title    = "Posterior parameter estimates by health zone",
    subtitle = "Half-eye densities with median, 80% and 95% credible intervals",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11)
    )
  )

ggsave(file.path(fig_dir, "fitted_parameters_comparison.png"),
  p_dist, width = 16, height = 12, dpi = 300)

# ---- 2. R0 estimates across HZs ----

if (nrow(r0_table) > 0) {
  p_r0 <- ggplot(r0_table, aes(x = R0_med, y = reorder(hz, R0_med))) +
    geom_point(size = 3, colour = "#2f6a4e") +
    geom_errorbarh(aes(xmin = R0_lo, xmax = R0_hi),
      height = 0.2, alpha = 0.6, colour = "#2f6a4e") +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "red", linewidth = 0.5) +
    labs(
      title    = "Estimated R0 Across Health Zones",
      subtitle = "Points show posterior medians; error bars show 95% credible intervals",
      x = expression(R[0]),
      y = "Health Zone"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(fig_dir, "fitted_r0_comparison.png"),
    p_r0, width = 10, height = 8, dpi = 300)
}

# ---- 3. trans_prob vs frac_neff (the identifiability ridge) ----

param_wide_plot <- all_params |>
  select(hz, parameter, median, total_cases) |>
  pivot_wider(names_from = parameter, values_from = median)

if (all(c("trans_prob", "frac_neff") %in% colnames(param_wide_plot))) {
  p_ridge <- ggplot(param_wide_plot, aes(x = trans_prob, y = frac_neff)) +
    geom_point(aes(size = total_cases, colour = total_cases)) +
    geom_text(aes(label = hz), hjust = -0.1, vjust = -0.1, size = 3) +
    scale_colour_viridis_c(trans = "log10") +
    scale_size_continuous(range = c(3, 10)) +
    labs(
      title  = "Transmission Probability vs Effective Population Fraction",
      subtitle = "Identifiability ridge: trans_prob and frac_neff are negatively correlated",
      x = "Transmission Probability (trans_prob)",
      y = "Effective Population Fraction (frac_neff)",
      colour = "Total Cases",
      size   = "Total Cases"
    ) +
    theme_bw(base_size = 12)

  ggsave(file.path(fig_dir, "fitted_params_trans_vs_frac_neff.png"),
    p_ridge, width = 10, height = 8, dpi = 300)
}

# ---- 4. obs_size vs total_cases ----

if ("obs_size" %in% colnames(param_wide_plot)) {
  p_obs <- ggplot(param_wide_plot, aes(x = total_cases, y = obs_size)) +
    geom_point(aes(colour = frac_neff), size = 4) +
    geom_text(aes(label = hz), hjust = -0.1, vjust = -0.1, size = 3) +
    scale_x_log10() +
    scale_colour_viridis_c() +
    labs(
      title  = "Observation Overdispersion vs Outbreak Size",
      x = "Total Cases (log scale)",
      y = "Observation Size (obs_size)",
      colour = "frac_neff"
    ) +
    theme_bw(base_size = 12)

  ggsave(file.path(fig_dir, "fitted_params_obssize_vs_cases.png"),
    p_obs, width = 10, height = 8, dpi = 300)
}

# ---- 5. Acceptance rates ----

accept_data <- all_params |>
  distinct(hz, acceptance_rate, total_cases)

p_accept <- ggplot(accept_data, aes(x = reorder(hz, acceptance_rate), y = acceptance_rate)) +
  geom_col(aes(fill = total_cases)) +
  geom_hline(yintercept = 0.15, linetype = "dashed", colour = "red", linewidth = 0.8) +
  scale_fill_viridis_c(trans = "log10", name = "Total\nCases") +
  labs(
    title    = "MCMC Acceptance Rates by Health Zone",
    subtitle = "Red dashed line indicates 0.15 threshold",
    x = "Health Zone",
    y = "Acceptance Rate"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "fitted_acceptance_rates.png"),
  p_accept, width = 10, height = 6, dpi = 300)

# ---- 6. Convergence diagnostics (R-hat) ----

if (nrow(diag_table) > 0) {
  rhat_cols <- diag_table |>
    select(hz, starts_with("rhat_")) |>
    pivot_longer(-hz, names_to = "parameter", values_to = "rhat") |>
    mutate(parameter = str_remove(parameter, "^rhat_"))

  p_rhat <- ggplot(rhat_cols, aes(x = rhat, y = reorder(hz, rhat))) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = 1.1, linetype = "dashed", colour = "red", linewidth = 0.5) +
    facet_wrap(~parameter, scales = "free_x", ncol = 2) +
    labs(
      title    = "R-hat Convergence Diagnostics",
      subtitle = "Values > 1.1 (red dashed) suggest poor convergence",
      x = expression(hat(R)),
      y = "Health Zone"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )

  ggsave(file.path(fig_dir, "fitted_rhat_diagnostics.png"),
    p_rhat, width = 12, height = 10, dpi = 300)
}

# ---- Summary ----

cat("\nSaved figures:\n")
cat("  ", file.path(fig_dir, "fitted_r0_density.png"), "\n")
cat("  ", file.path(fig_dir, "fitted_parameters_comparison.png"), "\n")
cat("  ", file.path(fig_dir, "fitted_r0_comparison.png"), "\n")
cat("  ", file.path(fig_dir, "fitted_params_trans_vs_frac_neff.png"), "\n")
cat("  ", file.path(fig_dir, "fitted_params_obssize_vs_cases.png"), "\n")
cat("  ", file.path(fig_dir, "fitted_acceptance_rates.png"), "\n")
cat("  ", file.path(fig_dir, "fitted_rhat_diagnostics.png"), "\n")

cat("\n", rep("=", 70), "\n", sep = "")
cat("Parameter collection complete!\n")
cat(rep("=", 70), "\n\n", sep = "")

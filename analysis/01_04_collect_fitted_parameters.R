# =========================================================================
# Collect and Compare Fitted Parameters Across Health Zones
# =========================================================================
#
# This script collects posterior estimates from all successful fits and
# generates comparative summaries and visualizations. Useful for identifying
# parameter variation across geographies and detecting outliers.
#
# Usage: Rscript collect_fitted_parameters.R
#
# =========================================================================

library(tidyverse)
library(ggplot2)

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- Setup ----

data_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/analysis/data"
output_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/output"
fig_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/figures"

# ---- Get list of successful fits ----

fit_files <- list.files(output_dir, pattern = "_fit\\.rds$", full.names = TRUE)

if (length(fit_files) == 0) {
  stop("No successful fits found in output directory.")
}

cat("\n", rep("=", 70), "\n", sep = "")
cat("Collecting Fitted Parameters\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Found", length(fit_files), "successful fits.\n\n")

# ---- Extract posterior summaries ----

parameter_summaries <- list()

for (fit_file in fit_files) {
  fit_data <- readRDS(fit_file)
  hz_name <- fit_data$hz_name

  post_summary <- fit_data$parameter_summary %||% fit_data$report$posterior_summary

  # Add metadata
  post_summary$hz <- hz_name
  post_summary$total_cases <- fit_data$total_cases
  post_summary$n_weeks <- fit_data$n_weeks
  post_summary$acceptance_rate <- fit_data$report$acceptance_rate

  parameter_summaries[[hz_name]] <- post_summary
}

# Combine all summaries
all_params <- bind_rows(parameter_summaries) %>%
  relocate(hz, total_cases, n_weeks, acceptance_rate)

cat("Successfully extracted parameters from", nrow(all_params) / 4, "health zones.\n\n")

# ---- Summary statistics by parameter ----

cat(rep("=", 70), "\n", sep = "")
cat("Parameter Summary Statistics\n")
cat(rep("=", 70), "\n\n", sep = "")

param_stats <- all_params %>%
  group_by(parameter) %>%
  summarise(
    n_hzs = n(),
    mean_median = mean(median),
    sd_median = sd(median),
    min_median = min(median),
    max_median = max(median),
    mean_q025 = mean(q025),
    mean_q975 = mean(q975),
    .groups = "drop"
  )

print(param_stats, n = Inf)
cat("\n")

# ---- Parameter estimates by HZ ----

cat(rep("=", 70), "\n", sep = "")
cat("Parameter Estimates by Health Zone\n")
cat(rep("=", 70), "\n\n", sep = "")

param_wide <- all_params %>%
  select(hz, parameter, median, q025, q975) %>%
  pivot_wider(
    names_from = parameter,
    values_from = c(median, q025, q975),
    names_glue = "{parameter}_{.value}"
  )

print(param_wide, n = Inf)
cat("\n")

# ---- Identify outliers ----

cat(rep("=", 70), "\n", sep = "")
cat("Potential Outliers (>2 SD from mean)\n")
cat(rep("=", 70), "\n\n", sep = "")

outliers_found <- FALSE

for (param in unique(all_params$parameter)) {
  param_data <- all_params %>% filter(parameter == !!param)

  mean_val <- mean(param_data$median)
  sd_val <- sd(param_data$median)

  outliers <- param_data %>%
    filter(abs(median - mean_val) > 2 * sd_val)

  if (nrow(outliers) > 0) {
    outliers_found <- TRUE
    cat(sprintf("%s:\n", param))
    for (i in 1:nrow(outliers)) {
      cat(sprintf(
        "  %s: %.4f (mean=%.4f, sd=%.4f)\n",
        outliers$hz[i], outliers$median[i], mean_val, sd_val
      ))
    }
    cat("\n")
  }
}

if (!outliers_found) {
  cat("No outliers detected.\n\n")
}

# ---- Correlation between outbreak size and parameters ----

cat(rep("=", 70), "\n", sep = "")
cat("Correlation Between Outbreak Size and Parameters\n")
cat(rep("=", 70), "\n\n", sep = "")

correlations <- all_params %>%
  group_by(parameter) %>%
  summarise(
    cor_with_cases = cor(median, total_cases, use = "complete.obs"),
    cor_with_weeks = cor(median, n_weeks, use = "complete.obs"),
    .groups = "drop"
  )

print(correlations, n = Inf)
cat("\n")

# ---- Export parameter table ----

write.csv(all_params,
  file.path(output_dir, "fitted_parameters_summary.csv"),
  row.names = FALSE
)

write.csv(param_wide,
  file.path(output_dir, "fitted_parameters_wide.csv"),
  row.names = FALSE
)

cat(rep("=", 70), "\n", sep = "")
cat("Exported tables:\n")
cat("  ", file.path(output_dir, "fitted_parameters_summary.csv"), "\n")
cat("  ", file.path(output_dir, "fitted_parameters_wide.csv"), "\n")
cat(rep("=", 70), "\n\n", sep = "")

# ---- Visualization: Parameter distributions ----

p_dist <- ggplot(all_params, aes(x = median, y = reorder(hz, median))) +
  geom_point(aes(colour = total_cases), size = 3) +
  geom_errorbarh(aes(xmin = q025, xmax = q975, colour = total_cases),
    height = 0.2, alpha = 0.6
  ) +
  facet_wrap(~parameter, scales = "free_x", ncol = 2) +
  scale_colour_viridis_c(name = "Total\nCases", trans = "log10") +
  labs(
    title = "Fitted Parameters Across Health Zones",
    subtitle = "Points show posterior medians; error bars show 95% credible intervals",
    x = "Parameter Value",
    y = "Health Zone"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fitted_parameters_comparison.png"),
  p_dist,
  width = 12, height = 10, dpi = 300
)

cat("Saved figure:", file.path(fig_dir, "fitted_parameters_comparison.png"), "\n\n")

# ---- Visualization: Parameter correlations ----

param_wide_plot <- all_params %>%
  select(hz, parameter, median, total_cases) %>%
  pivot_wider(names_from = parameter, values_from = median)

# trans_prob vs reporting_rate
p_cor1 <- ggplot(param_wide_plot, aes(x = trans_prob, y = reporting_rate)) +
  geom_point(aes(size = total_cases, colour = total_cases)) +
  geom_text(aes(label = hz), hjust = -0.1, vjust = -0.1, size = 3) +
  scale_colour_viridis_c(trans = "log10") +
  scale_size_continuous(range = c(3, 10)) +
  labs(
    title = "Transmission Probability vs Reporting Rate",
    x = "Transmission Probability (trans_prob)",
    y = "Reporting Rate",
    colour = "Total Cases",
    size = "Total Cases"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(fig_dir, "fitted_params_trans_vs_report.png"),
  p_cor1,
  width = 10, height = 8, dpi = 300
)

# obs_size vs total_cases
p_cor2 <- ggplot(param_wide_plot, aes(x = total_cases, y = obs_size)) +
  geom_point(aes(colour = reporting_rate), size = 4) +
  geom_text(aes(label = hz), hjust = -0.1, vjust = -0.1, size = 3) +
  scale_x_log10() +
  scale_colour_viridis_c() +
  labs(
    title = "Observation Overdispersion vs Outbreak Size",
    x = "Total Cases (log scale)",
    y = "Observation Size (obs_size)",
    colour = "Reporting\nRate"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(fig_dir, "fitted_params_obssize_vs_cases.png"),
  p_cor2,
  width = 10, height = 8, dpi = 300
)

cat("Saved figures:\n")
cat("  ", file.path(fig_dir, "fitted_params_trans_vs_report.png"), "\n")
cat("  ", file.path(fig_dir, "fitted_params_obssize_vs_cases.png"), "\n\n")

# ---- Visualization: Acceptance rates ----

p_accept <- all_params %>%
  distinct(hz, acceptance_rate, total_cases) %>%
  ggplot(aes(x = reorder(hz, acceptance_rate), y = acceptance_rate)) +
  geom_col(aes(fill = total_cases)) +
  geom_hline(yintercept = 0.15, linetype = "dashed", colour = "red", linewidth = 0.8) +
  scale_fill_viridis_c(trans = "log10", name = "Total\nCases") +
  labs(
    title = "MCMC Acceptance Rates by Health Zone",
    subtitle = "Red dashed line indicates 0.15 threshold (below may indicate identifiability issues)",
    x = "Health Zone",
    y = "Acceptance Rate"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "fitted_acceptance_rates.png"),
  p_accept,
  width = 10, height = 6, dpi = 300
)

cat("Saved figure:", file.path(fig_dir, "fitted_acceptance_rates.png"), "\n\n")

cat(rep("=", 70), "\n", sep = "")
cat("Parameter collection complete!\n")
cat(rep("=", 70), "\n\n", sep = "")

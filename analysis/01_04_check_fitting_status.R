# =========================================================================
# Check Fitting Status Across All Health Zones
# =========================================================================
#
# This script checks the status of PMCMC fitting jobs across all health zones.
# It reports which HZs have completed successfully, which failed, and which
# are still pending.
#
# Usage: Rscript check_fitting_status.R
#
# =========================================================================

library(tidyverse)

# ---- Setup ----

data_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/analysis/data"
output_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/output"

# ---- Get list of expected HZs ----

hz_params_long <- read.csv(file.path(data_dir, "hz_parameters.csv"),
                            stringsAsFactors = FALSE)
all_hzs <- sort(unique(hz_params_long$hz))

cat("\n", rep("=", 70), "\n", sep = "")
cat("PMCMC Fitting Status Report\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Total health zones:", length(all_hzs), "\n")
cat("Health zones:", paste(all_hzs, collapse = ", "), "\n\n")

# ---- Check status for each HZ ----

status <- data.frame(
  hz = all_hzs,
  status = character(length(all_hzs)),
  acceptance_rate = numeric(length(all_hzs)),
  total_cases = integer(length(all_hzs)),
  n_weeks = integer(length(all_hzs)),
  timestamp = character(length(all_hzs)),
  stringsAsFactors = FALSE
)

for (i in seq_along(all_hzs)) {
  hz_name <- all_hzs[i]

  # Check for successful fit
  success_file <- file.path(output_dir, sprintf("%s_fit.rds", hz_name))
  failed_file <- file.path(output_dir, sprintf("%s_FAILED.rds", hz_name))

  if (file.exists(success_file)) {
    fit_data <- readRDS(success_file)
    status$status[i] <- "SUCCESS"
    status$acceptance_rate[i] <- fit_data$report$acceptance_rate
    status$total_cases[i] <- fit_data$total_cases
    status$n_weeks[i] <- fit_data$n_weeks
    status$timestamp[i] <- as.character(fit_data$timestamp)

  } else if (file.exists(failed_file)) {
    fail_data <- readRDS(failed_file)
    status$status[i] <- "FAILED"
    status$acceptance_rate[i] <- NA
    status$total_cases[i] <- NA
    status$n_weeks[i] <- NA
    status$timestamp[i] <- as.character(fail_data$timestamp)

  } else {
    status$status[i] <- "PENDING"
    status$acceptance_rate[i] <- NA
    status$total_cases[i] <- NA
    status$n_weeks[i] <- NA
    status$timestamp[i] <- NA
  }
}

# ---- Summary statistics ----

n_success <- sum(status$status == "SUCCESS")
n_failed <- sum(status$status == "FAILED")
n_pending <- sum(status$status == "PENDING")

cat(rep("=", 70), "\n", sep = "")
cat("Summary\n")
cat(rep("=", 70), "\n", sep = "")
cat(sprintf("Successful: %d / %d (%.1f%%)\n", n_success, length(all_hzs), 100 * n_success / length(all_hzs)))
cat(sprintf("Failed:     %d / %d (%.1f%%)\n", n_failed, length(all_hzs), 100 * n_failed / length(all_hzs)))
cat(sprintf("Pending:    %d / %d (%.1f%%)\n", n_pending, length(all_hzs), 100 * n_pending / length(all_hzs)))
cat("\n")

# ---- Detailed status ----

cat(rep("=", 70), "\n", sep = "")
cat("Detailed Status\n")
cat(rep("=", 70), "\n\n", sep = "")

# Successful fits
if (n_success > 0) {
  cat("SUCCESSFUL FITS:\n")
  cat(rep("-", 70), "\n", sep = "")

  success_df <- status %>%
    filter(status == "SUCCESS") %>%
    arrange(hz)

  for (i in 1:nrow(success_df)) {
    cat(sprintf("  %s:\n", success_df$hz[i]))
    cat(sprintf("    Acceptance rate: %.3f\n", success_df$acceptance_rate[i]))
    cat(sprintf("    Total cases: %d\n", success_df$total_cases[i]))
    cat(sprintf("    Outbreak weeks: %d\n", success_df$n_weeks[i]))
    cat(sprintf("    Completed: %s\n", success_df$timestamp[i]))
  }
  cat("\n")
}

# Failed fits
if (n_failed > 0) {
  cat("FAILED FITS:\n")
  cat(rep("-", 70), "\n", sep = "")

  failed_df <- status %>%
    filter(status == "FAILED") %>%
    arrange(hz)

  for (i in 1:nrow(failed_df)) {
    cat(sprintf("  %s:\n", failed_df$hz[i]))

    # Try to get error message
    failed_file <- file.path(output_dir, sprintf("%s_FAILED.rds", failed_df$hz[i]))
    fail_data <- readRDS(failed_file)

    if (!is.null(fail_data$error)) {
      cat(sprintf("    Error: %s\n", fail_data$error))
    }
    cat(sprintf("    Failed at: %s\n", failed_df$timestamp[i]))
  }
  cat("\n")
}

# Pending fits
if (n_pending > 0) {
  cat("PENDING FITS:\n")
  cat(rep("-", 70), "\n", sep = "")

  pending_df <- status %>%
    filter(status == "PENDING") %>%
    arrange(hz)

  cat(sprintf("  %s\n", paste(pending_df$hz, collapse = ", ")))
  cat("\n")
}

# ---- Acceptance rate diagnostics ----

if (n_success > 0) {
  cat(rep("=", 70), "\n", sep = "")
  cat("Acceptance Rate Diagnostics\n")
  cat(rep("=", 70), "\n", sep = "")

  success_df <- status %>%
    filter(status == "SUCCESS") %>%
    arrange(acceptance_rate)

  cat(sprintf("Mean acceptance rate: %.3f\n", mean(success_df$acceptance_rate, na.rm = TRUE)))
  cat(sprintf("Median acceptance rate: %.3f\n", median(success_df$acceptance_rate, na.rm = TRUE)))
  cat(sprintf("Range: [%.3f, %.3f]\n\n",
              min(success_df$acceptance_rate, na.rm = TRUE),
              max(success_df$acceptance_rate, na.rm = TRUE)))

  # Flag low acceptance rates
  low_accept <- success_df %>% filter(acceptance_rate < 0.15)

  if (nrow(low_accept) > 0) {
    cat("WARNING: Low acceptance rates (<0.15) detected:\n")
    for (i in 1:nrow(low_accept)) {
      cat(sprintf("  %s: %.3f (may indicate identifiability issues)\n",
                  low_accept$hz[i], low_accept$acceptance_rate[i]))
    }
    cat("\n")
  }
}

# ---- Export status table ----

write.csv(status, file.path(output_dir, "fitting_status_summary.csv"), row.names = FALSE)

cat(rep("=", 70), "\n", sep = "")
cat("Status table saved to:", file.path(output_dir, "fitting_status_summary.csv"), "\n")
cat(rep("=", 70), "\n\n", sep = "")

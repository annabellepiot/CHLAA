#=========================================================================
#Collect and Compare Fitted Parameters Across Health Zones
#=========================================================================
#
#This script collects posterior estimates from all successful fits and
#generates comparative summaries and visualizations. Useful for identifying
#parameter variation across geographies and detecting outliers.
#
#Reads fit artifacts from: figures/.rds files/*_fit.rds
#Outputs tables to:        figures/tables/
#Outputs figures to:       figures/
#
#Usage: Rscript 01_05_collect_fitted_parameters.R
#
#=========================================================================

library(tidyverse)
library(chlaa)

`%||%` <- function(x, y) if (is.null(x)) y else x

#---- Setup ----

rds_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/figures/.rds files"
fig_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/figures"
tab_dir <- file.path(fig_dir, "tables")

dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

#---- Get list of successful fits ----

fit_files <- list.files(rds_dir, pattern = "_fit\\.rds$", full.names = TRUE)

if (length(fit_files) == 0) {
  stop("No successful fits found in: ", rds_dir)
}

cat("\n", rep("=", 70), "\n", sep = "")
cat("Collecting Fitted Parameters\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Found", length(fit_files), "fit artifacts.\n\n")

#---- Extract data from each artifact ----

all_posteriors  <- list()
all_traces      <- list()
all_diagnostics <- list()
all_r0          <- list()
all_budgets     <- list()

for (fit_file in fit_files) {
  art <- readRDS(fit_file)
  hz  <- art$hz_name

  #--- Posterior summary (log/logit scale from report) ---
  post <- art$report$posterior_summary
  post$hz           <- hz
  post$total_cases  <- art$total_cases
  post$n_weeks      <- art$n_weeks
  post$acceptance_rate <- art$report$acceptance_rate

  #--- Natural-scale posterior from trace (+ R0 per-draw) ---
  tr <- chlaa_fit_trace(art$fit, burnin = 0.25, scale = "natural")
  tr_wide <- tr |>
    pivot_wider(names_from = parameter, values_from = value)
  pars_for_r0 <- art$pars_warm
  pars_for_r0$trans_prob <- tr_wide$trans_prob
  pars_for_r0$N <- tr_wide$N
  tr_wide$R0 <- chlaa_r0(pars_for_r0)

  #Pivot back to long for summaries and plotting
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

  #--- Raw trace draws (for density plots) ---
  draws_long$hz <- hz
  all_traces[[hz]] <- draws_long

  #--- Pre-computed tables ---
  if (!is.null(art$diagnostics)) all_diagnostics[[hz]] <- art$diagnostics
  if (!is.null(art$r0_table))    all_r0[[hz]]          <- art$r0_table
  if (!is.null(art$budget))      all_budgets[[hz]]     <- art$budget

  cat(sprintf("  %-20s  cases=%5d  weeks=%3d  accept=%.3f\n",
              hz, art$total_cases, art$n_weeks, art$report$acceptance_rate))
}

#Combine
all_params <- bind_rows(all_posteriors) |>
  relocate(hz, total_cases, n_weeks, acceptance_rate)

all_draws <- bind_rows(all_traces)

diag_table <- bind_rows(all_diagnostics)
r0_table   <- bind_rows(all_r0)
budget_table <- bind_rows(all_budgets)

cat("\nSuccessfully extracted parameters from", length(all_posteriors), "health zones.\n\n")

#---- Summary statistics by parameter ----

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

#---- Parameter estimates by HZ (wide format) ----

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

#---- Identify outliers ----

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

#---- Correlation between outbreak size and parameters ----

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

#---- R0 summary (from pre-computed r0_table) ----

if (nrow(r0_table) > 0) {
  cat(rep("=", 70), "\n", sep = "")
  cat("R0 Estimates Across Health Zones\n")
  cat(rep("=", 70), "\n\n", sep = "")
  print(r0_table |> select(hz, pop, R0_med, R0_lo, R0_hi, N_eff_med), n = Inf)
  cat("\n")
}

#---- Export tables ----

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

#=========================================================================
#GOODNESS OF FIT  -  CRPS (free-running posterior predictive)
#=========================================================================
#
#For each health zone: draw posterior parameter samples, simulate the model
#FREELY from its seeded initial state over the observed weeks (no filtering,
#no conditioning on later data), apply the fitted observation model to get
#predicted REPORTED cases, and score those predictions against the observed
#weekly counts with the CRPS.
#
#This is a GOODNESS-OF-FIT / trajectory-replication measure, NOT forecast
#skill (the weeks scored are the weeks the model was fitted to). It is the
#free-running predictive (not a one-step-ahead filtered predictive) because
#the scenario, Shapley and cost-effectiveness analyses all run the model free
#from t = 0 - so this validates the mode in which the model is actually used.
#
#Observation model (matches fitting, see chlaa fit.R):
#   cases ~ NegBin(mu = reporting_rate * inc_symptoms_weekly, size = obs_size)
#
#HEADLINE NUMBER: crps_log - mean CRPS on the log1p scale (Bosse et al. 2023,
#PLoS Comput Biol). Raw CRPS is in weekly-case units so it scales with
#outbreak size and is NOT comparable across zones (Goma dwarfs Bumbu). The log
#transform makes the expected score approximately magnitude-independent, so
#crps_log behaves like a probabilistic relative error and IS comparable.
#Lower is better; expm1(crps_log) ~ typical relative error.

stopifnot(requireNamespace("scoringRules", quietly = TRUE))
stopifnot(requireNamespace("posterior",   quietly = TRUE))

#---- Settings ----
N_PRED <- 1000   #predictive samples per week (>=1000: the sample CRPS
                 #estimator is biased at small n). NOTE: drawn from an
                 #autocorrelated chain, so the effective sample size is
                 #smaller. Runtime is ~N_PRED simulations x ~12 zones - this
                 #is the tunable knob if the section is too slow.
BURNIN <- 0.25   #matches chlaa_fit_trace(burnin = 0.25) used above
DT     <- 0.25   #MUST match the dt used by the particle filter at fit time
SEED   <- 202

#---- Predictive samples: free-running simulation ----
#Uses the same engine helpers as chlaa's own predictive (chlaa_forecast_from_fit).
crps_predictive_cases <- function(art, n = N_PRED, seed = SEED) {
  dr <- chlaa:::.chlaa_fit_selected_draws_matrix(art$fit, burnin = BURNIN, thin = 1)
  tvec <- art$observed$time
  #weekly data -> inc_symptoms_weekly (robust to obs_interval)
  obs_var <- chlaa:::.chlaa_obs_incidence_var(attr(art$fit, "obs_interval") %||% 7)

  set.seed(seed)
  idx <- sample.int(nrow(dr), n, replace = n > nrow(dr))

  P <- matrix(NA_real_, n, length(tvec))
  for (i in seq_len(n)) {
    p <- chlaa:::.chlaa_update_pars_from_theta(dr[idx[i], ], art$pars_warm, art$fit)
    s <- tryCatch(
      chlaa_simulate(p, time = tvec, n_particles = 1, dt = DT,
                     deterministic = FALSE, seed = seed + i),
      error = function(e) NULL)
    if (is.null(s)) next
    mu <- pmax(1e-9, p$reporting_rate * s[[obs_var]])
    P[i, ] <- stats::rnbinom(length(mu), mu = mu, size = p$obs_size)
  }
  P[stats::complete.cases(P), , drop = FALSE]
}

#---- R-hat / ESS recomputed in-script from the raw posterior chains ----
#Recomputed here (NOT read from the stored art$rhat_ess table) so the numbers
#regenerate on any re-run from the fit object alone. Diagnostics are on the
#sampled (transformed) scale.
#
#IMPORTANT: we discard the same warm-up (BURNIN) as the scoring draws before
#computing R-hat / ESS. The fitting script (01_02) diagnoses the FULL chain
#including warm-up, which is (a) internally inconsistent with the posterior we
#actually score (post-burnin) and (b) sensitive to the deliberately dispersed
#R0-based start points, whose transient inflates between-chain variance. We
#subset to the retained iterations via the SAME index the scoring draws use
#(.chlaa_iteration_index), so diagnostics and scoring see identical draws.
#NOTE: this is a consistency fix, not a guaranteed improvement - dropping 25%
#of iterations also lowers ESS, so some zones' numbers get slightly worse.
crps_fit_rhat_ess <- function(fit, burnin = BURNIN) {
  pars_arr <- fit$pars                                          #(param, iter, chain)
  idx <- chlaa:::.chlaa_iteration_index(dim(pars_arr)[2], burnin = burnin, thin = 1)
  pars_arr <- pars_arr[, idx, , drop = FALSE]                   #drop warm-up
  dimnames(pars_arr) <- list(attr(fit, "packer")$names(), NULL, NULL)
  dd <- posterior::as_draws_array(aperm(pars_arr, c(2, 3, 1)))  #(iter, chain, param)
  posterior::summarise_draws(dd, "rhat", "ess_bulk", "ess_tail")
}

#---- Score one zone ----
#scoringRules::crps_sample expects dat with one column per observation and
#returns one CRPS per observation (week); we average over weeks.
crps_score_zone <- function(art) {
  obs <- as.numeric(art$observed$cases)
  P   <- crps_predictive_cases(art)
  if (!nrow(P)) return(NULL)

  #crps_predictive_cases drops failed simulations (complete.cases) silently.
  #Record how many draws actually survived and warn if too few remain, so a
  #zone whose predictive is thinned by simulation failures does not pass
  #unnoticed with an under-powered CRPS estimate.
  n_used <- nrow(P)
  if (n_used < 950) {
    warning(sprintf("%s: only %d/%d predictive simulations succeeded - CRPS may be imprecise",
                    art$hz_name, n_used, N_PRED), call. = FALSE)
  }

  #log scale: transform SAMPLES and OBSERVATIONS, then score
  #(never transform the score itself)
  crps_log <- mean(scoringRules::crps_sample(y = log1p(obs), dat = t(log1p(P))))

  #Prediction-interval coverage: fraction of observed weeks inside the 50% and
  #95% posterior-predictive intervals, computed HERE from the predictive draws
  #P. Reported together as "50% / 95%": adequate 95% with poor 50% coverage
  #indicates the predictive median is systematically displaced.
  pi_cover <- function(w) {
    lo <- apply(P, 2, stats::quantile, probs = 0.5 - w / 2, names = FALSE)
    hi <- apply(P, 2, stats::quantile, probs = 0.5 + w / 2, names = FALSE)
    mean(obs >= lo & obs <= hi)
  }
  cover_50 <- pi_cover(0.50)
  cover_95 <- pi_cover(0.95)

  #R-hat / ESS recomputed in-script from the fit chains, warm-up removed.
  #Min ESS = minimum across bulk AND tail across all four parameters.
  re <- tryCatch(crps_fit_rhat_ess(art$fit), error = function(e) NULL)
  max_rhat <- if (!is.null(re)) max(re$rhat) else NA_real_
  min_ess  <- if (!is.null(re)) min(re$ess_bulk, re$ess_tail) else NA_real_

  #Convergence: three tiers (Vehtari et al. 2021 thresholds). A binary flag is
  #too blunt - it would lump zones that are essentially fine (e.g. R-hat 1.03,
  #ESS ~150) in with genuinely broken chains (kokolo: R-hat 1.39, ESS 7).
  #   Pass     : R-hat <= 1.01 AND min ESS >= 400
  #   Marginal : R-hat <= 1.05 AND min ESS >= 100 (usable with caution)
  #   Fail     : worse than that (do not trust the posterior)
  status <- if (!is.finite(max_rhat) || !is.finite(min_ess)) {
    "Fail"
  } else if (max_rhat <= 1.01 && min_ess >= 400) {
    "Pass"
  } else if (max_rhat <= 1.05 && min_ess >= 100) {
    "Marginal"
  } else {
    "Fail"
  }

  data.frame(
    hz       = art$hz_name,
    n_weeks  = length(obs),
    max_rhat = round(max_rhat, 3),
    min_ess  = round(min_ess, 0),
    status   = status,
    crps_log = round(crps_log, 3),   #<-- HEADLINE, comparable across zones
    cover_50 = round(cover_50, 2),
    cover_95 = round(cover_95, 2)
  )
}

#---- Run (exclude the nyiragongo comparative debugging artifact) ----
cat(rep("=", 70), "\n", sep = "")
cat("CRPS goodness of fit (free-running posterior predictive)\n")
cat(rep("=", 70), "\n\n", sep = "")

crps_files <- fit_files[!grepl("comparative", basename(fit_files))]

crps_res <- do.call(rbind, lapply(crps_files, function(f) {
  art <- readRDS(f)
  cat("scoring:", art$hz_name, "\n")
  crps_score_zone(art)
}))
#Sort by convergence status (Pass -> Marginal -> Fail) then CRPS, so a
#low-CRPS zone that FAILED convergence (e.g. nyiragongo) no longer leads the
#table as if it were the best fit.
crps_res$status <- factor(crps_res$status, levels = c("Pass", "Marginal", "Fail"))
crps_res <- crps_res[order(crps_res$status, crps_res$crps_log), ]

#Export the seven table columns. check.names = FALSE keeps the readable
#headers (spaces, %, the R-hat symbol). The 50%/95% coverage pair is combined
#into one cell so both calibration figures fit without an extra column.
crps_export <- data.frame(
  "Health zone"           = crps_res$hz,
  "Weeks"                 = crps_res$n_weeks,
  "Max R̂"           = crps_res$max_rhat,
  "Min ESS"               = crps_res$min_ess,
  "Status"                = as.character(crps_res$status),
  "CRPS (log1p scale)"    = crps_res$crps_log,
  "PI coverage 50% / 95%" = sprintf("%.2f / %.2f", crps_res$cover_50, crps_res$cover_95),
  check.names = FALSE
)

write.csv(crps_export,
  file.path(tab_dir, "crps_goodness_of_fit.csv"),
  row.names = FALSE)

cat("\n=== CRPS goodness of fit ===\n")
print(crps_export, row.names = FALSE)
cat("\nSorted by convergence status then CRPS (log1p scale); lower CRPS = better fit,\n",
    "comparable across zones. This is CRPS scored on log1p-transformed data, not a\n",
    "log predictive score - no log-score/log-likelihood scoring rule is computed here.\n", sep = "")
cat("Convergence tiers (Vehtari et al. 2021):\n",
    "  Pass     R-hat <= 1.01 AND min ESS >= 400\n",
    "  Marginal R-hat <= 1.05 AND min ESS >= 100\n",
    "  Fail     worse - posterior not trustworthy\n", sep = "")
tier_counts <- table(crps_res$status)
cat(sprintf("Counts: Pass %d, Marginal %d, Fail %d\n",
    tier_counts[["Pass"]], tier_counts[["Marginal"]], tier_counts[["Fail"]]))
fail_zones <- crps_res$hz[crps_res$status == "Fail"]
if (length(fail_zones)) {
  cat("Fail zone(s) (scores not trustworthy):", paste(fail_zones, collapse = ", "), "\n")
}
cat("\n")

#=========================================================================
#VISUALIZATIONS
#=========================================================================

#---- 1. Posterior distributions across HZs (5 params, density curves) ----
#
#ggh4x::facetted_pos_scales is incompatible with ggdist::stat_halfeye
#(produces a blank plot), so we split into two sub-plots:
#   - linear-scale:  R0, frac_neff, obs_size   (3 panels)
#   - log10-scale:   trans_prob, E0             (2 panels)
#Combined with patchwork.

library(ggdist)
library(patchwork)

param_labels <- c(
  R0         = "Basic reproduction number (R0)",
  trans_prob = "Transmission probability",
  frac_neff  = "Effective population fraction",
  obs_size   = "Observation size",
  E0         = "Initial seed (E0)"
)

#Display-friendly health zone names (e.g. "ngiri_ngiri" -> "Ngiri Ngiri",
#"maluku_i" -> "Maluku I"). Capitalises every word, unlike
#tools::toTitleCase() which leaves single-letter words like "i" lowercase.
hz_label <- function(x) {
  x <- gsub("_", " ", x)
  gsub("(?<=^|\\s)([a-z])", "\\U\\1", x, perl = TRUE)
}

#Order HZs by median R0
hz_order <- all_draws |>
  filter(parameter == "R0") |>
  group_by(hz) |>
  summarise(med = median(value), .groups = "drop") |>
  arrange(med) |>
  pull(hz)

draws_all <- all_draws |>
  filter(parameter %in% names(param_labels)) |>
  mutate(hz = factor(hz, levels = hz_order))

#Per-facet summary labels (median [95% CrI])
summ <- draws_all |>
  group_by(hz, parameter) |>
  ggdist::median_qi(value, .width = 0.95) |>
  ungroup()

#Shared theme
theme_halfeye <- theme_bw(base_size = 14, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 13, colour = "black"),
    axis.text        = element_text(colour = "black"),
    axis.text.y      = element_text(size = 13, colour = "black"),
    axis.title       = element_text(colour = "black"),
    plot.title       = element_text(face = "bold", colour = "black")
  )

#--- Standalone R0 plot ---
draws_r0 <- draws_all |> filter(parameter == "R0")
summ_r0  <- summ |> filter(parameter == "R0")

p_r0_density <- ggplot(draws_r0, aes(x = value, y = hz)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey40") +
  stat_halfeye(
    .width = 0.95, fill = "#6baed6",
    normalize = "panels", slab_alpha = 0.8,
    point_size = 1.5, interval_size_range = c(0.6, 1.2)
  ) +
  geom_text(
    data = summ_r0,
    aes(x = 11.2, y = hz,
        label = sprintf("%.2f [%.2f, %.2f]", value, .lower, .upper)),
    hjust = 1, size = 4, colour = "grey30", family = "Helvetica"
  ) +
  scale_y_discrete(labels = hz_label) +
  scale_x_continuous(breaks = seq(0, 10, 2.5)) +
  coord_cartesian(xlim = c(0, 12)) +
  labs(
    x        = expression(R[0]),
    y        = NULL,
    title    = "R0 estimates by health zone",
    subtitle = "Densities with median and 95% uncertainty interval"
  ) +
  theme_halfeye +
  theme(
    panel.grid  = element_blank(),
    axis.text.y = element_text(size = 4 * .pt)
  )

ggsave(file.path(fig_dir, "fitted_r0_density.png"),
  p_r0_density, width = 8, height = 7, dpi = 600)

#--- 4-parameter comparison (2x2): linear pair + log pair via patchwork ---
#
#frac_neff and obs_size are built as SEPARATE single-parameter plots (rather
#than one facet_wrap, as before) so obs_size can get its own x-axis limit -
#ggh4x::facetted_pos_scales (the usual way to give one facet an independent
#scale) is incompatible with stat_halfeye (see note above), so per-facet
#scale control here means per-plot instead.
label_above <- function(summ_dat) {
  geom_label(
    data = summ_dat,
    aes(x = value, y = hz,
        label = sprintf("%.3g [%.3g, %.3g]", value, .lower, .upper)),
    nudge_y = 0.32, vjust = 0, size = 2.8, colour = "black", family = "Helvetica",
    fill = "white", label.size = 0, label.padding = unit(0.08, "lines")
  )
}

lin_params  <- c("frac_neff", "obs_size")
draws_lin   <- draws_all |>
  filter(parameter %in% lin_params) |>
  mutate(parameter = factor(parameter, levels = lin_params))
summ_lin    <- summ |>
  filter(parameter %in% lin_params) |>
  mutate(parameter = factor(parameter, levels = lin_params))

p_frac_neff <- ggplot(draws_lin |> filter(parameter == "frac_neff"), aes(x = value, y = hz)) +
  stat_halfeye(
    .width = 0.95, fill = "#6baed6",
    normalize = "panels", slab_alpha = 0.8,
    point_size = 1.5, interval_size = 0.6
  ) +
  label_above(summ_lin |> filter(parameter == "frac_neff")) +
  facet_wrap(~ parameter, scales = "free_x",
             labeller = labeller(parameter = param_labels)) +
  scale_y_discrete(labels = hz_label) +
  #Extra horizontal room so the value labels (centred over each median) are
  #not clipped when a median sits near a panel edge.
  scale_x_continuous(expand = expansion(mult = c(0.18, 0.08))) +
  labs(x = NULL, y = NULL) +
  theme_halfeye +
  theme(panel.grid = element_blank())

p_obs_size <- ggplot(draws_lin |> filter(parameter == "obs_size"), aes(x = value, y = hz)) +
  stat_halfeye(
    .width = 0.95, fill = "#6baed6",
    normalize = "panels", slab_alpha = 0.8,
    point_size = 1.5, interval_size = 0.6
  ) +
  label_above(summ_lin |> filter(parameter == "obs_size")) +
  facet_wrap(~ parameter, scales = "free_x",
             labeller = labeller(parameter = param_labels)) +
  coord_cartesian(xlim = c(0, 15)) +
  scale_y_discrete(labels = hz_label) +
  labs(x = NULL, y = NULL) +
  theme_halfeye +
  #Second HZ y-axis is redundant with the left panel (same zone order) - hide it
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

log_params  <- c("trans_prob", "E0")
draws_log   <- draws_all |>
  filter(parameter %in% log_params) |>
  mutate(parameter = factor(parameter, levels = log_params))
summ_log    <- summ |>
  filter(parameter %in% log_params) |>
  mutate(parameter = factor(parameter, levels = log_params))

p_log <- ggplot(draws_log, aes(x = value, y = hz)) +
  stat_halfeye(
    .width = 0.95, fill = "#6baed6",
    normalize = "panels", slab_alpha = 0.8,
    point_size = 1.5, interval_size = 0.6
  ) +
  label_above(summ_log) +
  facet_wrap(~ parameter, scales = "free_x", nrow = 1,
             labeller = labeller(parameter = param_labels)) +
  scale_x_log10(expand = expansion(mult = c(0.12, 0.13))) +
  scale_y_discrete(labels = hz_label) +
  labs(x = NULL, y = NULL) +
  theme_halfeye +
  theme(panel.grid = element_blank())

#Stack linear pair (top) and log pair (bottom) into a 2x2 grid
p_dist <- ((p_frac_neff | p_obs_size) / p_log) +
  plot_annotation(
    title    = "Posterior parameter estimates by health zone",
    subtitle = "Half-eye densities with median and 95% uncertainty intervals",
    theme    = theme(
      text          = element_text(family = "Helvetica"),
      plot.title    = element_text(face = "bold", size = 16, colour = "black"),
      plot.subtitle = element_text(size = 12.5, colour = "black")
    )
  )

ggsave(file.path(fig_dir, "fitted_parameters_comparison.png"),
  p_dist, width = 15, height = 13, dpi = 600)

#---- 2. R0 estimates across HZs ----

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
    theme_bw(base_size = 12, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(fig_dir, "fitted_r0_comparison.png"),
    p_r0, width = 10, height = 8, dpi = 300)
}

#---- 3. trans_prob vs frac_neff (the identifiability ridge) ----

param_wide_plot <- all_params |>
  select(hz, parameter, median, total_cases) |>
  pivot_wider(names_from = parameter, values_from = median)

if (all(c("trans_prob", "frac_neff") %in% colnames(param_wide_plot))) {
  p_ridge <- ggplot(param_wide_plot, aes(x = trans_prob, y = frac_neff)) +
    geom_point(aes(size = total_cases, colour = total_cases)) +
    geom_text(aes(label = hz), hjust = -0.1, vjust = -0.1, size = 3, family = "Helvetica") +
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
    theme_bw(base_size = 12, base_family = "Helvetica")

  ggsave(file.path(fig_dir, "fitted_params_trans_vs_frac_neff.png"),
    p_ridge, width = 10, height = 8, dpi = 300)
}

#---- 4. obs_size vs total_cases ----

if ("obs_size" %in% colnames(param_wide_plot)) {
  p_obs <- ggplot(param_wide_plot, aes(x = total_cases, y = obs_size)) +
    geom_point(aes(colour = frac_neff), size = 4) +
    geom_text(aes(label = hz), hjust = -0.1, vjust = -0.1, size = 3, family = "Helvetica") +
    scale_x_log10() +
    scale_colour_viridis_c() +
    labs(
      title  = "Observation Overdispersion vs Outbreak Size",
      x = "Total Cases (log scale)",
      y = "Observation Size (obs_size)",
      colour = "frac_neff"
    ) +
    theme_bw(base_size = 12, base_family = "Helvetica")

  ggsave(file.path(fig_dir, "fitted_params_obssize_vs_cases.png"),
    p_obs, width = 10, height = 8, dpi = 300)
}

#---- 5. Acceptance rates ----

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
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "fitted_acceptance_rates.png"),
  p_accept, width = 10, height = 6, dpi = 300)

#---- 6. Convergence diagnostics (R-hat) ----

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
    theme_bw(base_size = 12, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )

  ggsave(file.path(fig_dir, "fitted_rhat_diagnostics.png"),
    p_rhat, width = 12, height = 10, dpi = 300)
}

#---- 6b. Particle-count / log-likelihood variance check ----

if (all(c("n_particles_used", "loglik_var_at_n_prod") %in% names(diag_table))) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("Particle-count / log-likelihood variance check\n")
  cat(rep("=", 70), "\n\n", sep = "")
  print(diag_table |>
    select(hz, n_particles_used, loglik_var_at_n_prod, loglik_var_target),
  n = Inf)
}

#---- 7. CRPS goodness of fit by health zone ----
#
#Lollipop of crps_log (comparable across zones; lower = better fit), best at
#top. Points are coloured by the three-tier convergence status (Pass / Marginal
#/ Fail) so zones that are essentially fine are not lumped in with genuinely
#broken chains (kokolo). The 3-seed precision statement in the caption sets the
#smallest CRPS gap worth interpreting.

status_cols <- c("Pass" = "#2f6a4e", "Marginal" = "#e08214", "Fail" = "#c0392b")
status_shp  <- c("Pass" = 16,        "Marginal" = 17,         "Fail" = 1)

crps_plot_df <- crps_res |>
  mutate(hz = factor(hz, levels = crps_res$hz[order(-crps_res$crps_log)]))

p_crps <- ggplot(crps_plot_df, aes(x = crps_log, y = hz)) +
  geom_segment(aes(x = 0, xend = crps_log, yend = hz),
    colour = "grey70", linewidth = 0.5) +
  geom_point(aes(colour = status, shape = status), size = 3) +
  geom_text(aes(label = sprintf("%.2f", crps_log)),
    hjust = -0.35, size = 3, colour = "black", family = "Helvetica") +
  scale_colour_manual(values = status_cols, name = NULL, drop = FALSE) +
  scale_shape_manual(values = status_shp, name = NULL, drop = FALSE) +
  scale_y_discrete(labels = hz_label) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.16))) +
  labs(
    x        = expression(CRPS[log] ~ "(lower = better fit)"),
    y        = NULL,
    title    = "CRPS goodness of fit by health zone",
    subtitle = "Free-running posterior predictive, log1p scale (comparable across zones)",
    caption  = paste0(
      "Convergence (Vehtari et al. 2021): Pass R-hat<=1.01 & ESS>=400; ",
      "Marginal <=1.05 & >=100; Fail worse.\n",
      "CRPS from 1000 draws; 3-seed spread <0.01, so gaps below ~0.03 are not interpreted.")
  ) +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text          = element_text(colour = "black"),
    axis.title         = element_text(colour = "black"),
    plot.title         = element_text(face = "bold", colour = "black"),
    plot.subtitle      = element_text(colour = "grey30"),
    plot.caption       = element_text(colour = "grey30", hjust = 0, size = 9),
    legend.position    = "bottom"
  )

ggsave(file.path(fig_dir, "fitted_crps_goodness_of_fit.png"),
  p_crps, width = 9.5, height = max(2.4, 0.45 * nrow(crps_res) + 1.2), dpi = 300)

#---- Summary ----

cat("\nSaved figures:\n")
cat("  ", file.path(fig_dir, "fitted_crps_goodness_of_fit.png"), "\n")
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

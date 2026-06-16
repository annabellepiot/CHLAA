# =========================================================================
# kirotshe.R — Diagnostic fitting for Kirotshe health zone
# =========================================================================
#
# Uses chlaa_fit_pmcmc() with custom 4-parameter prior/packer.
# Fixed: contact_rate=0, incubation_time=4.845, duration_sym=14.48,
#        fatality_untreated=0.0043, fatality_treated=0.001
# Fitted: trans_prob, reporting_rate, obs_size, seek_severe
#
# Odin model uses zero_every = 7 (weekly accumulator), so IDSR data
# is passed directly as weekly observations at 7-day intervals.
#
# Two-stage approach:
#   1. Exploratory fit (100 particles, 1000 steps)
#   2. Production fit (200 particles, 10000 steps, warm-started from
#      exploratory posterior median)
#
# =========================================================================

library(chlaa)
library(ggplot2)
library(tidyverse)

data_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/CHLAA/chlaa/analysis/data"
output_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/CHLAA/output"
fig_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/CHLAA/figures"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load data ----

hz_params_long <- read.csv(file.path(data_dir, "hz_parameters.csv"),
    stringsAsFactors = FALSE
)
hz_params <- hz_params_long %>%
    pivot_wider(names_from = parameter, values_from = value)

idsr <- read.csv(file.path(data_dir, "IDSR_dataset.csv"))

# ---- 2. Extract Kirotshe parameters ----

i <- which(hz_params$hz == "kirotshe")
outbreak_start <- as.Date(hz_params$outbreak_start[i])
outbreak_end <- as.Date(hz_params$outbreak_end[i])

# Convert intervention dates to integer day offsets from outbreak_start.
# Returns 0L if the date is NA (intervention inactive).
safe_date_to_day <- function(date_str, origin) {
    d <- as.Date(date_str, format = "%Y-%m-%d")
    if (is.na(d)) {
        return(0L)
    }
    as.integer(d - origin)
}

orc_start_day <- safe_date_to_day(hz_params$orc_start[i], outbreak_start)
orc_end_day <- safe_date_to_day(hz_params$orc_end[i], outbreak_start)
ctc_start_day <- safe_date_to_day(hz_params$ctc_start[i], outbreak_start)
ctc_end_day <- safe_date_to_day(hz_params$ctc_end[i], outbreak_start)
chlor_start_day <- safe_date_to_day(hz_params$chlor_start[i], outbreak_start)
chlor_end_day <- safe_date_to_day(hz_params$chlor_end[i], outbreak_start)
hyg_start_day <- safe_date_to_day(hz_params$hyg_start[i], outbreak_start)
hyg_end_day <- safe_date_to_day(hz_params$hyg_end[i], outbreak_start)
cati_start_day <- safe_date_to_day(hz_params$cati_start[i], outbreak_start)
cati_end_day <- safe_date_to_day(hz_params$cati_end[i], outbreak_start)
lat_start_day <- safe_date_to_day(hz_params$lat_start[i], outbreak_start)
lat_end_day <- safe_date_to_day(hz_params$lat_end[i], outbreak_start)

chlor_effect_val <- as.numeric(hz_params$chlor_effect[i])
hyg_effect_val <- as.numeric(hz_params$hyg_effect[i])
cati_effect_val <- as.numeric(hz_params$cati_effect[i])
lat_effect_val <- as.numeric(hz_params$lat_effect[i])

cat("Outbreak window:", as.character(outbreak_start), "to", as.character(outbreak_end), "\n")

# ---- 3. Prepare weekly IDSR data ----

hz_weekly <- idsr %>%
    filter(hz == "kirotshe") %>%
    mutate(date = as.Date(date)) %>%
    select(date, year, week, cases, deaths, population)

hz_outbreak <- hz_weekly %>%
    filter(date >= outbreak_start, date <= outbreak_end) %>%
    arrange(date)

cat("Outbreak weeks:", nrow(hz_outbreak), "\n")
cat("Total weekly cases:", sum(hz_outbreak$cases), "\n")

# Weekly time points at 7-day intervals (matching zero_every = 7 accumulator)
hz_data_weekly <- hz_outbreak %>%
    mutate(time = seq_len(n()) * 7L) %>%
    select(time, date, cases) %>%
    arrange(time)

# ---- 4. Initial seeding from data ----

pop_hz <- hz_weekly$population[1]
seed_date <- outbreak_start - 14
seed_row <- hz_weekly %>%
    filter(date <= seed_date) %>%
    arrange(desc(date)) %>%
    slice(1)

# Scale observed cases by reporting rate to estimate true unobserved burden
expected_reporting_rate <- 0.3520
if (nrow(seed_row) > 0 && seed_row$cases > 0) {
    E0_val <- ceiling(seed_row$cases / expected_reporting_rate)
} else {
    E0_val <- ceiling(max(1, hz_outbreak$cases[1]) / expected_reporting_rate)
}

cat(sprintf("Initial seeding: E0=%d\n", E0_val))

# ---- 5. Set starting parameters ----

pars_args <- list(
    N = pop_hz,
    Sev0 = 0,
    E0 = E0_val,
    M0 = 0,
    immunity_asym = 280,
    contact_rate = 0, # setting to 0 for now as simplification, although chlaa is set to allow fitting this
    trans_prob = 0.003225,
    incubation_time = 4.845,
    duration_sym = 14.48,
    seek_mild = 0.1,
    seek_severe = 0.4086,
    vax2_doses_per_day = 0,
    vax2_total_doses = 0,
    reporting_rate = 0.3520,
    fatality_treated = 0.001,
    fatality_untreated = 0.0043,
    obs_size = 30,
    orc_start = orc_start_day,
    orc_end = orc_end_day,
    ctc_start = ctc_start_day,
    ctc_end = ctc_end_day,
    hyg_start = hyg_start_day,
    hyg_end = hyg_end_day,
    hyg_effect = hyg_effect_val,
    cati_start = cati_start_day,
    cati_end = cati_end_day,
    cati_effect = cati_effect_val
)

if (chlor_start_day > 0) {
    pars_args$chlor_start <- chlor_start_day
    pars_args$chlor_end <- chlor_end_day
    pars_args$chlor_effect <- chlor_effect_val
}
if (lat_start_day > 0) {
    pars_args$lat_start <- lat_start_day
    pars_args$lat_end <- lat_end_day
    pars_args$lat_effect <- lat_effect_val
}

pars <- do.call(chlaa_parameters, pars_args)

# ---- 6. Prepare weekly data for the fitter ----

fit_data <- data.frame(time = hz_data_weekly$time, cases = hz_data_weekly$cases)

# ---- 7. Custom 4-parameter prior and packer ----
# Fit only: trans_prob, reporting_rate, obs_size, seek_severe

custom_prior <- monty::monty_dsl({
    trans_prob ~ Gamma(shape = 2, rate = 200) # mean = 0.01, not flat
    reporting_rate ~ Beta(2, 20) # peaks around 0.1, should work for our calibration strating value found ~0.04.
    obs_size ~ Gamma(shape = 5, rate = 0.1) # shaped mean = 20, enough spread for weekly data
    seek_severe ~ Beta(2, 1.5) # mean = 0.57, a bit flat but still favors higher values (consistent with our calibration starting value of 1)
})

custom_packer <- function(pars) {
    names_fit <- c("trans_prob", "reporting_rate", "obs_size", "seek_severe")
    fixed <- pars[setdiff(names(pars), names_fit)]
    monty::monty_packer(names_fit, fixed = fixed)
}

# ---- 8. Exploratory fit: 100 particles, 1000 steps ----
# Full proposal covariance matrix from deterministic calibration (Hessian-derived, Gelman-scaled).
# BUT We must use a purely diagonal matrix for the exploratory run.
# The deterministic Hessian provided mathematically impossible correlations,
# so we discard just the off-diagonals and let the exploratory run learn them empirically.
# Note: Order must match custom_packer: trans_prob, reporting_rate, obs_size, seek_severe.
# The deterministic fit calibrated trans_prob, reporting_rate, seek_severe (indices 1,2,4).
# obs_size (index 3) was not in the deterministic fit, so it gets (0.1 * start)^2 on the diagonal.

explore_proposal <- matrix(0, 4, 4)

# Set the step sizes (variances) on the diagonal
explore_proposal[1, 1] <- 1.102081e-08 # trans_prob
explore_proposal[2, 2] <- 1.651216e-04 # reporting_rate
explore_proposal[3, 3] <- (0.1 * 30)^2 # obs_size (guess)
explore_proposal[4, 4] <- 1.669885e-03 # seek_severe

# Safeguard: Ensure no variance collapsed to exactly 0 or negative
mle_starts <- c(0.003225, 0.3520, 30, 0.4086)
for (k in 1:4) {
    if (explore_proposal[k, k] <= 0) {
        explore_proposal[k, k] <- (0.1 * mle_starts[k])^2
    }
}

# (No eigen decomposition needed; a matrix with positive values only on
# the diagonal is guaranteed to be positive definite!)

cat("\n=== EXPLORATORY FIT ===\n")
fit_explore <- chlaa_fit_pmcmc(
    data = fit_data,
    pars = pars,
    n_particles = 100,
    n_steps = 1000,
    seed = 42,
    prior = custom_prior,
    packer = custom_packer(pars),
    proposal_var = explore_proposal
)

# ---- 9. Exploratory diagnostics ----

report_explore <- chlaa_fit_report(fit_explore, burnin = 0.25, thin = 2)
cat("Exploratory acceptance rate:", report_explore$acceptance_rate, "\n")
print(report_explore$posterior_summary)

# ---- 10. Warm start from exploratory posterior ----

packer <- attr(fit_explore, "packer")
d <- length(packer$names())
n_samples <- dim(fit_explore$pars)[2]
start_idx <- floor(0.25 * n_samples) + 1

# Drop the chain dimension if present (3D array from monty with 1 chain)
pars_mat <- if (length(dim(fit_explore$pars)) == 3) {
    fit_explore$pars[, start_idx:n_samples, 1]
} else {
    fit_explore$pars[, start_idx:n_samples]
}
pooled <- t(pars_mat)
colnames(pooled) <- packer$names()

warm_vec <- apply(pooled, 2, median)
pars_warm <- pars
for (nm in names(warm_vec)) pars_warm[[nm]] <- warm_vec[[nm]]

# Compute full VCV from exploratory posterior and apply Gelman scaling
d <- length(packer$names()) # d = 4 parameters
prod_proposal <- cov(pooled) * (2.38^2 / d)

eig <- eigen(prod_proposal, symmetric = TRUE)
eig$values <- pmax(eig$values, 1e-12)
prod_proposal <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)

cat("\nWarm-start parameter values (posterior median from exploratory run):\n")
for (nm in packer$names()) {
    cat(sprintf("  %s = %.6f\n", nm, pars_warm[[nm]]))
}
cat("Production proposal covariance matrix:\n")
print(prod_proposal)

rm(fit_explore, report_explore, packer, pooled, pars_mat)

# ---- 11. Production fit: 200 particles, 10000 steps, warm start ----

cat("\n=== PRODUCTION FIT ===\n")
fit <- chlaa_fit_pmcmc(
    data = fit_data,
    pars = pars_warm,
    n_particles = 200,
    n_steps = 10000,
    seed = 123,
    prior = custom_prior,
    packer = custom_packer(pars_warm),
    proposal_var = prod_proposal
)

# ---- 12. Production diagnostics ----

packer_prod <- attr(fit, "packer")
fitted_names <- packer_prod$names()

report_prod <- chlaa_fit_report(fit, burnin = 0.25, thin = 2)
cat("Production acceptance rate:", report_prod$acceptance_rate, "\n")
print(report_prod$posterior_summary)

p_trace <- chlaa_plot_trace(
    fit,
    parameters = fitted_names,
    burnin = 0, thin = 1
)
ggsave(file.path(fig_dir, "diagnosis_kirotshe_production_trace.png"),
    p_trace,
    width = 14, height = 10, dpi = 150
)

p_lltrace <- chlaa_plot_likelihood_trace(fit, burnin = 0.25, thin = 2)
ggsave(file.path(fig_dir, "diagnosis_kirotshe_production_likelihood_trace.png"),
    p_lltrace,
    width = 8, height = 5, dpi = 150
)

p_pairs <- chlaa_plot_parameter_pairs(
    fit,
    parameters = fitted_names,
    burnin = 0.25, thin = 2
)
ggsave(file.path(fig_dir, "diagnosis_kirotshe_production_pairs.png"),
    p_pairs,
    width = 12, height = 12, dpi = 150
)

# ---- 13. Posterior predictive projection ----
# Forecast at weekly intervals (zero_every = 7), no aggregation needed.

n_weeks <- nrow(hz_data_weekly)
forecast_time <- seq(7, n_weeks * 7, by = 7)

fc <- chlaa_forecast_from_fit(
    fit = fit,
    pars = pars_warm,
    time = forecast_time,
    vars = "inc_symptoms",
    include_cases = TRUE,
    quantiles = c(0.025, 0.125, 0.25, 0.5, 0.75, 0.875, 0.975),
    n_draws = 100,
    burnin = 0.25,
    thin = 2,
    seed = 123,
    n_particles = 1
)

f_cases <- fc %>%
    filter(variable == "cases") %>%
    mutate(date = outbreak_start + time)

# ---- 14. Fit plot ----
# Plot weekly observed data against weekly forecast

p_fit <- ggplot() +
    geom_col(
        data = hz_data_weekly, aes(x = date, y = cases),
        fill = "grey70", width = 6, alpha = 0.6
    ) +
    geom_ribbon(
        data = f_cases,
        aes(x = date, ymin = q0p025, ymax = q0p975, fill = "95% CI"),
        alpha = 0.10
    ) +
    geom_ribbon(
        data = f_cases,
        aes(x = date, ymin = q0p25, ymax = q0p75, fill = "50% CI"),
        alpha = 0.25
    ) +
    geom_line(
        data = f_cases, aes(x = date, y = mean, colour = "Mean fit"),
        linewidth = 0.8
    ) +
    geom_point(
        data = hz_data_weekly, aes(x = date, y = cases, colour = "IDSR data"),
        size = 1.5
    ) +
    scale_colour_manual(values = c("IDSR data" = "black", "Mean fit" = "#2c7f62")) +
    scale_fill_manual(values = c("50% CI" = "#b88a66", "95% CI" = "#9a3b32")) +
    labs(
        title = "Kirotshe — chlaa_fit_pmcmc (10000 steps, 200 particles, warm start, 4-par: trans_prob, reporting_rate, obs_size, seek_severe)",
        x = NULL, y = "Cases/week",
        colour = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 12)

ggsave(file.path(fig_dir, "diagnosis_kirotshe_production_fit.png"),
    p_fit,
    width = 10, height = 6, dpi = 150
)

cat("\nDone. Fit plot saved to:", file.path(fig_dir, "diagnosis_kirotshe_production_fit.png"), "\n")

# =========================================================================
# All_HZ_fitting.r — Production fitting for all health zones
# =========================================================================
#
# The fitting workflow follows the vignette at:
#   https://ojwatson.github.io/chlaa/articles/fitting_workflow.html
# You will need:
# 1) A CSV file dataset equivalent to our real IDSR data (idsr_clean.csv)
#    (not made available in this repo) with the following columns:
#   - hz_name, year, epi week number, date (YYYY-MM-DD),
#     cases, deaths, population size
# 2) A CSV file (hz_parameters.csv) with intervention dates and effects.
# _______NOTE:___must add Nyiragongo May-Dec 2025,
# and choose outbreak dates and add for Kalemie, Nyemba, Kabalo, Kongolo____###
#
# 3) We also apply the parameter decisions (A05-A22) listed in
#    Model Assumptions.qmd. Verify the assumptions, and set custom
#    parameters below as needed.
#
# This script loops over all target health zones, reading their
# intervention dates and outbreak windows from hz_parameters.csv.
# For each HZ it produces:
#   1. The production fit curve (fitted model vs IDSR data)
#   2. Diagnostic plots (trace, likelihood trace, parameter pairs)
#   3. Saved fitted parameters ({hz_name}_fitted_pars.rds) for use by
#      Scenario_Analysis_Interventions_All_HZ.R
#
# =========================================================================
# KEY METHODOLOGY CHOICES
# =========================================================================
#
# A) DAILY FITTING (zero_every = 1 in odin model):
#    The odin accumulator inc_symptoms resets every day.
#    Weekly IDSR cases are disaggregated to daily pseudo-observations
#    (cases / 7, rounded with remainder distributed across first days).
#
# B) PURE ENVIRONMENTAL (WATERBORNE) TRANSMISSION:
#    contact_rate is fixed at 0, removing person-to-person transmission.
#    trans_prob is kept as a fitted parameter because it also scales the
#    environmental force of infection:
#      env_force = trans_prob * C / (C + contam_half_sat)
#
# C) DATA-DRIVEN INITIAL SEEDING:
#    outbreak_start in hz_parameters.csv is the first nonzero case week
#    within 3 months before the first intervention date.
#    E0 is set to cases at outbreak_start (the first nonzero week).
#
# D) FOUR FITTED PARAMETERS:
#    trans_prob, reporting_rate, obs_size, seek_severe.
#    All others fixed at literature/regional values.
#
# E) PER-PARAMETER PROPOSAL VARIANCES:
#    chlaa_fit_pmcmc accepts a vector of per-parameter proposal
#    variances (fit.R was modified to support this).
#    The exploratory run uses hand-tuned diagonal proposals.
#    The production run uses posterior variance from the exploratory fit,
#    floored by the exploratory proposals.
#
# F) TWO-STAGE WARM START:
#    Exploratory (10000 steps, 60 particles) → extract posterior median →
#    Production (20000 steps, 200 particles) with learned proposals.
#
# G) USES chlaa_fit_pmcmc():
#    The package's built-in fitting function handles filter creation,
#    sampling, and chlaa_fit object tagging internally.
# =========================================================================

library(chlaa)
library(ggplot2)
library(tidyverse)

data_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/CHLAA/chlaa/analysis/data"
output_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/CHLAA/output"
fig_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/CHLAA/figures"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load data ----

hz_params_long <- read.csv(file.path(data_dir, "hz_parameters.csv"),
    stringsAsFactors = FALSE
)
hz_params <- hz_params_long %>%
    pivot_wider(names_from = parameter, values_from = value)

idsr <- read.csv(file.path(data_dir, "IDSR_dataset.csv"))

# ---- Helpers ----

safe_date <- function(x) as.Date(x, format = "%Y-%m-%d")

# Convert a date string to integer day offset from outbreak_start.
# Returns 0L if the date is NA (intervention inactive).
date_to_day <- function(date_str, outbreak_start) {
    d <- safe_date(date_str)
    if (is.na(d)) return(0L)
    as.integer(d - outbreak_start)
}

# =========================================================================
# FITTING CONFIGURATION
# =========================================================================

# ---- Custom prior: 4 fitted parameters ----
# Gamma(5, 0.25) on obs_size prevents the degenerate mode where the
# sampler exploits maximum overdispersion (obs_size → 1) instead of
# fitting the epidemic curve. Mean = 20, mode = 16, P(< 3) ≈ 0.2%.
custom_prior <- monty::monty_dsl({
    trans_prob ~ Uniform(0.001, 0.2)
    reporting_rate ~ Uniform(0.01, 1.0)
    obs_size ~ Gamma(shape = 5, rate = 0.25)
    seek_severe ~ Uniform(0.1, 1.0)
})

# ---- Custom packer (4 parameters, matching the prior) ----
custom_packer <- function(pars) {
    names_fit <- c("trans_prob", "reporting_rate", "obs_size", "seek_severe")
    fixed <- pars[setdiff(names(pars), names_fit)]
    monty::monty_packer(names_fit, fixed = fixed)
}

# ---- Per-parameter exploratory proposal variances ----
# Default proposals (used if no pre-tuned values available).
# Tuned to parameter scales; trans_prob is very sensitive (~0.001)
# and needs tiny steps (variance 1e-8, sd = 1e-4).
default_proposal <- c(
    trans_prob = 1e-8,
    reporting_rate = 1e-4,
    obs_size = 1.0,
    seek_severe = 0.01
)

fitted_names <- c("trans_prob", "reporting_rate", "obs_size", "seek_severe")

# ---- Load pre-tuned proposals if available ----
# Generated by pretune_proposals.R. Falls back to defaults if missing.
proposals_path <- file.path(data_dir, "Copy_hz_proposals.csv")
if (file.exists(proposals_path)) {
    hz_proposals <- read.csv(proposals_path, stringsAsFactors = FALSE)
    cat("Loaded pre-tuned proposals from:", proposals_path, "\n")
    cat("HZs with proposals:", paste(hz_proposals$hz, collapse = ", "), "\n\n")
} else {
    hz_proposals <- NULL
    cat("No pre-tuned proposals found — using defaults for all HZs.\n\n")
}

# =========================================================================
# MAIN LOOP — iterate over each health zone
# =========================================================================

for (i in seq_len(nrow(hz_params))) {
    hz_name <- hz_params$hz[i]
    cat("\n=========================================================================\n")
    cat("Processing health zone:", hz_name, "\n")
    cat("=========================================================================\n\n")

    # ---- 0. Look up per-HZ best starting values from pretune ----
    start_trans_prob      <- 0.001
    start_reporting_rate  <- 0.20
    start_obs_size        <- 30
    start_seek_severe     <- 0.5
    if (!is.null(hz_proposals) && hz_name %in% hz_proposals$hz) {
        row <- hz_proposals[hz_proposals$hz == hz_name, ]
        if (!is.na(row$best_trans_prob))      start_trans_prob      <- row$best_trans_prob
        if (!is.na(row$best_reporting_rate))  start_reporting_rate  <- row$best_reporting_rate
        if (!is.na(row$best_obs_size))        start_obs_size        <- row$best_obs_size
        if (!is.na(row$best_seek_severe))     start_seek_severe     <- row$best_seek_severe
        cat("Pre-tuned starting values: trans_prob =", start_trans_prob,
            ", reporting_rate =", start_reporting_rate,
            ", obs_size =", start_obs_size,
            ", seek_severe =", start_seek_severe, "\n")
    } else {
        cat("No pre-tune available — using default starting values\n")
    }
    explore_proposal <- unname(default_proposal)

    # ---- 1. Extract parameters from the CSV row ----

    outbreak_start <- as.Date(hz_params$outbreak_start[i])
    outbreak_end <- as.Date(hz_params$outbreak_end[i])

    orc_start_day <- date_to_day(hz_params$orc_start[i], outbreak_start)
    orc_end_day <- date_to_day(hz_params$orc_end[i], outbreak_start)
    ctc_start_day <- date_to_day(hz_params$ctc_start[i], outbreak_start)
    ctc_end_day <- date_to_day(hz_params$ctc_end[i], outbreak_start)
    chlor_start_day <- date_to_day(hz_params$chlor_start[i], outbreak_start)
    chlor_end_day <- date_to_day(hz_params$chlor_end[i], outbreak_start)
    hyg_start_day <- date_to_day(hz_params$hyg_start[i], outbreak_start)
    hyg_end_day <- date_to_day(hz_params$hyg_end[i], outbreak_start)
    cati_start_day <- date_to_day(hz_params$cati_start[i], outbreak_start)
    cati_end_day <- date_to_day(hz_params$cati_end[i], outbreak_start)
    lat_start_day <- date_to_day(hz_params$lat_start[i], outbreak_start)
    lat_end_day <- date_to_day(hz_params$lat_end[i], outbreak_start)

    chlor_effect_val <- as.numeric(hz_params$chlor_effect[i])
    hyg_effect_val <- as.numeric(hz_params$hyg_effect[i])
    cati_effect_val <- as.numeric(hz_params$cati_effect[i])
    lat_effect_val <- as.numeric(hz_params$lat_effect[i])

    # ---- 2. Load and prepare IDSR data ----

    hz_weekly <- idsr %>%
        filter(hz == hz_name) %>%
        mutate(date = as.Date(date)) %>%
        select(date, year, week, cases, deaths, population)

    hz_outbreak <- hz_weekly %>%
        filter(date >= outbreak_start, date <= outbreak_end) %>%
        arrange(date)

    cat("Outbreak window:", as.character(outbreak_start), "to", as.character(outbreak_end), "\n")
    cat("Outbreak weeks:", nrow(hz_outbreak), "\n")
    cat("Total weekly cases:", sum(hz_outbreak$cases), "\n")

    if (nrow(hz_outbreak) == 0) {
        cat("WARNING: No data found for", hz_name, "— skipping.\n")
        next
    }

    # ---- 3. Disaggregate weekly cases to daily pseudo-observations ----
    # With zero_every = 1 in the odin model, inc_symptoms resets daily.
    # Each weekly count is spread evenly across 7 days (rounded to integers,
    # remainder assigned to the first days of the week).

    hz_daily <- hz_outbreak %>%
        rowwise() %>%
        mutate(daily_list = list({
            base <- cases %/% 7L
            remainder <- cases %% 7L
            daily_cases <- rep(base, 7)
            if (remainder > 0) daily_cases[1:remainder] <- daily_cases[1:remainder] + 1L
            data.frame(
                date = date + 0:6,
                cases = as.integer(daily_cases)
            )
        })) %>%
        ungroup() %>%
        select(daily_list) %>%
        tidyr::unnest(daily_list) %>%
        mutate(time = as.integer(date - min(date))) %>%
        select(time, date, cases) %>%
        arrange(time)

    cat("Daily observation rows:", nrow(hz_daily), "\n")
    cat("Total daily cases:", sum(hz_daily$cases), "\n")

    # Keep weekly data for the fit plot
    hz_data_weekly <- hz_outbreak %>%
        mutate(time = as.integer(date - min(date))) %>%
        select(time, date, cases) %>%
        arrange(time)

    # ---- 4. Initial seeding from data ----
    # E0 = cases at outbreak_start (the first nonzero case week).
    # outbreak_start in hz_parameters.csv is already set to the first
    # nonzero case week within 3 months of the first intervention date.

    pop_hz <- hz_weekly$population[1]
    E0_val <- max(1, hz_outbreak$cases[1])
    cat(sprintf("Initial seeding: E0=%d (from outbreak_start week)\n", E0_val))

    # ---- 5. Set starting parameters ----
    # contact_rate = 0: pure environmental transmission
    # trans_prob, reporting_rate, obs_size, seek_severe: from pre-tuned
    #   best starting values (grid search in pretune), falls back to defaults.
    # Fixed at literature/regional values:
    #   incubation_time = 4.845 days, duration_sym = 14.48 days,
    #   fatality_treated = 0.001 (North-Kivu), fatality_untreated = 0.0043

    pars_args <- list(
        N = pop_hz,
        Sev0 = 0,
        E0 = E0_val,
        M0 = 0,
        immunity_asym = 280,
        contact_rate = 0,
        trans_prob = start_trans_prob,
        incubation_time = 4.845,
        duration_sym = 14.48,
        seek_mild = 0.1,
        seek_severe = start_seek_severe,
        vax2_doses_per_day = 0,
        vax2_total_doses = 0,
        reporting_rate = start_reporting_rate,
        fatality_treated = 0.001,
        fatality_untreated = 0.0043,
        obs_size = start_obs_size,
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

    # ---- 6. Prepare daily data for the fitter ----

    fit_data <- data.frame(time = hz_daily$time, cases = hz_daily$cases)

    # ---- 7. Exploratory fit: 60 particles, 10000 steps ----

    cat("\n=== EXPLORATORY FIT ===\n")
    fit_explore <- tryCatch(
        chlaa_fit_pmcmc(
            data = fit_data,
            pars = pars,
            n_particles = 60,
            n_steps = 10000,
            seed = 42,
            prior = custom_prior,
            packer = custom_packer(pars),
            proposal_var = explore_proposal
        ),
        error = function(e) {
            cat("EXPLORATORY FIT FAILED for", hz_name, ":", conditionMessage(e), "\n")
            NULL
        }
    )

    if (is.null(fit_explore)) {
        cat("Skipping", hz_name, "due to exploratory fit failure.\n")
        next
    }

    # ---- 8. Exploratory diagnostics ----

    report_explore <- chlaa_fit_report(fit_explore, burnin = 0.25, thin = 2)
    cat("Exploratory acceptance rate:", report_explore$acceptance_rate, "\n")
    print(report_explore$posterior_summary)

    # Gate: skip production if exploratory acceptance is too low
    if (report_explore$acceptance_rate < 0.001) {
        cat("WARNING: Exploratory acceptance < 0.1% for", hz_name,
            "— skipping production fit.\n")
        rm(fit_explore, report_explore)
        next
    }

    # ---- 9. Warm start from exploratory posterior ----

    packer_ex <- attr(fit_explore, "packer")
    d <- length(packer_ex$names())
    n_samples <- dim(fit_explore$pars)[2]
    start_idx <- floor(0.25 * n_samples) + 1

    # Drop chain dimension if present (3D array from monty with 1 chain)
    pars_mat <- if (length(dim(fit_explore$pars)) == 3) {
        fit_explore$pars[, start_idx:n_samples, 1]
    } else {
        fit_explore$pars[, start_idx:n_samples]
    }
    pooled <- t(pars_mat)
    colnames(pooled) <- packer_ex$names()

    # Posterior median as warm-start values
    warm_vec <- apply(pooled, 2, median)
    pars_warm <- pars
    for (nm in names(warm_vec)) pars_warm[[nm]] <- warm_vec[[nm]]

    # Adaptive per-parameter production proposals from exploratory posterior.
    # Proposal variance = posterior variance, floored by the exploratory
    # proposal so parameters that barely moved still get reasonable steps.
    post_vars <- apply(pooled, 2, var)
    prod_proposal <- pmax(post_vars, explore_proposal)

    cat("\nWarm-start parameter values (posterior median):\n")
    for (nm in packer_ex$names()) {
        cat(sprintf("  %s = %.6f\n", nm, pars_warm[[nm]]))
    }
    cat("Production proposal_var (per-parameter):\n")
    print(setNames(prod_proposal, packer_ex$names()))

    rm(fit_explore, report_explore, packer_ex, pooled, pars_mat)

    # ---- 10. Production fit: 200 particles, 20000 steps, warm start ----

    cat("\n=== PRODUCTION FIT ===\n")
    fit <- tryCatch(
        chlaa_fit_pmcmc(
            data = fit_data,
            pars = pars_warm,
            n_particles = 200,
            n_steps = 20000,
            seed = 123,
            prior = custom_prior,
            packer = custom_packer(pars_warm),
            proposal_var = prod_proposal
        ),
        error = function(e) {
            cat("PRODUCTION FIT FAILED for", hz_name, ":", conditionMessage(e), "\n")
            NULL
        }
    )

    if (is.null(fit)) {
        cat("Skipping diagnostics/output for", hz_name,
            "due to production fit failure.\n")
        next
    }

    # ---- 11. Production diagnostics ----

    report_prod <- chlaa_fit_report(fit, burnin = 0.25, thin = 2)
    cat("Production acceptance rate:", report_prod$acceptance_rate, "\n")
    print(report_prod$posterior_summary)

    p_trace <- chlaa_plot_trace(fit, parameters = fitted_names, burnin = 0, thin = 1)
    ggsave(file.path(fig_dir, paste0(hz_name, "_production_trace.png")),
        p_trace, width = 14, height = 10, dpi = 150
    )

    p_lltrace <- chlaa_plot_likelihood_trace(fit, burnin = 0.25, thin = 2)
    ggsave(file.path(fig_dir, paste0(hz_name, "_production_likelihood_trace.png")),
        p_lltrace, width = 8, height = 5, dpi = 150
    )

    p_pairs <- chlaa_plot_parameter_pairs(fit, parameters = fitted_names, burnin = 0.25, thin = 2)
    ggsave(file.path(fig_dir, paste0(hz_name, "_production_pairs.png")),
        p_pairs, width = 12, height = 12, dpi = 150
    )

    # ---- 12. Posterior predictive projection ----
    # Forecast daily (zero_every = 1), then aggregate to weekly for plot.

    max_day <- as.integer(outbreak_end - outbreak_start)
    forecast_time <- seq(0, max_day, by = 1)

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

    # Aggregate daily forecast to weekly for comparison with IDSR data
    f_cases <- fc %>%
        filter(variable == "cases") %>%
        mutate(date = outbreak_start + time) %>%
        mutate(week_start = as.Date(cut(date, breaks = "7 days", start.on.monday = FALSE))) %>%
        group_by(week_start) %>%
        summarise(
            mean = sum(mean),
            q0p025 = sum(q0p025),
            q0p25 = sum(q0p25),
            q0p75 = sum(q0p75),
            q0p975 = sum(q0p975),
            .groups = "drop"
        ) %>%
        rename(date = week_start)

    # ---- 13. Production fit plot ----
    # Show 3 months of IDSR data before and after the outbreak window
    # for context. Model overlay only covers the outbreak period.

    plot_start <- outbreak_start - 90
    plot_end   <- outbreak_end + 90

    hz_data_wide <- hz_weekly %>%
        filter(date >= plot_start, date <= plot_end) %>%
        select(date, cases) %>%
        arrange(date)

    p_fit <- ggplot() +
        geom_col(
            data = hz_data_wide, aes(x = date, y = cases),
            fill = "grey70", width = 6, alpha = 0.6
        ) +
        geom_vline(xintercept = c(outbreak_start, outbreak_end),
            linetype = "dashed", colour = "grey40", linewidth = 0.5) +
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
            data = hz_data_wide, aes(x = date, y = cases, colour = "IDSR data"),
            size = 1.5
        ) +
        scale_colour_manual(values = c("IDSR data" = "black", "Mean fit" = "#2c7f62")) +
        scale_fill_manual(values = c("50% CI" = "#b88a66", "95% CI" = "#9a3b32")) +
        scale_x_date(date_breaks = "1 month", date_labels = "%Y-%m",
            limits = c(plot_start, plot_end)) +
        labs(
            title = paste0(hz_name, " — production fit (20000 steps, 200 particles, warm start, 4-par)"),
            subtitle = paste0("Dashed lines = outbreak window: ",
                format(outbreak_start, "%d %b %Y"), " – ",
                format(outbreak_end, "%d %b %Y")),
            x = NULL, y = "Cases/week",
            colour = NULL, fill = NULL
        ) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(fig_dir, paste0(hz_name, "_production_fit.png")),
        p_fit, width = 12, height = 6, dpi = 150
    )

    # ---- 14. Save fitted parameters for scenario analysis ----

    pars_fitted <- chlaa_update_from_fit(
        fit = fit, pars = pars, draw = "median", burnin = 0.25, thin = 2
    )
    saveRDS(pars_fitted, file.path(output_dir, paste0(hz_name, "_fitted_pars.rds")))
    cat("Saved fitted parameters to:", file.path(output_dir, paste0(hz_name, "_fitted_pars.rds")), "\n")

    # ---- Cleanup ----
    rm(
        fit, report_prod, fc, f_cases, p_trace, p_lltrace, p_pairs, p_fit,
        hz_weekly, hz_outbreak, hz_daily, hz_data_weekly, hz_data_wide, fit_data,
        pars, pars_warm, pars_fitted, warm_vec, post_vars, prod_proposal
    )
    gc()

    cat("\nFinished fitting", hz_name, "\n")
}

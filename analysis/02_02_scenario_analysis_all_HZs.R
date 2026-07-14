# =========================================================================
# Multi-Site Scenario Analysis for All Health Zones
# =========================================================================
#
# This script applies the scenario workflow from 03_01_scenario_workflow.R
# to all health zones using the fitted outputs from 01_02_fitting_all_HZs.R.
#
# For each HZ it:
#   1. Loads the pMCMC fit artifact from figures/.rds files/{hz_name}_fit.rds
#   2. Reconstructs observed data from the IDSR dataset
#   3. Defines three scenarios: no_intervention, aa_response, aa_response+vaccine
#   4. Runs posterior scenario forecasts (paired with baseline)
#   5. Saves 2 figures per HZ:
#      - scenarios_{hz_name}_absolute_forecasts_cases.png
#      - scenarios_{hz_name}_difference_cumulative_deaths.png
#   6. Saves decision summary tables (RDS + CSV) per HZ
#
# Designed for PBS array job submission (one HZ per job) or sequential
# interactive use. Mirrors the array job pattern of 01_02_fitting_all_HZs.R.
#
# =========================================================================

library(chlaa)
library(ggplot2)
library(dplyr)
library(tidyr)

# ---- Setup ----

data_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/analysis/data"
output_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/output"
fig_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/figures"
rds_dir <- file.path(fig_dir, ".rds files")
tables_dir <- file.path(fig_dir, "tables")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Core Scenario Function ----

run_scenarios_hz <- function(hz_name,
                             trigger_threshold = 50,
                             campaign_days = 28,
                             vax_coverage = 0.20,
                             n_draws = 60,
                             n_particles_summary = 50,
                             burnin = 0.25,
                             seed = 12,
                             verbose = TRUE) {
    if (verbose) cat("\n========================================\n")
    if (verbose) cat("Scenario analysis for:", hz_name, "\n")
    if (verbose) cat("========================================\n")

    # ---- 1. Load fit artifact ----

    fit_path <- file.path(rds_dir, sprintf("%s_fit.rds", hz_name))
    if (!file.exists(fit_path)) {
        stop(sprintf("Fit artifact not found: %s", fit_path))
    }

    fit_obj <- readRDS(fit_path)
    fit <- fit_obj$fit
    base_pars <- fit_obj$pars_warm
    outbreak_start <- fit_obj$outbreak_start
    outbreak_end <- fit_obj$outbreak_end

    if (verbose) cat("Outbreak window:", as.character(outbreak_start), "to", as.character(outbreak_end), "\n")

    # ---- 2. Reconstruct observed data ----

    idsr <- read.csv(file.path(data_dir, "IDSR_dataset.csv"))

    observed <- idsr %>%
        filter(hz == hz_name) %>%
        mutate(date = as.Date(date)) %>%
        filter(date >= outbreak_start, date <= outbreak_end) %>%
        arrange(date) %>%
        mutate(time = seq_len(n()) * 7L) %>%
        select(time, date, cases, deaths, population)

    if (nrow(observed) == 0) {
        stop(sprintf("No outbreak data found for %s", hz_name))
    }

    pop_hz <- observed$population[1]
    if (verbose) cat("Population:", pop_hz, "\n")
    if (verbose) cat("Outbreak weeks:", nrow(observed), "\n")
    if (verbose) cat("Total cases:", sum(observed$cases), "\n")

    # ---- 3. Update parameters from posterior ----

    pars_fit <- chlaa_update_from_fit(
        fit = fit,
        pars = base_pars,
        draw = "median",
        burnin = burnin
    )

    # ---- 4. Define scenario parameters ----

    horizon <- max(observed$time) + 182
    scenario_time <- seq(7, horizon, by = 7)
    response_end <- horizon + 1
    vaccine_doses <- floor(vax_coverage * pars_fit$N)

    # Find trigger time (first week with cases >= threshold)
    trigger_weeks <- observed$time[observed$cases >= trigger_threshold]

    if (length(trigger_weeks) == 0) {
        if (verbose) cat("WARNING: Trigger threshold (", trigger_threshold,
            " cases) never reached. Skipping vaccination scenario.\n")
        trigger_time <- NA_real_
        include_vax_scenario <- FALSE
    } else {
        trigger_time <- min(trigger_weeks)
        include_vax_scenario <- TRUE
        if (verbose) cat("Trigger time:", trigger_time, "(day)\n")
    }

    # ---- 5. Build vaccination schedule helper ----

    make_scenario_vax_schedule <- function(total_doses, start_day, end_day) {
        n_days <- max(end_day - start_day, 1L)
        daily_doses <- total_doses / n_days
        sched_time <- as.integer(c(start_day, end_day))
        sched_doses <- c(daily_doses, 0)
        if (min(sched_time) > 0) {
            sched_time <- c(0L, sched_time)
            sched_doses <- c(0, sched_doses)
        }
        list(
            vax1_schedule_time = sched_time,
            vax1_schedule_doses = sched_doses,
            n_vax1_schedule = length(sched_time)
        )
    }

    # Empty vaccination schedule arrays
    empty_vax <- list(
        vax1_schedule_time = c(0L, 1L),
        vax1_schedule_doses = c(0, 0),
        n_vax1_schedule = 2L,
        vax2_schedule_time = c(0L, 1L),
        vax2_schedule_doses = c(0, 0),
        n_vax2_schedule = 2L
    )

    # ---- 6. Define scenarios ----

    # No intervention counterfactual
    no_intervention <- c(list(
        chlor_start = 0, chlor_end = 0, chlor_effect = 0,
        hyg_start = 0, hyg_end = 0, hyg_effect = 0,
        lat_start = 0, lat_end = 0, lat_effect = 0,
        cati_start = 0, cati_end = 0, cati_effect = 0,
        orc_start = 0, orc_end = 0, orc_capacity = 0,
        ctc_start = 0, ctc_end = 0, ctc_capacity = 0,
        vax1_start = 0, vax1_end = 0, vax1_total_doses = 0,
        vax2_start = 0, vax2_end = 0, vax2_total_doses = 0
    ), empty_vax)

    scenarios <- list(
        chlaa_scenario("no_interventions", no_intervention)
    )

    if (!is.na(trigger_time)) {
        # Anticipatory-action response: WASH + treatment at trigger
        aa_response <- c(list(
            chlor_start = trigger_time, chlor_end = response_end,
            chlor_effect = pars_fit$chlor_effect,
            hyg_start = trigger_time, hyg_end = response_end,
            hyg_effect = pars_fit$hyg_effect,
            lat_start = trigger_time, lat_end = response_end,
            lat_effect = max(pars_fit$lat_effect, 0.10),
            cati_start = trigger_time, cati_end = response_end,
            cati_effect = pars_fit$cati_effect,
            orc_start = trigger_time, orc_end = response_end,
            orc_capacity = pars_fit$orc_capacity,
            ctc_start = trigger_time, ctc_end = response_end,
            ctc_capacity = pars_fit$ctc_capacity,
            vax1_start = 0, vax1_end = 0, vax1_total_doses = 0,
            vax2_start = 0, vax2_end = 0, vax2_total_doses = 0
        ), empty_vax)

        scenarios <- c(scenarios, list(
            chlaa_scenario("aa_response", aa_response)
        ))

        if (include_vax_scenario) {
            # Anticipatory-action + vaccination campaign
            vax1_start_day <- trigger_time + 14
            vax1_end_day <- vax1_start_day + campaign_days
            vax_sched <- make_scenario_vax_schedule(
                vaccine_doses, vax1_start_day, vax1_end_day
            )
            aa_vaccination <- modifyList(aa_response, c(list(
                vax1_start = vax1_start_day,
                vax1_end = vax1_end_day,
                vax1_total_doses = vaccine_doses
            ), vax_sched))

            scenarios <- c(scenarios, list(
                chlaa_scenario("aa_response_plus_vaccine", aa_vaccination)
            ))
        }
    }

    if (verbose) {
        cat("Scenarios:", paste(vapply(scenarios, `[[`, character(1), "name"), collapse = ", "), "\n")
    }

    # ---- 7. Posterior scenario forecasts ----

    if (verbose) cat("\n=== POSTERIOR SCENARIO FORECASTS ===\n")

    scenario_fc <- chlaa_forecast_scenarios_from_fit(
        fit = fit,
        pars = base_pars,
        scenarios = scenarios,
        baseline_name = "fitted_response",
        time = scenario_time,
        vars = c("inc_symptoms_weekly", "cum_symptoms", "cum_deaths"),
        include_cases = TRUE,
        obs_interval = 7,
        n_draws = n_draws,
        burnin = burnin,
        seed = seed,
        dt = 1
    )

    # ---- 8. Scenario figures ----

    # Figure 1: Absolute scenario forecasts (weekly cases)
    p_absolute <- chlaa_plot_scenario_forecasts(
        scenario_fc,
        var = "cases",
        type = "absolute",
        data = observed,
        data_y = "cases"
    ) +
        labs(title = sprintf("%s — Scenario forecasts (weekly cases)",
            tools::toTitleCase(gsub("_", " ", hz_name))))

    ggsave(
        file.path(fig_dir, sprintf("scenarios_%s_absolute_forecasts_cases.png", hz_name)),
        plot = p_absolute, width = 12, height = 7, dpi = 300
    )

    # Figure 2: Difference in cumulative deaths relative to baseline
    p_diff_deaths <- chlaa_plot_scenario_forecasts(
        scenario_fc,
        var = "cum_deaths",
        type = "difference",
        include_baseline = FALSE
    ) +
        labs(title = sprintf("%s — Cumulative deaths averted vs fitted response",
            tools::toTitleCase(gsub("_", " ", hz_name))))

    ggsave(
        file.path(fig_dir, sprintf("scenarios_%s_difference_cumulative_deaths.png", hz_name)),
        plot = p_diff_deaths, width = 12, height = 7, dpi = 300
    )

    if (verbose) cat("Figures saved.\n")

    # ---- 9. Decision summary ----

    if (verbose) cat("\n=== DECISION SUMMARY ===\n")

    scenario_runs <- chlaa_run_scenarios(
        pars = pars_fit,
        scenarios = c(list(chlaa_scenario("fitted_response", list())), scenarios),
        time = scenario_time,
        n_particles = n_particles_summary,
        dt = 1,
        seed = seed + 1
    )

    scenario_summary <- chlaa_scenario_summary(
        scenario_runs,
        baseline = "fitted_response",
        incidence_var = "inc_symptoms_weekly"
    )

    scenario_comparison <- chlaa_compare_scenarios(
        scenario_runs,
        baseline = "fitted_response",
        include_econ = TRUE,
        wtp = 1500
    )

    if (verbose) {
        cat("\nScenario Summary:\n")
        print(scenario_summary)
        cat("\nScenario Comparison:\n")
        print(scenario_comparison)
    }

    # ---- 10. Save outputs ----

    # RDS
    scenario_output <- list(
        hz_name = hz_name,
        scenario_summary = scenario_summary,
        scenario_comparison = scenario_comparison,
        scenario_forecasts = scenario_fc,
        scenarios_defined = vapply(scenarios, `[[`, character(1), "name"),
        trigger_time = trigger_time,
        trigger_threshold = trigger_threshold,
        vaccine_doses = vaccine_doses,
        outbreak_start = outbreak_start,
        outbreak_end = outbreak_end,
        timestamp = Sys.time()
    )
    saveRDS(scenario_output, file.path(rds_dir, sprintf("%s_scenarios.rds", hz_name)))

    # CSV
    summary_df <- as.data.frame(scenario_summary)
    summary_df$hz <- hz_name
    write.csv(summary_df,
        file.path(tables_dir, sprintf("%s_scenario_summary.csv", hz_name)),
        row.names = FALSE
    )

    comparison_df <- as.data.frame(scenario_comparison)
    comparison_df$hz <- hz_name
    write.csv(comparison_df,
        file.path(tables_dir, sprintf("%s_scenario_comparison.csv", hz_name)),
        row.names = FALSE
    )

    if (verbose) {
        cat("\nRDS saved to:", rds_dir, "\n")
        cat("\nTables saved to:", tables_dir, "\n")
        cat("Figures saved to:", fig_dir, "\n")
    }

    return(scenario_output)
}

# ---- Main Execution ----

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
    # Array job mode: run single HZ specified by argument
    hz_to_run <- args[1]

    cat("\n", rep("=", 60), "\n", sep = "")
    cat("Running scenario analysis for:", hz_to_run, "\n")
    cat(rep("=", 60), "\n\n", sep = "")

    result <- tryCatch(
        {
            run_scenarios_hz(hz_to_run, verbose = TRUE)
        },
        error = function(e) {
            cat("\nERROR in scenario analysis for", hz_to_run, ":\n")
            cat(conditionMessage(e), "\n")
            saveRDS(list(
                hz_name = hz_to_run,
                error = conditionMessage(e),
                timestamp = Sys.time()
            ), file.path(rds_dir, sprintf("%s_scenarios_FAILED.rds", hz_to_run)))
            NULL
        }
    )

    if (!is.null(result)) {
        cat("\nSUCCESS: Completed scenario analysis for", hz_to_run, "\n")
    } else {
        cat("\nFAILURE: Could not complete scenario analysis for", hz_to_run, "\n")
        quit(status = 1)
    }
} else {
    # Interactive mode: run all HZs sequentially
    cat("No HZ specified. Running in sequential mode for all HZs.\n")
    cat("For production use, submit as array job with HZ name as argument.\n\n")

    # Discover which HZs have completed fit artifacts (saved in rds_dir by 01_02)
    all_fits <- list.files(rds_dir, pattern = "_fit\\.rds$", full.names = FALSE)
    all_hzs <- sort(gsub("_fit\\.rds$", "", all_fits))

    # Exclude any failed fits
    failed_fits <- list.files(rds_dir, pattern = "_FAILED\\.rds$", full.names = FALSE)
    failed_hzs <- gsub("_FAILED\\.rds$", "", failed_fits)
    all_hzs <- setdiff(all_hzs, failed_hzs)

    if (length(all_hzs) == 0) {
        stop("No completed fit artifacts found in: ", rds_dir)
    }

    cat("Health zones with completed fits:", paste(all_hzs, collapse = ", "), "\n\n")

    results <- list()
    for (hz_name in all_hzs) {
        results[[hz_name]] <- tryCatch(
            {
                run_scenarios_hz(hz_name, verbose = TRUE)
            },
            error = function(e) {
                cat("\nERROR in scenario analysis for", hz_name, ":\n")
                cat(conditionMessage(e), "\n\n")
                NULL
            }
        )
    }

    # ---- Summary ----

    cat("\n", rep("=", 60), "\n", sep = "")
    cat("SCENARIO ANALYSIS SUMMARY\n")
    cat(rep("=", 60), "\n", sep = "")

    n_success <- sum(sapply(results, function(x) !is.null(x)))
    n_failed <- length(results) - n_success

    cat("Total HZs:", length(results), "\n")
    cat("Successful:", n_success, "\n")
    cat("Failed:", n_failed, "\n")

    if (n_failed > 0) {
        failed_hzs <- names(results)[sapply(results, is.null)]
        cat("Failed HZs:", paste(failed_hzs, collapse = ", "), "\n")
    }

    # Combine all summaries into a single CSV
    if (n_success > 0) {
        all_summaries <- list.files(tables_dir, pattern = "_scenario_summary\\.csv$", full.names = TRUE)
        if (length(all_summaries) > 0) {
            combined <- do.call(rbind, lapply(all_summaries, read.csv))
            write.csv(combined,
                file.path(tables_dir, "all_hz_scenario_summary.csv"),
                row.names = FALSE
            )
            cat("\nCombined summary saved to:", file.path(tables_dir, "all_hz_scenario_summary.csv"), "\n")
        }

        all_comparisons <- list.files(tables_dir, pattern = "_scenario_comparison\\.csv$", full.names = TRUE)
        if (length(all_comparisons) > 0) {
            combined_comp <- do.call(rbind, lapply(all_comparisons, read.csv))
            write.csv(combined_comp,
                file.path(tables_dir, "all_hz_scenario_comparison.csv"),
                row.names = FALSE
            )
            cat("Combined comparison saved to:", file.path(tables_dir, "all_hz_scenario_comparison.csv"), "\n")
        }
    }
}

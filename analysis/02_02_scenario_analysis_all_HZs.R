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

# ---- Global constants (must match 01_02 — needed by deserialized packer closures) ----

H_REF <- 1.0
POP_REF <- 516000
RR_FIXED <- 0.30

seed_state_names <- c("E0", "A0", "M0", "Sev0", "Mu0", "Mt0", "Sevu0", "Sevt0", "C0")

seed_state <- function(E0, p) {
    newI <- E0 / p$incubation_time
    nsymp <- newI * (1 - p$prop_asym)
    s <- list(
        E0    = E0,
        A0    = newI * p$prop_asym * p$duration_asym,
        M0    = nsymp * (1 - p$p_progress_severe) * p$time_to_next_stage,
        Sev0  = nsymp * p$p_progress_severe * p$time_to_next_stage,
        Mu0   = nsymp * (1 - p$p_progress_severe) * (1 - p$seek_mild) * p$duration_sym,
        Mt0   = nsymp * (1 - p$p_progress_severe) * p$seek_mild * p$duration_sym,
        Sevu0 = nsymp * p$p_progress_severe * (1 - p$seek_severe) * p$duration_sym,
        Sevt0 = nsymp * p$p_progress_severe * p$seek_severe * p$duration_sym
    )
    shed <- p$shed_asym * s$A0 +
        p$shed_mild * (s$M0 + s$Mu0) + p$shed_mild * p$treated_shed_mult_orc * s$Mt0 +
        p$shed_severe * (s$Sev0 + s$Sevu0) + p$shed_severe * p$treated_shed_mult_ctc * s$Sevt0
    s$C0 <- (shed / p$contam_scale / p$time_to_contaminate) * p$water_clearance_time
    s
}

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
        dt = 0.25
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
        dt = 0.25,
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

# =========================================================================
# Composite scenario figure across all health zones
# =========================================================================
#
# One figure, 3 scenario rows x 2 outcome columns (Cases | Deaths):
#   row 1  AA response
#   row 2  AA response + vaccination
#   row 3  No interventions
# Each facet holds one horizontal box+whisker (median, 50% and 95%
# uncertainty intervals) per health zone, showing excess cumulative
# cases/deaths relative to the fitted response (the "difference" scenario
# type is scenario-minus-baseline on paired posterior draws, so the fitted
# response is identically 0 and is shown as the dashed reference line)
# at ~365 days since the start of the modelled outbreak. Style mirrors
# build_excess_plot() in 02_01_scenario_workflow.R, rotated into a
# forest-plot layout (one row per HZ) like fitted_r0_density.png in
# 01_05_collect_fitted_parameters.R.
#
# Health zones that never crossed the response trigger (and therefore
# have no aa_response / aa_response_plus_vaccine scenario, only
# no_interventions) are excluded from the whole figure so every facet
# covers the same set of zones.

cat("\n", rep("=", 60), "\n", sep = "")
cat("Building composite all-HZ scenario figure\n")
cat(rep("=", 60), "\n\n", sep = "")

# Display-friendly health zone names (e.g. "ngiri_ngiri" -> "Ngiri Ngiri",
# "maluku_i" -> "Maluku I"). Capitalises every word, unlike
# tools::toTitleCase() which leaves single-letter words like "i" lowercase.
# Shared by both the composite scenario figure and the intervention
# contribution figures further below.
hz_label <- function(x) {
    x <- gsub("_", " ", x)
    gsub("(?<=^|\\s)([a-z])", "\\U\\1", x, perl = TRUE)
}

scenario_rds_files <- list.files(rds_dir, pattern = "_scenarios\\.rds$", full.names = TRUE)
scenario_rds_files <- scenario_rds_files[!grepl("FAILED", scenario_rds_files)]

if (length(scenario_rds_files) == 0) {
    cat("No *_scenarios.rds files found - skipping composite figure.\n")
} else {
    scenario_objs <- lapply(scenario_rds_files, readRDS)
    names(scenario_objs) <- vapply(scenario_objs, `[[`, character(1), "hz_name")

    # Keep only HZs where all three scenarios were defined (i.e. the
    # response trigger was reached, so aa_response/+vaccine exist)
    full_scenario_set <- c("no_interventions", "aa_response", "aa_response_plus_vaccine")
    keep_hz <- names(scenario_objs)[vapply(scenario_objs, function(x) {
        all(full_scenario_set %in% x$scenarios_defined)
    }, logical(1))]

    dropped_hz <- setdiff(names(scenario_objs), keep_hz)
    if (length(dropped_hz) > 0) {
        cat("Excluding HZs without a full scenario set (response trigger never reached):",
            paste(dropped_hz, collapse = ", "), "\n")
    }
    cat("HZs included in composite figure:", paste(keep_hz, collapse = ", "), "\n\n")

    if (length(keep_hz) == 0) {
        cat("No HZs have a full scenario set - skipping composite figure.\n")
    } else {
        # For each retained HZ, snap to the scenario_time grid point closest
        # to 365 days since outbreak start (grids are weekly and horizons
        # vary by HZ, so the exact snapped day can differ slightly by HZ)
        target_day <- 365

        all_fc <- lapply(keep_hz, function(hz) {
            fc <- scenario_objs[[hz]]$scenario_forecasts
            snap_time <- fc$time[which.min(abs(unique(fc$time) - target_day))]
            fc %>%
                filter(
                    type == "difference",
                    scenario %in% full_scenario_set,
                    variable %in% c("cum_symptoms", "cum_deaths"),
                    time == snap_time
                ) %>%
                mutate(hz = hz, snap_time = snap_time)
        })

        composite_dat <- bind_rows(all_fc)

        snap_report <- composite_dat %>% distinct(hz, snap_time)
        cat("Snapped time point (days since outbreak start) per HZ:\n")
        for (i in seq_len(nrow(snap_report))) {
            cat(sprintf("  %-14s day %s\n", snap_report$hz[i], snap_report$snap_time[i]))
        }

        # ---- Health zone ordering: by no-intervention excess cases,
        #      reused across every facet for comparability. Ascending order
        #      so the largest excess ends up plotted at the top. ----
        hz_rank <- composite_dat %>%
            filter(scenario == "no_interventions", variable == "cum_symptoms") %>%
            arrange(q0p5) %>%
            pull(hz)

        composite_dat <- composite_dat %>%
            mutate(hz = factor(hz, levels = hz_rank))

        # ---- Facet ordering/labels and scenario colours (matches
        #      02_01_scenario_workflow.R's palette) ----
        scenario_order_all <- c("aa_response", "aa_response_plus_vaccine", "no_interventions")
        scenario_facet_labels <- c(
            "aa_response"              = "AA response",
            "aa_response_plus_vaccine" = "AA response + vaccination",
            "no_interventions"         = "No interventions"
        )
        scenario_colours_all <- c(
            "aa_response"              = "#cb86ff",
            "aa_response_plus_vaccine" = "#88b517",
            "no_interventions"         = "#f7776d"
        )
        baseline_colour_all <- "#0abfc6"
        variable_facet_labels <- c(cum_symptoms = "Cases", cum_deaths = "Deaths")

        composite_dat <- composite_dat %>%
            mutate(
                scenario = factor(scenario, levels = scenario_order_all),
                variable = factor(variable, levels = c("cum_symptoms", "cum_deaths"))
            )

        # Numeric label: "median (lower to upper)", positioned just beyond
        # the whisker end on whichever side the bar points (mirrors
        # build_excess_plot()'s num_label/label_y/label_vjust logic, using
        # hjust instead of vjust since the boxes are now horizontal)
        composite_dat <- composite_dat %>%
            mutate(
                num_label   = paste0(round(q0p5), " (", round(q0p025), " to ", round(q0p975), ")"),
                label_x     = ifelse(q0p5 >= 0, q0p975, q0p025),
                label_hjust = ifelse(q0p5 >= 0, -0.08, 1.08)
            )

        # Pad each variable's (column's) x-range so the numeric labels have
        # room to sit beyond the whiskers without being clipped by the
        # panel edge. facet_grid(scales = "free_x") shares one x-scale down
        # each column, and unions in any data mapped to x within that
        # column - including these invisible padding points.
        x_pad <- composite_dat %>%
            group_by(variable) %>%
            summarise(
                lo = min(q0p025, na.rm = TRUE),
                hi = max(q0p975, na.rm = TRUE),
                .groups = "drop"
            ) %>%
            mutate(range = hi - lo) %>%
            tidyr::pivot_longer(c(lo, hi), names_to = "side", values_to = "x") %>%
            mutate(x = ifelse(side == "lo", x - range * 0.55, x + range * 0.55)) %>%
            mutate(
                scenario = factor("aa_response", levels = scenario_order_all),
                hz = factor(hz_rank[1], levels = hz_rank)
            )

        # "Fitted response" annotation next to the dashed baseline, shown
        # once (top-left facet) rather than in every panel
        label_df <- data.frame(
            scenario = factor("aa_response", levels = scenario_order_all),
            variable = factor("cum_symptoms", levels = c("cum_symptoms", "cum_deaths")),
            x = 0, y = Inf, label = "Fitted response "
        )

        trigger_threshold_used <- scenario_objs[[keep_hz[1]]]$trigger_threshold
        caption_txt <- if (length(dropped_hz) > 0) {
            sprintf(
                "Excludes %s: response trigger (%s+ cases) never reached, so no AA response scenario was modelled.",
                paste(hz_label(dropped_hz), collapse = ", "), trigger_threshold_used
            )
        } else {
            NULL
        }

        p_composite <- ggplot(composite_dat, aes(y = hz)) +
            geom_vline(xintercept = 0, linetype = "dashed", colour = baseline_colour_all,
                       linewidth = 0.8) +
            geom_blank(data = x_pad, aes(x = x, y = hz)) +
            geom_text(
                data = label_df, aes(x = x, y = y, label = label),
                inherit.aes = FALSE, hjust = 0, vjust = 1.6, size = 3.2,
                colour = baseline_colour_all, fontface = "bold", family = "Helvetica"
            ) +
            geom_crossbar(
                aes(x = q0p5, xmin = q0p25, xmax = q0p75, fill = scenario),
                orientation = "y", width = 0.65, alpha = 0.8, colour = "grey30",
                linewidth = 0.3, middle.linewidth = 0.6
            ) +
            geom_errorbar(
                aes(xmin = q0p025, xmax = q0p975),
                orientation = "y", width = 0.35, linewidth = 0.4, colour = "grey30"
            ) +
            geom_text(
                aes(x = label_x, label = num_label, hjust = label_hjust),
                size = 2.6, family = "Helvetica"
            ) +
            facet_grid(
                scenario ~ variable,
                labeller = labeller(scenario = scenario_facet_labels, variable = variable_facet_labels),
                scales = "free_x"
            ) +
            scale_y_discrete(labels = hz_label) +
            scale_fill_manual(values = scenario_colours_all, guide = "none") +
            labs(
                x = "Excess cumulative count (vs fitted response)",
                y = NULL,
                title = "Scenario impact across health zones",
                subtitle = sprintf(
                    "Excess cumulative cases and deaths ~365 days since outbreak start (n = %d health zones)",
                    length(keep_hz)
                ),
                caption = caption_txt
            ) +
            theme_minimal(base_family = "Helvetica", base_size = 12) +
            theme(
                panel.grid       = element_blank(),
                panel.background = element_rect(fill = "white", colour = "grey70"),
                panel.border     = element_rect(fill = NA, colour = "grey70", linewidth = 0.5),
                plot.background  = element_rect(fill = "white", colour = NA),
                strip.text       = element_text(face = "bold", size = 11),
                axis.ticks.x     = element_line(colour = "grey40"),
                axis.ticks.y     = element_blank(),
                plot.title       = element_text(face = "bold", size = 15),
                plot.subtitle    = element_text(size = 10, colour = "grey40"),
                plot.caption     = element_text(size = 8.5, colour = "grey40", hjust = 0),
                panel.spacing    = unit(1, "lines")
            )

        ggsave(file.path(fig_dir, "scenario_excess_all_hz.png"),
            plot = p_composite, width = 12, height = 10, dpi = 300)
        ggsave(file.path(fig_dir, "scenario_excess_all_hz.pdf"),
            plot = p_composite, width = 12, height = 10)

        cat("\nComposite figure saved to:\n")
        cat("  ", file.path(fig_dir, "scenario_excess_all_hz.png"), "\n")
        cat("  ", file.path(fig_dir, "scenario_excess_all_hz.pdf"), "\n")
    }
}

# =========================================================================
# Relative contribution of individual interventions to the fitted response
# =========================================================================
#
# "Add-one-in" design: for each of 6 intervention types (CTC, ORC, CATI,
# Hygiene, Chlorination, Vaccination) simulate a scenario that starts from
# a fully zeroed-out (no_interventions) parameter set and restores ONLY
# that one intervention's real historical timing/effect - read straight
# off the fit's base_pars, i.e. the true fitted schedule (e.g. the actual
# daily vaccine dose array, not a synthetic one). Each is compared back to
# that SAME zeroed no_interventions baseline, sharing posterior draws in
# one chlaa_forecast_scenarios_from_fit() call, so the resulting
# cumulative cases/deaths difference is a statistically valid isolated
# effect for that lever alone.
#
# Each lever's isolated cases/deaths averted (relative to no_interventions)
# is expressed as a % of the TOTAL cases/deaths averted by the full
# historical response - the no_interventions-vs-fitted_response difference
# already computed and saved in <hz>_scenarios.rds - i.e. "what share of
# the total historical intervention effect can be attributed to this one
# lever". Because levers are evaluated one at a time from a common zero
# baseline (not by sequentially removing them from the full response),
# their shares can overlap and are NOT expected to sum to 100%.
#
# Vaccination 2nd dose (never used in any zone) and Latrines (used in only
# 1/12 zones) are excluded - too sparse to compare meaningfully across HZs.

cat("\n", rep("=", 60), "\n", sep = "")
cat("Computing individual intervention contributions\n")
cat(rep("=", 60), "\n\n", sep = "")

intervention_defs <- list(
    ctc_only = c("ctc_start", "ctc_end", "ctc_capacity"),
    orc_only = c("orc_start", "orc_end", "orc_capacity"),
    cati_only = c("cati_start", "cati_end", "cati_effect"),
    hygiene_only = c("hyg_start", "hyg_end", "hyg_effect"),
    chlorination_only = c("chlor_start", "chlor_end", "chlor_effect"),
    vaccination_only = c(
        "vax1_start", "vax1_end", "vax1_total_doses",
        "n_vax1_schedule", "vax1_schedule_time", "vax1_schedule_doses"
    )
)

intervention_labels <- c(
    ctc_only = "CTC", orc_only = "ORC", cati_only = "CATI",
    hygiene_only = "Hygiene", chlorination_only = "Chlorination",
    vaccination_only = "Vaccination"
)
intervention_order <- names(intervention_defs)

# Fully zeroed intervention parameter set - the "no_interventions" baseline
# used as the reference point for every add-one-in scenario
zero_intervention_pars <- function(base_pars) {
    modifyList(base_pars, list(
        chlor_start = 0, chlor_end = 0, chlor_effect = 0,
        hyg_start = 0, hyg_end = 0, hyg_effect = 0,
        lat_start = 0, lat_end = 0, lat_effect = 0,
        cati_start = 0, cati_end = 0, cati_effect = 0,
        orc_start = 0, orc_end = 0, orc_capacity = 0,
        ctc_start = 0, ctc_end = 0, ctc_capacity = 0,
        vax1_start = 0, vax1_end = 0, vax1_total_doses = 0,
        vax1_schedule_time = c(0L, 1L), vax1_schedule_doses = c(0, 0), n_vax1_schedule = 2L,
        vax2_start = 0, vax2_end = 0, vax2_total_doses = 0,
        vax2_schedule_time = c(0L, 1L), vax2_schedule_doses = c(0, 0), n_vax2_schedule = 2L
    ))
}

# Computes (and caches to rds_dir) the 6 add-one-in scenario forecasts for
# one HZ. Set force = TRUE to recompute even if a cached file exists.
compute_intervention_contributions_hz <- function(hz_name, n_draws = 60, burnin = 0.25,
                                                   seed = 21, verbose = TRUE, force = FALSE) {
    out_path <- file.path(rds_dir, sprintf("%s_intervention_contributions.rds", hz_name))
    if (file.exists(out_path) && !force) {
        if (verbose) cat("Using cached contributions for", hz_name, "\n")
        return(readRDS(out_path))
    }

    if (verbose) cat("\n--- Intervention contributions:", hz_name, "---\n")

    fit_obj <- readRDS(file.path(rds_dir, sprintf("%s_fit.rds", hz_name)))
    fit <- fit_obj$fit
    base_pars <- fit_obj$pars_warm
    outbreak_start <- fit_obj$outbreak_start
    outbreak_end <- fit_obj$outbreak_end

    idsr <- read.csv(file.path(data_dir, "IDSR_dataset.csv"))
    observed <- idsr %>%
        filter(hz == hz_name) %>%
        mutate(date = as.Date(date)) %>%
        filter(date >= outbreak_start, date <= outbreak_end) %>%
        arrange(date) %>%
        mutate(time = seq_len(n()) * 7L)

    horizon <- max(observed$time) + 182
    scenario_time <- seq(7, horizon, by = 7)

    pars_zeroed <- zero_intervention_pars(base_pars)

    scenarios <- lapply(names(intervention_defs), function(nm) {
        fields <- intervention_defs[[nm]]
        modify <- setNames(lapply(fields, function(f) base_pars[[f]]), fields)
        chlaa_scenario(nm, modify)
    })

    contrib_fc <- chlaa_forecast_scenarios_from_fit(
        fit = fit,
        pars = pars_zeroed,
        scenarios = scenarios,
        baseline_name = "no_interventions",
        time = scenario_time,
        vars = c("cum_symptoms", "cum_deaths"),
        include_cases = FALSE,
        obs_interval = 7,
        n_draws = n_draws,
        burnin = burnin,
        seed = seed,
        dt = 0.25
    )

    result <- list(
        hz_name = hz_name,
        contribution_forecasts = contrib_fc,
        scenario_time = scenario_time,
        timestamp = Sys.time()
    )
    saveRDS(result, out_path)
    if (verbose) cat("Saved:", out_path, "\n")
    result
}

all_fit_hzs <- list.files(rds_dir, pattern = "_fit\\.rds$")
all_fit_hzs <- sort(gsub("_fit\\.rds$", "", all_fit_hzs[!grepl("comparative", all_fit_hzs)]))

contribution_objs <- list()
for (hz in all_fit_hzs) {
    contribution_objs[[hz]] <- tryCatch(
        compute_intervention_contributions_hz(hz),
        error = function(e) {
            cat("ERROR computing contributions for", hz, ":", conditionMessage(e), "\n")
            NULL
        }
    )
}
contribution_objs <- contribution_objs[!vapply(contribution_objs, is.null, logical(1))]

# ---- Build the relative-contribution figures (one per outcome) ----

if (length(contribution_objs) == 0) {
    cat("No intervention contribution results available - skipping figure.\n")
} else {
    # Re-derive per-HZ total "no interventions" excess (the full historical
    # response's total preventable burden) as the normalising denominator -
    # already computed and saved in <hz>_scenarios.rds by run_scenarios_hz()
    scenario_rds_files2 <- list.files(rds_dir, pattern = "_scenarios\\.rds$", full.names = TRUE)
    scenario_rds_files2 <- scenario_rds_files2[!grepl("FAILED", scenario_rds_files2)]
    scenario_objs2 <- lapply(scenario_rds_files2, readRDS)
    names(scenario_objs2) <- vapply(scenario_objs2, `[[`, character(1), "hz_name")

    denom_dat <- lapply(names(scenario_objs2), function(hz) {
        fc <- scenario_objs2[[hz]]$scenario_forecasts
        snap_time <- fc$time[which.min(abs(unique(fc$time) - 365))]
        fc %>%
            filter(
                type == "difference", scenario == "no_interventions",
                variable %in% c("cum_symptoms", "cum_deaths"), time == snap_time
            ) %>%
            transmute(hz = hz, variable, denom = q0p5)
    }) %>% bind_rows()

    cat("Denominator check - HZs with a small no-intervention excess (< 20):\n")
    print(denom_dat %>% filter(abs(denom) < 20))

    contrib_dat <- lapply(names(contribution_objs), function(hz) {
        fc <- contribution_objs[[hz]]$contribution_forecasts
        snap_time <- fc$time[which.min(abs(unique(fc$time) - 365))]
        fc %>%
            filter(
                type == "difference", scenario %in% names(intervention_defs),
                variable %in% c("cum_symptoms", "cum_deaths"), time == snap_time
            ) %>%
            mutate(hz = hz, snap_time = snap_time)
    }) %>% bind_rows()

    # Cases/deaths averted BY that lever alone = -(scenario - no_interventions).
    # Negating a set of quantiles swaps their order (q0p025 <-> q0p975, etc.)
    contrib_dat <- contrib_dat %>%
        mutate(
            averted_q0p025 = -q0p975, averted_q0p25 = -q0p75,
            averted_q0p5   = -q0p5,
            averted_q0p75  = -q0p25, averted_q0p975 = -q0p025
        ) %>%
        left_join(denom_dat, by = c("hz", "variable")) %>%
        mutate(
            pct_q0p025 = 100 * averted_q0p025 / denom,
            pct_q0p25  = 100 * averted_q0p25  / denom,
            pct_q0p5   = 100 * averted_q0p5   / denom,
            pct_q0p75  = 100 * averted_q0p75  / denom,
            pct_q0p975 = 100 * averted_q0p975 / denom,
            intervention = factor(intervention_labels[scenario], levels = unname(intervention_labels[intervention_order])),
            variable = factor(variable, levels = c("cum_symptoms", "cum_deaths")),
            num_label = paste0(round(pct_q0p5), "% (", round(pct_q0p025), " to ", round(pct_q0p975), ")"),
            # Anchored at the 25/75% box edge, not the 95% whisker end: a few
            # HZs (small effective population -> noisy add-one-in estimates)
            # have 95% intervals wide enough to overflow any sensible shared
            # x-axis, dragging their label off-panel with them. The box edge
            # stays comfortably inside the axis range in every case, while
            # num_label still reports the true 95% interval.
            label_x = ifelse(pct_q0p5 >= 0, pct_q0p75, pct_q0p25),
            label_hjust = ifelse(pct_q0p5 >= 0, -0.08, 1.08)
        )

    # HZ ordering: reuse the same "no-intervention excess" magnitude ranking
    # as the composite scenario figure, now over all 12 HZs
    hz_rank_all <- denom_dat %>%
        filter(variable == "cum_symptoms") %>%
        arrange(denom) %>%
        pull(hz)

    contrib_dat <- contrib_dat %>%
        mutate(hz = factor(hz, levels = hz_rank_all))

    build_contribution_plot <- function(var_name, plot_title, xlim_clip = NULL) {
        dat <- contrib_dat %>% filter(variable == var_name)

        n_clipped <- 0
        if (!is.null(xlim_clip)) {
            n_clipped <- sum(dat$pct_q0p025 < xlim_clip[1] | dat$pct_q0p975 > xlim_clip[2], na.rm = TRUE)
            x_pad <- data.frame(
                x = xlim_clip,
                hz = factor(hz_rank_all[1], levels = hz_rank_all),
                intervention = factor(levels(dat$intervention)[1], levels = levels(dat$intervention))
            )
            # Keep each label anchored just inside the clipped view even when
            # its true whisker end (label_x) falls outside it - otherwise the
            # text itself gets clipped/garbled at the panel edge. The printed
            # numbers (num_label) still show the true, unclamped values. When
            # clamped, the label must also grow INWARD (hjust flips) instead
            # of outward, or it would just push back past the same edge.
            inner_margin <- diff(xlim_clip) * 0.02
            dat <- dat %>%
                mutate(
                    label_x_clamped = ifelse(
                        pct_q0p5 >= 0,
                        pmin(label_x, xlim_clip[2] - inner_margin),
                        pmax(label_x, xlim_clip[1] + inner_margin)
                    ),
                    was_clamped = label_x_clamped != label_x,
                    label_hjust = ifelse(was_clamped, -label_hjust + 1, label_hjust),
                    label_x = label_x_clamped
                )
        } else {
            x_range <- range(c(dat$pct_q0p025, dat$pct_q0p975), na.rm = TRUE)
            x_pad_amt <- diff(x_range) * 0.4
            x_pad <- data.frame(
                x = c(x_range[1] - x_pad_amt, x_range[2] + x_pad_amt),
                hz = factor(hz_rank_all[1], levels = hz_rank_all),
                intervention = factor(levels(dat$intervention)[1], levels = levels(dat$intervention))
            )
        }

        caption_txt <- paste(
            "Excludes Latrines (used in 1/12 zones) and 2nd vaccine dose (used in 0/12 zones).",
            "Contribution = (add-one-in scenario) vs (no interventions), as % of (no interventions vs fitted response).\n",
            "KNOWN ISSUE: Vaccination currently shows 0% for every zone due to a bug in how vax1_start/vax1_end/vax1_total_doses",
            "were computed at the fitting stage (fixed in 01_02_fitting_all_HZs.R but not yet reflected in the saved fits) -\n",
            "re-run the fit for the 7 affected zones before trusting this panel."
        )
        if (n_clipped > 0) {
            caption_txt <- paste0(
                caption_txt, sprintf(
                    " %d interval(s) extend beyond the axis range shown (x-axis clipped to %g to %g%% for readability).",
                    n_clipped, xlim_clip[1], xlim_clip[2]
                )
            )
        }

        p <- ggplot(dat, aes(y = hz)) +
            geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
            geom_blank(data = x_pad, aes(x = x, y = hz)) +
            geom_crossbar(
                aes(x = pct_q0p5, xmin = pct_q0p25, xmax = pct_q0p75, fill = intervention),
                orientation = "y", width = 0.65, alpha = 0.8, colour = "grey30",
                linewidth = 0.3, middle.linewidth = 0.6
            ) +
            geom_errorbar(
                aes(xmin = pct_q0p025, xmax = pct_q0p975),
                orientation = "y", width = 0.35, linewidth = 0.4, colour = "grey30"
            ) +
            geom_label(
                aes(x = label_x, label = num_label, hjust = label_hjust),
                size = 2.3, family = "Helvetica", colour = "black",
                fill = "white", linewidth = 0, label.padding = unit(0.12, "lines")
            ) +
            facet_wrap(~intervention, ncol = 2) +
            scale_y_discrete(labels = hz_label) +
            scale_fill_brewer(palette = "Set2", guide = "none") +
            labs(
                x = "Cases/deaths averted by this intervention alone, as % of total averted by the full historical response",
                y = NULL,
                title = plot_title,
                subtitle = "Each lever simulated in isolation from a no-intervention baseline; shares can overlap and need not sum to 100%",
                caption = caption_txt
            )

        if (!is.null(xlim_clip)) p <- p + coord_cartesian(xlim = xlim_clip)

        p +
            theme_minimal(base_family = "Helvetica", base_size = 12) +
            theme(
                panel.grid       = element_blank(),
                panel.background = element_rect(fill = "white", colour = "grey70"),
                panel.border     = element_rect(fill = NA, colour = "grey70", linewidth = 0.5),
                plot.background  = element_rect(fill = "white", colour = NA),
                strip.text       = element_text(face = "bold", size = 11),
                axis.ticks.x     = element_line(colour = "grey40"),
                axis.ticks.y     = element_blank(),
                plot.title       = element_text(face = "bold", size = 15),
                plot.subtitle    = element_text(size = 9.5, colour = "grey40"),
                plot.caption     = element_text(size = 7.5, colour = "grey40", hjust = 0),
                panel.spacing    = unit(1, "lines")
            )
    }

    p_contrib_cases <- build_contribution_plot(
        "cum_symptoms", "Relative contribution of individual interventions - Cases",
        xlim_clip = c(-150, 250)
    )
    ggsave(file.path(fig_dir, "intervention_contribution_cases.png"), p_contrib_cases, width = 12, height = 10, dpi = 300)
    ggsave(file.path(fig_dir, "intervention_contribution_cases.pdf"), p_contrib_cases, width = 12, height = 10)

    p_contrib_deaths <- build_contribution_plot("cum_deaths", "Relative contribution of individual interventions - Deaths")
    ggsave(file.path(fig_dir, "intervention_contribution_deaths.png"), p_contrib_deaths, width = 12, height = 10, dpi = 300)
    ggsave(file.path(fig_dir, "intervention_contribution_deaths.pdf"), p_contrib_deaths, width = 12, height = 10)

    cat("\nIntervention contribution figures saved to:\n")
    cat("  ", file.path(fig_dir, "intervention_contribution_cases.png"), "\n")
    cat("  ", file.path(fig_dir, "intervention_contribution_deaths.png"), "\n")
}


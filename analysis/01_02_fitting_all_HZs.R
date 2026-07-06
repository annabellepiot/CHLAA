# =========================================================================
# Multi-Site PMCMC Fitting Pipeline for All Health Zones
# =========================================================================
#
# This script wraps the Kirotshe fitting methodology into a function that
# can be applied to all health zones independently. Each HZ is treated as
# a spatially independent outbreak, allowing natural parameter variation
# across geographies without forcing shrinkage to a global mean.
#
# Key features:
# - Dynamic E0 initialization based on outbreak size
# - Robust error handling for small/sparse outbreaks
# - Vaccination campaign handling (including multiple campaigns per HZ)
# - Adaptive proposal covariance with fallback for failed pilots
# - Individual output files per HZ to prevent cascade failures
#
# =========================================================================

library(chlaa)
library(ggplot2)
library(tidyverse)

# ---- Setup ----

data_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/analysis/data"
output_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/output"
fig_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/figures"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Shared plotting helper ----

plot_case_fit <- function(fit, observed, title, seed, burnin = 0.25) {
    fc <- chlaa_forecast_from_fit(
        fit = fit,
        time = observed$time,
        vars = "inc_symptoms_weekly",
        include_cases = TRUE,
        obs_model = "nbinom",
        n_draws = 200,
        burnin = burnin,
        seed = seed,
        dt = 1,
        deterministic = FALSE
    )

    fit_cases <- fc |>
        filter(variable == "cases") |>
        left_join(observed |> select(time, date), by = "time")

    ggplot() +
        geom_ribbon(
            data = fit_cases,
            aes(date, ymin = q0p025, ymax = q0p975),
            fill = "#6baed6",
            alpha = 0.25
        ) +
        geom_ribbon(
            data = fit_cases,
            aes(date, ymin = q0p25, ymax = q0p75),
            fill = "#6baed6",
            alpha = 0.45
        ) +
        geom_line(data = fit_cases, aes(date, q0p5), colour = "#08519c", linewidth = 0.8) +
        geom_point(data = observed, aes(date, cases), size = 1.6) +
        labs(x = NULL, y = "Weekly reported cases", title = title)
}

# ---- Core Fitting Function ----

fit_hz <- function(hz_name,
                   n_explore = 100,
                   n_explore_steps = 1000,
                   n_prod = 200,
                   n_prod_steps = 10000,
                   seed_explore = 42,
                   seed_prod = 123,
                   verbose = TRUE) {
    if (verbose) cat("\n========================================\n")
    if (verbose) cat("Fitting health zone:", hz_name, "\n")
    if (verbose) cat("========================================\n")

    # ---- 1. Load and prepare data ----

    hz_params_long <- read.csv(file.path(data_dir, "hz_parameters.csv"),
        stringsAsFactors = FALSE
    )

    idsr <- read.csv(file.path(data_dir, "IDSR_dataset.csv"))
    vax_dates <- read.csv(file.path(data_dir, "ocv_vaccination_dates.csv"),
        stringsAsFactors = FALSE
    )

    # ---- 2. Extract parameters for this HZ ----
    # Keep in long format to handle multiple vaccination campaigns

    hz_rows_long <- hz_params_long %>% filter(hz == hz_name)

    if (nrow(hz_rows_long) == 0) {
        stop(sprintf("Health zone '%s' not found in hz_parameters.csv", hz_name))
    }

    # Extract outbreak dates (should be consistent across all rows for this HZ)
    outbreak_start <- hz_rows_long %>%
        filter(parameter == "outbreak_start") %>%
        pull(value) %>%
        first() %>%
        as.Date()

    outbreak_end <- hz_rows_long %>%
        filter(parameter == "outbreak_end") %>%
        pull(value) %>%
        first() %>%
        as.Date()

    if (verbose) cat("Outbreak window:", as.character(outbreak_start), "to", as.character(outbreak_end), "\n")

    # ---- 3. Helper functions ----

    safe_date_to_day <- function(date_str, origin) {
        d <- as.Date(date_str, format = "%Y-%m-%d")
        if (is.na(d)) {
            return(0L)
        }
        as.integer(d - origin)
    }

    get_param <- function(param_name) {
        val <- hz_rows_long %>%
            filter(parameter == param_name) %>%
            pull(value) %>%
            first()
        if (is.null(val) || length(val) == 0) {
            return(NA)
        }
        val
    }

    # ---- 4. Extract intervention dates ----

    orc_start_day <- safe_date_to_day(get_param("orc_start"), outbreak_start)
    orc_end_day <- safe_date_to_day(get_param("orc_end"), outbreak_start)
    ctc_start_day <- safe_date_to_day(get_param("ctc_start"), outbreak_start)
    ctc_end_day <- safe_date_to_day(get_param("ctc_end"), outbreak_start)
    chlor_start_day <- safe_date_to_day(get_param("chlor_start"), outbreak_start)
    chlor_end_day <- safe_date_to_day(get_param("chlor_end"), outbreak_start)
    hyg_start_day <- safe_date_to_day(get_param("hyg_start"), outbreak_start)
    hyg_end_day <- safe_date_to_day(get_param("hyg_end"), outbreak_start)
    cati_start_day <- safe_date_to_day(get_param("cati_start"), outbreak_start)
    cati_end_day <- safe_date_to_day(get_param("cati_end"), outbreak_start)
    lat_start_day <- safe_date_to_day(get_param("lat_start"), outbreak_start)
    lat_end_day <- safe_date_to_day(get_param("lat_end"), outbreak_start)

    chlor_effect_val <- as.numeric(get_param("chlor_effect"))
    hyg_effect_val <- as.numeric(get_param("hyg_effect"))
    cati_effect_val <- as.numeric(get_param("cati_effect"))
    lat_effect_val <- as.numeric(get_param("lat_effect"))

    # ---- 5. Handle vaccination campaigns with delivery profiles ----
    # Load vaccination data for this HZ and generate schedules with
    # empirical delivery profile (30.5%, 37.7%, 22.7%, 7.4%, 1.4%, 0.3%)

    # Vaccination helper functions
    # Note: odin2 interpolate("constant") requires the time array to cover
    # the full simulation range, starting at or before time 0.
    generate_vax_schedule <- function(total_doses, start_date, end_date, outbreak_start) {
        profile <- c(0.305, 0.377, 0.227, 0.074, 0.014, 0.003)
        n_days <- as.integer(end_date - start_date) + 1L
        n_days <- max(n_days, 1L)
        # Stretch or truncate profile to campaign duration
        day_weights <- approx(seq(0, 1, length.out = length(profile)),
            profile,
            xout = seq(0, 1, length.out = n_days)
        )$y
        day_weights <- day_weights / sum(day_weights)
        daily_doses <- round(total_doses * day_weights)
        # Correct rounding residual
        daily_doses[which.max(daily_doses)] <- daily_doses[which.max(daily_doses)] +
            (total_doses - sum(daily_doses))
        start_day <- as.integer(start_date - outbreak_start)
        sched <- data.frame(
            time = seq(start_day, by = 1L, length.out = n_days),
            doses = as.numeric(daily_doses)
        )
        # Ensure schedule covers time 0 for interpolation
        if (min(sched$time) > 0L) {
            sched <- rbind(data.frame(time = 0L, doses = 0), sched)
        }
        sched
    }

    generate_empty_vax_schedule <- function(outbreak_start) {
        # Need at least 2 time points for odin2 interpolate
        data.frame(time = c(0L, 1L), doses = c(0, 0))
    }

    prepare_vax_arrays <- function(schedule) {
        list(time = as.integer(schedule$time), doses = as.numeric(schedule$doses))
    }

    # Extract vaccination data from ocv_vaccination_dates.csv
    # Handle name variations (underscore vs hyphen vs space)
    vax_hz <- vax_dates %>%
        filter(healthzone == hz_name |
            healthzone == gsub("_", " ", hz_name) |
            healthzone == gsub("_", "-", hz_name))

    # Extract vax1 campaign data with total doses
    vax1_total_doses_hz <- hz_rows_long %>%
        filter(parameter == "vax1_total_doses") %>%
        pull(value) %>%
        as.numeric() %>%
        first()

    vax1_campaigns <- vax_hz %>%
        filter(parameter %in% c("vax1_start", "vax1_end", "vax1_total_doses")) %>%
        distinct(dose_round, parameter, .keep_all = TRUE) %>%
        pivot_wider(names_from = parameter, values_from = value, id_cols = c(healthzone, dose_round))

    if (!("vax1_start" %in% names(vax1_campaigns))) {
        vax1_campaigns <- vax1_campaigns[0, ]
    } else {
        vax1_campaigns <- vax1_campaigns %>%
            filter(!is.na(vax1_start) & vax1_start != "NA")
    }

    # Generate vax1 schedule
    if (nrow(vax1_campaigns) > 0 && !is.na(vax1_total_doses_hz) && vax1_total_doses_hz > 0) {
        # Use first campaign in outbreak period
        vax1_camp <- vax1_campaigns[1, ]
        vax1_start_date <- as.Date(vax1_camp$vax1_start)
        vax1_end_date <- as.Date(vax1_camp$vax1_end)

        if (!is.na(vax1_start_date) && !is.na(vax1_end_date) &&
            vax1_start_date >= outbreak_start && vax1_start_date <= outbreak_end) {
            vax1_schedule <- generate_vax_schedule(
                total_doses = vax1_total_doses_hz,
                start_date = vax1_start_date,
                end_date = vax1_end_date,
                outbreak_start = outbreak_start
            )
            vax1_arrays <- prepare_vax_arrays(vax1_schedule)
            vax1_start_day <- min(vax1_schedule$time)
            vax1_end_day <- max(vax1_schedule$time) + 1

            if (verbose) {
                cat(sprintf("Vaccination campaign 1: days %d to %d\n", vax1_start_day, vax1_end_day))
                cat(sprintf("  Total doses: %d over %d days\n", vax1_total_doses_hz, nrow(vax1_schedule)))
                cat(sprintf(
                    "  Daily range: %.0f - %.0f doses\n",
                    min(vax1_schedule$doses), max(vax1_schedule$doses)
                ))
            }
        } else {
            vax1_schedule <- generate_empty_vax_schedule(outbreak_start)
            vax1_arrays <- prepare_vax_arrays(vax1_schedule)
            vax1_start_day <- 0L
            vax1_end_day <- 0L
        }
    } else {
        vax1_schedule <- generate_empty_vax_schedule(outbreak_start)
        vax1_arrays <- prepare_vax_arrays(vax1_schedule)
        vax1_start_day <- 0L
        vax1_end_day <- 0L
    }

    # Extract vax2 campaign data with total doses
    vax2_total_doses_hz <- hz_rows_long %>%
        filter(parameter == "vax2_total_doses") %>%
        pull(value) %>%
        as.numeric() %>%
        first()

    vax2_campaigns <- vax_hz %>%
        filter(parameter %in% c("vax2_start", "vax2_end", "vax2_total_doses")) %>%
        distinct(dose_round, parameter, .keep_all = TRUE) %>%
        pivot_wider(names_from = parameter, values_from = value, id_cols = c(healthzone, dose_round))

    if (!("vax2_start" %in% names(vax2_campaigns))) {
        vax2_campaigns <- vax2_campaigns[0, ]
    } else {
        vax2_campaigns <- vax2_campaigns %>%
            filter(!is.na(vax2_start) & vax2_start != "NA")
    }

    # Generate vax2 schedule
    if (nrow(vax2_campaigns) > 0 && !is.na(vax2_total_doses_hz) && vax2_total_doses_hz > 0) {
        vax2_camp <- vax2_campaigns[1, ]
        vax2_start_date <- as.Date(vax2_camp$vax2_start)
        vax2_end_date <- as.Date(vax2_camp$vax2_end)

        if (!is.na(vax2_start_date) && !is.na(vax2_end_date) &&
            vax2_start_date >= outbreak_start && vax2_start_date <= outbreak_end) {
            vax2_schedule <- generate_vax_schedule(
                total_doses = vax2_total_doses_hz,
                start_date = vax2_start_date,
                end_date = vax2_end_date,
                outbreak_start = outbreak_start
            )
            vax2_arrays <- prepare_vax_arrays(vax2_schedule)
            vax2_start_day <- min(vax2_schedule$time)
            vax2_end_day <- max(vax2_schedule$time) + 1

            if (verbose) {
                cat(sprintf("Vaccination campaign 2: days %d to %d\n", vax2_start_day, vax2_end_day))
                cat(sprintf("  Total doses: %d over %d days\n", vax2_total_doses_hz, nrow(vax2_schedule)))
            }
        } else {
            vax2_schedule <- generate_empty_vax_schedule(outbreak_start)
            vax2_arrays <- prepare_vax_arrays(vax2_schedule)
            vax2_start_day <- 0L
            vax2_end_day <- 0L
        }
    } else {
        vax2_schedule <- generate_empty_vax_schedule(outbreak_start)
        vax2_arrays <- prepare_vax_arrays(vax2_schedule)
        vax2_start_day <- 0L
        vax2_end_day <- 0L
    }

    # ---- 6. Prepare weekly case data ----

    hz_weekly <- idsr %>%
        filter(hz == hz_name) %>%
        mutate(
            date = as.Date(date),
            deaths = replace_na(deaths, 0L)
        ) %>%
        select(date, year, week, cases, deaths, population)

    hz_outbreak <- hz_weekly %>%
        filter(date >= outbreak_start, date <= outbreak_end) %>%
        arrange(date)

    if (nrow(hz_outbreak) == 0) {
        stop(sprintf("No outbreak data found for %s in the specified outbreak period", hz_name))
    }

    total_cases <- sum(hz_outbreak$cases)
    if (verbose) cat("Outbreak weeks:", nrow(hz_outbreak), "\n")
    if (verbose) cat("Total weekly cases:", total_cases, "\n")

    # Weekly time points at 7-day intervals (matching zero_every = 7 accumulator)
    hz_data_weekly <- hz_outbreak %>%
        mutate(time = seq_len(n()) * 7L) %>%
        select(time, date, cases, deaths) %>%
        arrange(time)

    natural_fit_names <- c("trans_prob", "reporting_rate", "obs_size", "E0", "seek_severe", "beta_p2p")
    fit_names <- c("log_trans_prob", "logit_reporting_rate", "log_obs_size", "log_E0", "logit_seek_severe", "log_beta_p2p")
    freeze_reporting_rate <- TRUE

    parameter_summary <- function(fit, burnin = 0.25) {
        chlaa_fit_trace(fit, burnin = burnin, scale = "natural") %>%
            select(chain, iteration, parameter, value) %>%
            group_by(parameter) %>%
            summarise(
                q025 = quantile(value, 0.025),
                median = median(value),
                q975 = quantile(value, 0.975),
                .groups = "drop"
            )
    }

    draws_wide <- function(fit, burnin = 0.25, scale = c("sampled", "natural")) {
        scale <- match.arg(scale)
        chlaa_fit_trace(fit, burnin = burnin, scale = scale) %>%
            pivot_wider(names_from = parameter, values_from = value) %>%
            arrange(chain, iteration)
    }

    chain_median_starts <- function(fit, template_pars, burnin = 0.25) {
        dr <- draws_wide(fit, burnin = burnin, scale = "sampled")
        packer <- attr(fit, "packer", exact = TRUE)
        fit_names_local <- packer$names()

        dr %>%
            group_split(chain) %>%
            map(function(d) {
                theta <- vapply(fit_names_local, function(nm) median(d[[nm]], na.rm = TRUE), numeric(1))
                out <- template_pars
                for (nm in names(theta)) out[[nm]] <- theta[[nm]]
                out
            })
    }

    # ---- 7. Dynamic E0 initialization ----
    # Scale E0 based on outbreak size to improve identifiability
    # For small outbreaks, use conservative seeding; for large outbreaks,
    # seed from pre-outbreak data

    pop_hz <- hz_weekly$population[1]
    seed_date <- outbreak_start - 14
    seed_row <- hz_weekly %>%
        filter(date <= seed_date) %>%
        arrange(desc(date)) %>%
        slice(1)

    expected_reporting_rate <- 0.3520

    if (total_cases < 50) {
        # Small outbreak: conservative seeding
        E0_val <- max(5, ceiling(mean(hz_outbreak$cases[1:min(3, nrow(hz_outbreak))]) / expected_reporting_rate))
        if (verbose) cat("Small outbreak detected. Using conservative E0 initialization.\n")
    } else if (nrow(seed_row) > 0 && seed_row$cases > 0) {
        # Large outbreak with pre-outbreak data
        E0_val <- ceiling(seed_row$cases / expected_reporting_rate)
    } else {
        # Fallback: use early outbreak data
        E0_val <- ceiling(max(5, hz_outbreak$cases[1]) / expected_reporting_rate)
    }

    E0_val <- max(5, E0_val) # Ensure minimum of 5 (avoids prior boundary)

    if (verbose) cat(sprintf("Initial seeding: E0=%d\n", E0_val))

    # ---- 8. Set starting parameters ----

    pars_args <- list(
        N = pop_hz,
        Sev0 = 0,
        E0 = E0_val,
        M0 = 0,
        immunity_asym = 280,
        beta_p2p = 0.05,
        contam_half_sat = 1.0,
        trans_prob = 0.003225,
        incubation_time = 4.845,
        duration_sym = 14.48,
        seek_mild = 0.1,
        seek_severe = 0.4,
        reporting_rate = 0.3520,
        fatality_treated = 0.0021,
        fatality_untreated = 0.5,
        obs_size = 30,
        death_reporting_rate = 0.5,
        obs_size_deaths = 1.0,
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

    # Add optional interventions only if active
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
    if (vax1_start_day > 0) {
        pars_args$vax1_start <- vax1_start_day
        pars_args$vax1_end <- vax1_end_day
        pars_args$vax1_total_doses <- sum(vax1_schedule$doses)
    }
    if (vax2_start_day > 0) {
        pars_args$vax2_start <- vax2_start_day
        pars_args$vax2_end <- vax2_end_day
        pars_args$vax2_total_doses <- sum(vax2_schedule$doses)
    }

    pars <- do.call(chlaa_parameters, pars_args)

    # Append vaccination schedule arrays (odin-level parameters, not in chlaa_parameter_info)
    pars$vax1_schedule_time <- vax1_arrays$time
    pars$vax1_schedule_doses <- vax1_arrays$doses
    pars$n_vax1_schedule <- length(vax1_arrays$time)
    pars$vax2_schedule_time <- vax2_arrays$time
    pars$vax2_schedule_doses <- vax2_arrays$doses
    pars$n_vax2_schedule <- length(vax2_arrays$time)

    add_transformed_values <- function(pars, freeze_reporting_rate = FALSE) {
        pars$log_trans_prob <- log(pars$trans_prob)
        if (!freeze_reporting_rate) pars$logit_reporting_rate <- qlogis(pars$reporting_rate)
        pars$log_obs_size <- log(pars$obs_size)
        pars$log_E0 <- log(pars$E0)
        pars$logit_seek_severe <- qlogis(pars$seek_severe)
        pars$log_beta_p2p <- log(pars$beta_p2p)
        pars
    }

    make_packer <- function(pars, freeze_reporting_rate = FALSE) {
        fit_names_use <- if (freeze_reporting_rate) {
            c("log_trans_prob", "log_obs_size", "log_E0", "logit_seek_severe", "log_beta_p2p")
        } else {
            fit_names
        }
        fixed <- pars[setdiff(names(pars), c(fit_names_use, natural_fit_names))]
        monty::monty_packer(
            scalar = fit_names_use,
            fixed = fixed,
            process = function(p) {
                out <- list(
                    trans_prob = exp(p$log_trans_prob),
                    obs_size = exp(p$log_obs_size),
                    E0 = exp(p$log_E0),
                    seek_severe = plogis(p$logit_seek_severe),
                    beta_p2p = exp(p$log_beta_p2p)
                )
                if (!freeze_reporting_rate) out$reporting_rate <- plogis(p$logit_reporting_rate)
                out
            }
        )
    }

    make_start <- function(trans_prob, reporting_rate, obs_size, E0,
                           seek_severe = 0.4, beta_p2p = 0.05,
                           freeze_reporting_rate = FALSE) {
        start_args <- c(
            list(
                N = pop_hz,
                Sev0 = 0,
                E0 = E0,
                M0 = 0,
                immunity_asym = 280,
                beta_p2p = beta_p2p,
                contam_half_sat = 1.0,
                trans_prob = trans_prob,
                incubation_time = 4.845,
                duration_sym = 14.48,
                seek_mild = 0.1,
                seek_severe = seek_severe,
                reporting_rate = reporting_rate,
                fatality_treated = 0.0021,
                fatality_untreated = 0.5,
                obs_size = obs_size,
                death_reporting_rate = 0.5,
                obs_size_deaths = 1.0,
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
            ),
            if (chlor_start_day > 0) {
                list(
                    chlor_start = chlor_start_day,
                    chlor_end = chlor_end_day,
                    chlor_effect = chlor_effect_val
                )
            } else {
                NULL
            },
            if (lat_start_day > 0) {
                list(
                    lat_start = lat_start_day,
                    lat_end = lat_end_day,
                    lat_effect = lat_effect_val
                )
            } else {
                NULL
            },
            if (vax1_start_day > 0) {
                list(
                    vax1_start = vax1_start_day,
                    vax1_end = vax1_end_day,
                    vax1_total_doses = sum(vax1_schedule$doses)
                )
            } else {
                NULL
            },
            if (vax2_start_day > 0) {
                list(
                    vax2_start = vax2_start_day,
                    vax2_end = vax2_end_day,
                    vax2_total_doses = sum(vax2_schedule$doses)
                )
            } else {
                NULL
            }
        )

        out <- do.call(chlaa_parameters, start_args)
        # Append vaccination schedule arrays (odin-level, not in chlaa_parameter_info)
        out$vax1_schedule_time <- vax1_arrays$time
        out$vax1_schedule_doses <- vax1_arrays$doses
        out$n_vax1_schedule <- length(vax1_arrays$time)
        out$vax2_schedule_time <- vax2_arrays$time
        out$vax2_schedule_doses <- vax2_arrays$doses
        out$n_vax2_schedule <- length(vax2_arrays$time)
        add_transformed_values(out, freeze_reporting_rate = freeze_reporting_rate)
    }

    fit_prior_full <- monty::monty_dsl(
        {
            log_trans_prob ~ Uniform(-9.210340, -4.605170)
            logit_reporting_rate ~ Uniform(-2.944439, 1.386294)
            log_obs_size ~ Uniform(0, 5.703782)
            log_E0 ~ Uniform(1.386294, 7.600902)
            logit_seek_severe ~ Uniform(-2.944439, 2.944439)
            log_beta_p2p ~ Uniform(-6.907755, -0.693147)
        },
        gradient = FALSE
    )

    fit_prior_freeze <- monty::monty_dsl(
        {
            log_trans_prob ~ Uniform(-9.210340, -4.605170)
            log_obs_size ~ Uniform(0, 5.703782)
            log_E0 ~ Uniform(1.386294, 7.600902)
            logit_seek_severe ~ Uniform(-2.944439, 2.944439)
            log_beta_p2p ~ Uniform(-6.907755, -0.693147)
        },
        gradient = FALSE
    )

    fit_starts_full <- list(
        make_start(0.003225, 0.3520, 30, E0_val, seek_severe = 0.4, beta_p2p = 0.05),
        make_start(1.6e-3, 0.12, 20, max(5, round(E0_val * 0.5)), seek_severe = 0.2, beta_p2p = 0.01),
        make_start(4.0e-4, 0.70, 200, max(5, round(E0_val * 1.5)), seek_severe = 0.6, beta_p2p = 0.15)
    )

    fit_starts_freeze <- list(
        make_start(0.003225, 0.35, 30, E0_val, seek_severe = 0.4, beta_p2p = 0.05, freeze_reporting_rate = TRUE),
        make_start(1.6e-3, 0.35, 20, max(5, round(E0_val * 0.5)), seek_severe = 0.2, beta_p2p = 0.01, freeze_reporting_rate = TRUE),
        make_start(4.0e-4, 0.35, 200, max(5, round(E0_val * 1.5)), seek_severe = 0.6, beta_p2p = 0.15, freeze_reporting_rate = TRUE)
    )

    fit_prior_stage1 <- if (freeze_reporting_rate) fit_prior_freeze else fit_prior_full
    fit_starts_stage1 <- if (freeze_reporting_rate) fit_starts_freeze else fit_starts_full
    fit_packer_stage1 <- make_packer(fit_starts_stage1[[1]], freeze_reporting_rate = freeze_reporting_rate)

    # ---- 9. Prepare fit data ----

    fit_data <- data.frame(time = hz_data_weekly$time, cases = hz_data_weekly$cases, deaths = hz_data_weekly$deaths)

    # ---- 10. Exploratory fit with robust covariance ----

    # Start with diagonal proposal (size depends on freeze_reporting_rate)
    n_params <- if (freeze_reporting_rate) 5 else 6
    explore_proposal <- matrix(0, n_params, n_params)

    if (freeze_reporting_rate) {
        # 5 parameters: log_trans_prob, log_obs_size, log_E0, logit_seek_severe, log_beta_p2p
        explore_proposal[1, 1] <- 0.02 # log_trans_prob
        explore_proposal[2, 2] <- 0.08 # log_obs_size
        explore_proposal[3, 3] <- 0.08 # log_E0
        explore_proposal[4, 4] <- 0.05 # logit_seek_severe
        explore_proposal[5, 5] <- 0.05 # log_beta_p2p
    } else {
        # 6 parameters: log_trans_prob, logit_reporting_rate, log_obs_size, log_E0, logit_seek_severe, log_beta_p2p
        explore_proposal[1, 1] <- 0.02 # log_trans_prob
        explore_proposal[2, 2] <- 0.05 # logit_reporting_rate
        explore_proposal[3, 3] <- 0.08 # log_obs_size
        explore_proposal[4, 4] <- 0.08 # log_E0
        explore_proposal[5, 5] <- 0.05 # logit_seek_severe
        explore_proposal[6, 6] <- 0.05 # log_beta_p2p
    }

    # Safeguard: ensure positive variances
    stopifnot(all(diag(explore_proposal) > 0))

    if (verbose) cat("\n=== EXPLORATORY FIT ===\n")

    # Wrap in tryCatch to handle potential failures in small outbreaks
    fit_explore <- tryCatch(
        {
            chlaa_fit_pmcmc(
                data = fit_data,
                pars = fit_starts_stage1[[1]],
                chain_pars = fit_starts_stage1,
                n_chains = length(fit_starts_stage1),
                n_particles = n_explore,
                n_steps = n_explore_steps,
                seed = seed_explore,
                prior = fit_prior_stage1,
                packer = fit_packer_stage1,
                proposal_var = explore_proposal,
                obs_interval = 7,
                time_start = 0
            )
        },
        error = function(e) {
            if (verbose) cat("Exploratory fit failed:", conditionMessage(e), "\n")
            if (verbose) cat("Attempting with reduced particles...\n")

            # Fallback: try with fewer particles
            chlaa_fit_pmcmc(
                data = fit_data,
                pars = fit_starts_stage1[[1]],
                chain_pars = fit_starts_stage1,
                n_chains = length(fit_starts_stage1),
                n_particles = max(n_explore, 200),
                n_steps = n_explore_steps,
                seed = seed_explore,
                prior = fit_prior_stage1,
                packer = fit_packer_stage1,
                proposal_var = explore_proposal,
                obs_interval = 7,
                time_start = 0
            )
        }
    )

    # ---- 11. Exploratory diagnostics ----

    report_explore <- chlaa_fit_report(fit_explore, burnin = 0.25, thin = 2)
    if (verbose) cat("Exploratory acceptance rate:", report_explore$acceptance_rate, "\n")
    if (verbose) print(report_explore$posterior_summary)

    # ---- 12. Learn covariance from exploratory run ----

    packer <- attr(fit_explore, "packer")
    d <- length(packer$names())
    n_samples <- dim(fit_explore$pars)[2]
    start_idx <- floor(0.25 * n_samples) + 1

    pars_mat <- if (length(dim(fit_explore$pars)) == 3) {
        fit_explore$pars[, start_idx:n_samples, 1]
    } else {
        fit_explore$pars[, start_idx:n_samples]
    }
    pooled <- t(pars_mat)
    colnames(pooled) <- packer$names()

    warm_vec <- apply(pooled, 2, median)
    pars_warm <- fit_starts_stage1[[1]]
    for (nm in names(warm_vec)) pars_warm[[nm]] <- warm_vec[[nm]]

    fit_starts_stage2 <- if (freeze_reporting_rate) {
        chain_median_starts(fit_explore, template_pars = fit_starts_full[[1]])
    } else {
        fit_starts_full
    }

    fit_prior_stage2 <- fit_prior_full
    fit_packer_stage2 <- make_packer(fit_starts_stage2[[1]], freeze_reporting_rate = FALSE)

    # Robust covariance calculation with fallback
    prod_proposal <- tryCatch(
        {
            cov_mat <- cov(pooled) * (2.38^2 / d)
            eig <- eigen(cov_mat, symmetric = TRUE)

            # Check for degenerate covariance (can happen in small outbreaks)
            if (any(eig$values < 1e-12) || any(!is.finite(eig$values))) {
                if (verbose) cat("Warning: Degenerate covariance matrix. Using diagonal fallback.\n")
                stop("Degenerate covariance")
            }

            eig$values <- pmax(eig$values, 1e-12)
            eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
        },
        error = function(e) {
            # Fallback to scaled diagonal from exploratory run
            if (verbose) cat("Using diagonal proposal (covariance learning failed).\n")
            diag(apply(pooled, 2, var)) * (2.38^2 / d)
        }
    )

    # If exploratory was 5-param (frozen reporting_rate) but production is 6-param,
    # expand the proposal matrix by inserting logit_reporting_rate at position 2
    if (freeze_reporting_rate && ncol(prod_proposal) == 5) {
        prop6 <- matrix(0, 6, 6)
        # Map: explore indices 1,2,3,4,5 -> production indices 1,3,4,5,6
        idx_map <- c(1L, 3L, 4L, 5L, 6L)
        prop6[idx_map, idx_map] <- prod_proposal
        # Default variance for logit_reporting_rate (position 2)
        prop6[2, 2] <- 0.05 * (2.38^2 / 6)
        prod_proposal <- prop6
    }

    if (verbose) {
        cat("\nWarm-start parameter values:\n")
        for (nm in packer$names()) {
            cat(sprintf("  %s = %.6f\n", nm, pars_warm[[nm]]))
        }
    }

    # Clean up exploratory objects
    rm(fit_explore, report_explore, packer, pooled, pars_mat)
    gc()

    # ---- 13. Production fit ----

    if (verbose) cat("\n=== PRODUCTION FIT ===\n")

    fit <- tryCatch(
        {
            chlaa_fit_pmcmc(
                data = fit_data,
                pars = pars_warm,
                chain_pars = fit_starts_stage2,
                n_chains = length(fit_starts_stage2),
                n_particles = n_prod,
                n_steps = n_prod_steps,
                seed = seed_prod,
                prior = fit_prior_stage2,
                packer = fit_packer_stage2,
                proposal_var = prod_proposal,
                obs_interval = 7,
                time_start = 0
            )
        },
        error = function(e) {
            if (verbose) cat("Production fit failed:", conditionMessage(e), "\n")
            if (verbose) cat("Attempting with reduced particles and steps...\n")

            # Fallback for small outbreaks
            chlaa_fit_pmcmc(
                data = fit_data,
                pars = pars_warm,
                chain_pars = fit_starts_stage2,
                n_chains = length(fit_starts_stage2),
                n_particles = max(n_prod, 200),
                n_steps = max(2000, n_prod_steps / 2),
                seed = seed_prod,
                prior = fit_prior_stage2,
                packer = fit_packer_stage2,
                proposal_var = prod_proposal,
                obs_interval = 7,
                time_start = 0
            )
        }
    )

    # ---- 15. Production diagnostics ----

    packer_prod <- attr(fit, "packer")
    fitted_names <- packer_prod$names()

    report_prod <- chlaa_fit_report(fit, burnin = 0.25, thin = 2)
    if (verbose) cat("Production acceptance rate:", report_prod$acceptance_rate, "\n")
    if (verbose) print(report_prod$posterior_summary)

    # Check for poor mixing or identifiability issues
    if (report_prod$acceptance_rate < 0.10) {
        if (verbose) cat("WARNING: Low acceptance rate may indicate identifiability issues.\n")
    }

    # ---- 15. Generate diagnostic plots (without printing to avoid duplication) ----

    p_trace <- chlaa_plot_trace(
        fit,
        parameters = natural_fit_names,
        burnin = 0.25,
        scale = "natural"
    )
    ggsave(file.path(fig_dir, sprintf("diagnosis_%s_production_trace.png", hz_name)),
        p_trace,
        width = 12, height = 8, dpi = 300
    )

    p_lltrace <- chlaa_plot_likelihood_trace(fit, burnin = 0.25, thin = 2)
    ggsave(file.path(fig_dir, sprintf("diagnosis_%s_production_likelihood_trace.png", hz_name)),
        p_lltrace,
        width = 12, height = 8, dpi = 300
    )

    p_pairs <- chlaa_plot_parameter_pairs(
        fit,
        parameters = natural_fit_names,
        burnin = 0.25,
        scale = "natural",
        max_points = 3000
    )

    ggsave(
        file.path(fig_dir, sprintf("diagnosis_%s_production_pairs.png", hz_name)),
        plot = p_pairs,
        width = 10,
        height = 10,
        dpi = 300
    )

    p_dist <- chlaa_plot_parameter_distributions(
        fit,
        parameters = natural_fit_names,
        burnin = 0.25,
        scale = "natural"
    )
    ggsave(file.path(fig_dir, sprintf("diagnosis_%s_production_distributions.png", hz_name)),
        p_dist,
        width = 12, height = 8, dpi = 300
    )

    # ---- 16. Posterior predictive forecast ----

    n_weeks <- nrow(hz_data_weekly)

    # ---- 18. Fit plot ----

    p_fit <- plot_case_fit(
        fit,
        hz_data_weekly,
        sprintf("%s weekly reported cases", tools::toTitleCase(hz_name)),
        seed = seed_prod
    )

    ggsave(file.path(fig_dir, sprintf("diagnosis_%s_production_fit.png", hz_name)),
        p_fit,
        width = 10, height = 6, dpi = 300
    )

    # ---- 19. Save fit object ----

    fit_output <- list(
        hz_name = hz_name,
        fit = fit,
        pars_start = fit_starts_stage2,
        pars_warm = pars_warm,
        report = report_prod,
        parameter_summary = parameter_summary(fit),
        outbreak_start = outbreak_start,
        outbreak_end = outbreak_end,
        total_cases = total_cases,
        n_weeks = n_weeks,
        fitted_parameters = natural_fit_names,
        timestamp = Sys.time()
    )

    saveRDS(fit_output, file.path(output_dir, sprintf("%s_fit.rds", hz_name)))

    if (verbose) cat("\nFit saved to:", file.path(output_dir, sprintf("%s_fit.rds", hz_name)), "\n")
    if (verbose) cat("Figures saved to:", fig_dir, "\n")

    return(fit_output)
}

# ---- Main Execution ----

# Get command line arguments for array job submission
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
    # Array job mode: fit single HZ specified by argument
    hz_to_fit <- args[1]

    cat("\n", rep("=", 60), "\n", sep = "")
    cat("Running array job for:", hz_to_fit, "\n")
    cat(rep("=", 60), "\n\n", sep = "")

    result <- tryCatch(
        {
            fit_hz(hz_to_fit, verbose = TRUE)
        },
        error = function(e) {
            cat("\nERROR fitting", hz_to_fit, ":\n")
            cat(conditionMessage(e), "\n")
            saveRDS(list(
                hz_name = hz_to_fit,
                error = conditionMessage(e),
                timestamp = Sys.time()
            ), file.path(output_dir, sprintf("%s_FAILED.rds", hz_to_fit)))
            NULL
        }
    )

    if (!is.null(result)) {
        cat("\nSUCCESS: Completed fitting for", hz_to_fit, "\n")
    } else {
        cat("\nFAILURE: Could not fit", hz_to_fit, "\n")
        quit(status = 1)
    }
} else {
    # Interactive mode: fit all HZs sequentially (for testing)
    cat("No HZ specified. Running in test mode for all HZs sequentially.\n")
    cat("For production use, submit as array job with HZ name as argument.\n\n")

    hz_params_long <- read.csv(file.path(data_dir, "hz_parameters.csv"),
        stringsAsFactors = FALSE
    )
    all_hzs <- sort(unique(hz_params_long$hz))

    cat("Health zones to fit:", paste(all_hzs, collapse = ", "), "\n\n")

    results <- list()
    for (hz_name in all_hzs) {
        results[[hz_name]] <- tryCatch(
            {
                fit_hz(hz_name, verbose = TRUE)
            },
            error = function(e) {
                cat("\nERROR fitting", hz_name, ":\n")
                cat(conditionMessage(e), "\n\n")
                NULL
            }
        )
    }

    # Summary
    cat("\n", rep("=", 60), "\n", sep = "")
    cat("FITTING SUMMARY\n")
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
}

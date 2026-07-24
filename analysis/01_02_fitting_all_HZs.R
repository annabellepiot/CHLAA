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
# - Dynamic E0 initialization with seed_state quasi-equilibrium
# - frac_neff: effective population fraction (N_eff = frac * pop_hz)
# - contam_half_sat decoupled from N_eff (uses census pop)
# - reporting_rate fixed at 0.30, not fitted
# - R0-based starting points
# - Proper R-hat / ESS diagnostics via posterior package
# - Individual output files per HZ to prevent cascade failures
#
# =========================================================================

library(chlaa)
library(ggplot2)
library(tidyverse)
library(posterior)

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

# ---- Global constants ----

H_REF <- 1.0
POP_REF <- 516000
RR_FIXED <- 0.30

# ---- Quasi-equilibrium seeding ----

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
        dt = 0.25,
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
                   variance_check_reps = 20,
                   variance_target = 2, # target Var[log-lik], not sd
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

    hz_rows_long <- hz_params_long %>% filter(hz == hz_name)

    if (nrow(hz_rows_long) == 0) {
        stop(sprintf("Health zone '%s' not found in hz_parameters.csv", hz_name))
    }

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

    # Burn-in: model starts this many days before first observation (day 7)
    time_start <- -21L

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

    generate_vax_schedule <- function(total_doses, start_date, end_date, outbreak_start) {
        profile <- c(0.305, 0.377, 0.227, 0.074, 0.014, 0.003)
        n_days <- as.integer(end_date - start_date) + 1L
        n_days <- max(n_days, 1L)
        day_weights <- approx(seq(0, 1, length.out = length(profile)),
            profile,
            xout = seq(0, 1, length.out = n_days)
        )$y
        day_weights <- day_weights / sum(day_weights)
        daily_doses <- round(total_doses * day_weights)
        daily_doses[which.max(daily_doses)] <- daily_doses[which.max(daily_doses)] +
            (total_doses - sum(daily_doses))
        start_day <- as.integer(start_date - outbreak_start)
        sched <- data.frame(
            time = seq(start_day, by = 1L, length.out = n_days),
            doses = as.numeric(daily_doses)
        )
        if (min(sched$time) > time_start) {
            sched <- rbind(data.frame(time = time_start, doses = 0), sched)
        }
        sched
    }

    generate_empty_vax_schedule <- function() {
        data.frame(time = c(time_start, time_start + 1L), doses = c(0, 0))
    }

    prepare_vax_arrays <- function(schedule) {
        list(time = as.integer(schedule$time), doses = as.numeric(schedule$doses))
    }

    vax_hz <- vax_dates %>%
        filter(healthzone == hz_name |
            healthzone == gsub("_", " ", hz_name) |
            healthzone == gsub("_", "-", hz_name))

    # vax1
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

    if (nrow(vax1_campaigns) > 0 && !is.na(vax1_total_doses_hz) && vax1_total_doses_hz > 0) {
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
            # NOTE: min(vax1_schedule$time) would incorrectly pick up the
            # zero-dose interpolation anchor point that generate_vax_schedule()
            # prepends at time_start (see line ~219-221) when the real campaign
            # starts after time_start - that anchor is only there so the odin
            # model has an interpolation point before the campaign, not a
            # signal that dosing began there. Using the first day with
            # doses > 0 gives the true campaign start (this previously made
            # vax1_start_day negative for every HZ with a real campaign,
            # which failed the `vax1_start_day > 0` gate below and silently
            # zeroed out vax1_start/vax1_end/vax1_total_doses - the fields
            # the simulator actually uses to switch vaccination on).
            vax1_start_day <- min(vax1_schedule$time[vax1_schedule$doses > 0])
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
            vax1_schedule <- generate_empty_vax_schedule()
            vax1_arrays <- prepare_vax_arrays(vax1_schedule)
            vax1_start_day <- 0L
            vax1_end_day <- 0L
        }
    } else {
        vax1_schedule <- generate_empty_vax_schedule()
        vax1_arrays <- prepare_vax_arrays(vax1_schedule)
        vax1_start_day <- 0L
        vax1_end_day <- 0L
    }

    # vax2
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
            # See the matching NOTE above vax1_start_day: exclude the
            # zero-dose interpolation anchor point when finding the true
            # campaign start.
            vax2_start_day <- min(vax2_schedule$time[vax2_schedule$doses > 0])
            vax2_end_day <- max(vax2_schedule$time) + 1

            if (verbose) {
                cat(sprintf("Vaccination campaign 2: days %d to %d\n", vax2_start_day, vax2_end_day))
                cat(sprintf("  Total doses: %d over %d days\n", vax2_total_doses_hz, nrow(vax2_schedule)))
            }
        } else {
            vax2_schedule <- generate_empty_vax_schedule()
            vax2_arrays <- prepare_vax_arrays(vax2_schedule)
            vax2_start_day <- 0L
            vax2_end_day <- 0L
        }
    } else {
        vax2_schedule <- generate_empty_vax_schedule()
        vax2_arrays <- prepare_vax_arrays(vax2_schedule)
        vax2_start_day <- 0L
        vax2_end_day <- 0L
    }

    # ---- 6. Prepare weekly case data ----

    hz_weekly <- idsr %>%
        filter(hz == hz_name) %>%
        mutate(
            date = as.Date(date)
        ) %>%
        select(date, year, week, cases, deaths, population)

    hz_outbreak <- hz_weekly %>%
        filter(date >= outbreak_start, date <= outbreak_end) %>%
        arrange(date)

    if (nrow(hz_outbreak) == 0) {
        stop(sprintf("No outbreak data found for %s in the specified outbreak period", hz_name))
    }

    pop_hz <- hz_weekly$population[1]
    total_cases <- sum(hz_outbreak$cases)
    if (verbose) cat("Population:", pop_hz, "\nOutbreak weeks:", nrow(hz_outbreak), "\nTotal cases:", total_cases, "\n")

    hz_data_weekly <- hz_outbreak %>%
        mutate(time = seq_len(n()) * 7L) %>%
        select(time, date, cases, deaths) %>%
        arrange(time)

    # ---- 7. Fitting infrastructure ----

    natural_fit_names <- c("trans_prob", "obs_size", "E0", "frac_neff")
    fit_names <- c("log_trans_prob", "log_obs_size", "log_E0", "logit_frac_neff")

    add_transformed_values <- function(pars) {
        pars$log_trans_prob <- log(pars$trans_prob)
        pars$log_obs_size <- log(pars$obs_size)
        pars$log_E0 <- log(pars$E0)
        pars$logit_frac_neff <- qlogis(pars$frac_neff)
        pars
    }

    make_packer <- function(pars) {
        # Capture globals locally so closure survives saveRDS/readRDS
        h_ref <- H_REF
        pop_ref <- POP_REF
        rr_fixed <- RR_FIXED
        pop <- pop_hz
        seed_fn <- seed_state

        fixed <- pars[setdiff(names(pars), c(
            fit_names, natural_fit_names,
            seed_state_names, "N", "contam_half_sat", "reporting_rate"
        ))]
        monty::monty_packer(
            scalar = fit_names,
            fixed = fixed,
            process = function(p) {
                frac <- plogis(p$logit_frac_neff)
                N_eff <- frac * pop
                out <- c(
                    list(
                        trans_prob = exp(p$log_trans_prob),
                        obs_size = exp(p$log_obs_size),
                        frac_neff = frac,
                        N = N_eff,
                        contam_half_sat = h_ref * (pop / pop_ref),
                        reporting_rate = rr_fixed
                    ),
                    seed_fn(exp(p$log_E0), p)
                )
                out
            }
        )
    }

    make_start <- function(trans_prob, obs_size, E0, frac_neff = 0.10) {
        N_eff <- frac_neff * pop_hz
        h <- H_REF * (pop_hz / POP_REF)
        pars_args <- list(
            N = N_eff, E0 = E0, M0 = 0, Sev0 = 0,
            immunity_asym = 280, contact_rate = 0, contam_half_sat = h,
            trans_prob = trans_prob, incubation_time = 4.845, duration_sym = 14.48,
            seek_mild = 0.1, seek_severe = 0.85,
            reporting_rate = RR_FIXED, fatality_treated = 0.0021, fatality_untreated = 0.5,
            obs_size = obs_size, death_reporting_rate = 0.5, obs_size_deaths = 1.0,
            orc_start = orc_start_day, orc_end = orc_end_day,
            ctc_start = ctc_start_day, ctc_end = ctc_end_day,
            hyg_start = hyg_start_day, hyg_end = hyg_end_day, hyg_effect = hyg_effect_val,
            cati_start = cati_start_day, cati_end = cati_end_day, cati_effect = cati_effect_val
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
        out <- do.call(chlaa_parameters, pars_args)
        ss <- seed_state(E0, out)
        for (nm in names(ss)) out[[nm]] <- ss[[nm]]
        out$frac_neff <- frac_neff
        out$vax1_schedule_time <- vax1_arrays$time
        out$vax1_schedule_doses <- vax1_arrays$doses
        out$n_vax1_schedule <- length(vax1_arrays$time)
        out$vax2_schedule_time <- vax2_arrays$time
        out$vax2_schedule_doses <- vax2_arrays$doses
        out$n_vax2_schedule <- length(vax2_arrays$time)
        add_transformed_values(out)
    }

    draws_wide <- function(fit, burnin = 0.25, scale = c("sampled", "natural")) {
        scale <- match.arg(scale)
        chlaa_fit_trace(fit, burnin = burnin, scale = scale) %>%
            pivot_wider(names_from = parameter, values_from = value) %>%
            arrange(chain, iteration)
    }

    chain_median_starts <- function(fit, template_pars, burnin = 0.25) {
        dr <- draws_wide(fit, burnin = burnin, scale = "sampled")
        packer_local <- attr(fit, "packer", exact = TRUE)
        fit_names_local <- packer_local$names()
        dr %>%
            group_split(chain) %>%
            map(function(d) {
                theta <- vapply(fit_names_local, function(nm) median(d[[nm]], na.rm = TRUE), numeric(1))
                out <- template_pars
                for (nm in names(theta)) out[[nm]] <- theta[[nm]]
                out
            })
    }

    # ---- 8. E0 initialization ----

    expected_reporting_rate <- 0.10
    seed_date <- outbreak_start - 14
    seed_row <- hz_weekly %>%
        filter(date <= seed_date) %>%
        arrange(desc(date)) %>%
        slice(1)

    if (total_cases < 50) {
        E0_val <- max(5, ceiling(mean(hz_outbreak$cases[1:min(3, nrow(hz_outbreak))]) / expected_reporting_rate))
        if (verbose) cat("Small outbreak detected. Using conservative E0 initialization.\n")
    } else if (nrow(seed_row) > 0 && seed_row$cases > 0) {
        E0_val <- ceiling(seed_row$cases / expected_reporting_rate)
    } else {
        E0_val <- ceiling(max(5, hz_outbreak$cases[1]) / expected_reporting_rate)
    }

    E0_MAX <- 800
    E0_val <- min(E0_val, 0.9 * E0_MAX) # start at most 720, strictly inside
    E0_val <- max(10, E0_val)

    if (verbose) {
        cat(sprintf("Population: %d\n", pop_hz))
        cat(sprintf(
            "contam_half_sat (census): %.6f  (pop ratio: %.3f)\n",
            H_REF * (pop_hz / POP_REF), pop_hz / POP_REF
        ))
        cat(sprintf("Initial seeding: E0=%d\n", E0_val))
    }

    # ---- 9. Prior (4 params, reporting_rate fixed at %.2f) ----

    fit_prior <- monty::monty_dsl(
        {
            log_trans_prob ~ Uniform(-9.21034, -2.995732) # log(c(1e-4, 5e-2))
            log_obs_size ~ Uniform(0, 5.703782)
            log_E0 ~ Uniform(2.302585, 6.684612) # log(c(10, 800))
            logit_frac_neff ~ Uniform(-4.6, 2.944439) # ~qlogis(c(0.01, 0.95))
        },
        gradient = FALSE
    )

    # ---- 10. R0-based starting points ----
    # trans_prob is exactly linear in R0 (contact_rate = 0 throughout this
    # pipeline), so a unit-trans_prob probe through chlaa_r0() gives an exact
    # per-chain inversion without needing a separate closed-form constant.

    r0_targets <- c(1.5, 2.5, 4.0)
    frac_starts <- c(0.10, 0.05, 0.20)
    tp_starts <- vapply(seq_along(frac_starts), function(i) {
        N_i <- frac_starts[i] * pop_hz
        h_i <- H_REF * (pop_hz / POP_REF)
        probe_pars <- chlaa_parameters(
            N = N_i, contact_rate = 0, contam_half_sat = h_i,
            incubation_time = 4.845, duration_sym = 14.48,
            seek_mild = 0.1, seek_severe = 0.85, trans_prob = 1
        )
        r0_targets[i] / chlaa_r0(probe_pars)
    }, numeric(1))

    fit_starts <- list(
        make_start(tp_starts[1], 30, E0_val, frac_neff = frac_starts[1]),
        make_start(tp_starts[2], 20, max(10, round(E0_val * 0.5)), frac_neff = frac_starts[2]),
        make_start(tp_starts[3], 100, min(0.9 * E0_MAX, max(10, round(E0_val * 1.5))), frac_neff = frac_starts[3])
    )

    fit_packer_stage1 <- make_packer(fit_starts[[1]])

    # ---- 11. Prepare fit data ----

    fit_data <- data.frame(
        time = hz_data_weekly$time,
        cases = hz_data_weekly$cases,
        deaths = hz_data_weekly$deaths
    )

    # ---- 12. Exploratory fit ----

    explore_proposal <- matrix(0, 4, 4)
    explore_proposal[1, 1] <- 0.02 # log_trans_prob
    explore_proposal[2, 2] <- 0.08 # log_obs_size
    explore_proposal[3, 3] <- 0.08 # log_E0
    explore_proposal[4, 4] <- 0.10 # logit_frac_neff

    if (verbose) cat("\n=== EXPLORATORY FIT ===\n")

    fit_explore <- tryCatch(
        {
            chlaa_fit_pmcmc(
                data = fit_data,
                pars = fit_starts[[1]],
                chain_pars = fit_starts,
                n_chains = length(fit_starts),
                n_particles = n_explore,
                dt = 0.25,
                n_steps = n_explore_steps,
                seed = seed_explore,
                prior = fit_prior,
                packer = fit_packer_stage1,
                proposal_var = explore_proposal,
                obs_interval = 7,
                time_start = time_start
            )
        },
        error = function(e) {
            if (verbose) cat("Exploratory fit failed:", conditionMessage(e), "\n")
            if (verbose) cat("Attempting with more particles...\n")
            chlaa_fit_pmcmc(
                data = fit_data,
                pars = fit_starts[[1]],
                chain_pars = fit_starts,
                n_chains = length(fit_starts),
                n_particles = max(n_explore, 200),
                dt = 0.25,
                n_steps = n_explore_steps,
                seed = seed_explore,
                prior = fit_prior,
                packer = fit_packer_stage1,
                proposal_var = explore_proposal,
                obs_interval = 7,
                time_start = time_start
            )
        }
    )

    # ---- 13. Exploratory diagnostics ----

    report_explore <- chlaa_fit_report(fit_explore, burnin = 0.25, thin = 2)
    if (verbose) cat("Exploratory acceptance rate:", report_explore$acceptance_rate, "\n")
    if (verbose) print(report_explore$posterior_summary)

    # ---- 14. Learn covariance and prepare production ----

    packer <- attr(fit_explore, "packer")
    d <- length(packer$names())
    n_samples <- dim(fit_explore$pars)[2]
    start_idx <- floor(0.25 * n_samples) + 1

    n_ch <- dim(fit_explore$pars)[3]
    pooled <- do.call(rbind, lapply(seq_len(n_ch), function(k) {
        t(fit_explore$pars[, start_idx:n_samples, k])
    }))
    colnames(pooled) <- packer$names()

    warm_vec <- apply(pooled, 2, median)
    pars_warm <- fit_starts[[1]]
    for (nm in names(warm_vec)) pars_warm[[nm]] <- warm_vec[[nm]]

    fit_starts_stage2 <- chain_median_starts(fit_explore, template_pars = fit_starts[[1]])
    fit_packer_stage2 <- make_packer(fit_starts_stage2[[1]])

    prod_proposal <- tryCatch(
        {
            cov_mat <- cov(pooled) * (2.38^2 / d)
            eig <- eigen(cov_mat, symmetric = TRUE)
            if (any(eig$values < 1e-12) || any(!is.finite(eig$values))) {
                if (verbose) cat("Warning: Degenerate covariance matrix. Using diagonal fallback.\n")
                stop("Degenerate covariance")
            }
            eig$values <- pmax(eig$values, 1e-12)
            eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
        },
        error = function(e) {
            if (verbose) cat("Using diagonal proposal (covariance learning failed).\n")
            diag(pmax(apply(pooled, 2, var), 1e-6)) * (2.38^2 / d)
        }
    )
    diag(prod_proposal) <- pmax(diag(prod_proposal), 1e-6)

    if (verbose) {
        cat("\nWarm-start parameter values:\n")
        for (nm in packer$names()) {
            cat(sprintf("  %s = %.6f\n", nm, pars_warm[[nm]]))
        }
    }

    rm(fit_explore, report_explore, packer, pooled)
    gc()

    # ---- 15. Production fit ----

    if (verbose) cat("\n=== PRODUCTION FIT ===\n")

    fit <- tryCatch(
        {
            chlaa_fit_pmcmc(
                data = fit_data,
                pars = pars_warm,
                chain_pars = fit_starts_stage2,
                n_chains = length(fit_starts_stage2),
                n_particles = n_prod,
                dt = 0.25,
                n_steps = n_prod_steps,
                seed = seed_prod,
                prior = fit_prior,
                packer = fit_packer_stage2,
                proposal_var = prod_proposal,
                obs_interval = 7,
                time_start = time_start
            )
        },
        error = function(e) {
            if (verbose) cat("Production fit failed:", conditionMessage(e), "\n")
            if (verbose) cat("Attempting with fewer steps...\n")
            chlaa_fit_pmcmc(
                data = fit_data,
                pars = pars_warm,
                chain_pars = fit_starts_stage2,
                n_chains = length(fit_starts_stage2),
                n_particles = n_prod,
                dt = 0.25,
                n_steps = max(2000, n_prod_steps %/% 2),
                seed = seed_prod,
                prior = fit_prior,
                packer = fit_packer_stage2,
                proposal_var = prod_proposal,
                obs_interval = 7,
                time_start = time_start
            )
        }
    )

    # ---- 15b. Particle-count / log-likelihood variance check ----
    # pMCMC targets the exact posterior only if Var[log p-hat(y|theta)] at the
    # posterior mode is small (rule of thumb: approx 1-2). n_prod is fixed
    # across all HZs, so check it here and bump particles + refit if needed.

    packer_prod <- attr(fit, "packer")
    n_samples_prod <- dim(fit$pars)[2]
    start_idx_prod <- floor(0.25 * n_samples_prod) + 1
    n_ch_prod <- dim(fit$pars)[3]
    pooled_prod <- do.call(rbind, lapply(seq_len(n_ch_prod), function(k) {
        t(fit$pars[, start_idx_prod:n_samples_prod, k])
    }))
    colnames(pooled_prod) <- packer_prod$names()
    theta_median_prod <- apply(pooled_prod, 2, median)
    pars_at_median <- packer_prod$unpack(theta_median_prod)

    particle_grid <- unique(pmax(1, round(n_prod * c(0.5, 1, 2, 4))))

    variance_by_particles <- vapply(particle_grid, function(np) {
        ll <- vapply(seq_len(variance_check_reps), function(r) {
            chlaa_loglik_at(
                data = fit_data, pars = pars_at_median,
                n_particles = np, seed = 10000 * np + r, dt = 0.25,
                obs_interval = 7, time_start = time_start
            )
        }, numeric(1))
        var(ll)
    }, numeric(1))
    names(variance_by_particles) <- as.character(particle_grid)

    if (verbose) {
        cat("\n=== Var[log-likelihood] vs particle count (at posterior median) ===\n")
        print(data.frame(n_particles = particle_grid, var_loglik = variance_by_particles))
    }

    var_at_n_prod <- unname(variance_by_particles[as.character(n_prod)])
    n_particles_used <- n_prod

    if (!is.na(var_at_n_prod) && var_at_n_prod > variance_target) {
        meets_target <- particle_grid[variance_by_particles <= variance_target]
        if (length(meets_target) > 0) {
            n_particles_used <- min(meets_target)
            if (verbose) {
                cat(sprintf(
                    "WARNING: Var[log-lik] at n_prod=%d = %.2f > target %.2f. Refitting production stage at n_particles=%d.\n",
                    n_prod, var_at_n_prod, variance_target, n_particles_used
                ))
            }
        } else {
            n_particles_used <- max(particle_grid)
            if (verbose) {
                cat(sprintf(
                    "WARNING: Var[log-lik] exceeds target %.2f at ALL tested particle counts. Using largest tested (%d) and proceeding with a warning.\n",
                    variance_target, n_particles_used
                ))
            }
        }

        fit <- chlaa_fit_pmcmc(
            data = fit_data,
            pars = pars_warm,
            chain_pars = fit_starts_stage2,
            n_chains = length(fit_starts_stage2),
            n_particles = n_particles_used,
            dt = 0.25,
            n_steps = n_prod_steps,
            seed = seed_prod,
            prior = fit_prior,
            packer = fit_packer_stage2,
            proposal_var = prod_proposal,
            obs_interval = 7,
            time_start = time_start
        )
    }

    variance_check_table <- tibble::tibble(
        hz = hz_name,
        n_particles = particle_grid,
        var_loglik = as.numeric(variance_by_particles)
    )

    # ---- 16. Production diagnostics ----

    report_prod <- chlaa_fit_report(fit, burnin = 0.25, thin = 2)
    if (verbose) cat("Production acceptance rate:", report_prod$acceptance_rate, "\n")
    if (verbose) print(report_prod$posterior_summary)

    if (report_prod$acceptance_rate < 0.10) {
        if (verbose) cat("WARNING: Low acceptance rate may indicate identifiability issues.\n")
    }

    # R-hat and ESS (proper multi-chain diagnostic)
    pars_arr <- fit$pars # (param, iter, chain)
    dimnames(pars_arr) <- list(attr(fit, "packer")$names(), NULL, NULL)
    draws_diag <- posterior::as_draws_array(aperm(pars_arr, c(2, 3, 1))) # (iter, chain, param)
    rhat_ess <- posterior::summarise_draws(draws_diag, "rhat", "ess_bulk", "ess_tail")
    if (verbose) {
        cat("\n=== R-hat and ESS ===\n")
        print(rhat_ess)
    }
    rm(pars_arr, draws_diag)

    # ---- 17. Diagnostic plots ----

    p_trace <- chlaa_plot_trace(fit, parameters = natural_fit_names, burnin = 0.25, scale = "natural")
    ggsave(file.path(fig_dir, sprintf("fitting_%s_production_trace.png", hz_name)),
        p_trace,
        width = 12, height = 8, dpi = 300
    )

    p_lltrace <- chlaa_plot_likelihood_trace(fit, burnin = 0.25, thin = 2)
    ggsave(file.path(fig_dir, sprintf("fitting_%s_production_likelihood_trace.png", hz_name)),
        p_lltrace,
        width = 12, height = 8, dpi = 300
    )

    p_pairs <- chlaa_plot_parameter_pairs(
        fit,
        parameters = natural_fit_names,
        burnin = 0.25, scale = "natural", max_points = 3000
    )
    ggsave(file.path(fig_dir, sprintf("fitting_%s_production_pairs.png", hz_name)),
        plot = p_pairs, width = 10, height = 10, dpi = 300
    )

    p_dist <- chlaa_plot_parameter_distributions(
        fit,
        parameters = natural_fit_names, burnin = 0.25, scale = "natural"
    )
    ggsave(file.path(fig_dir, sprintf("fitting_%s_production_distributions.png", hz_name)),
        p_dist,
        width = 12, height = 8, dpi = 300
    )

    # ---- 18. Fit plot ----

    n_weeks <- nrow(hz_data_weekly)

    p_fit <- plot_case_fit(
        fit, hz_data_weekly,
        sprintf("%s weekly reported cases", tools::toTitleCase(hz_name)),
        seed = seed_prod
    )
    ggsave(file.path(fig_dir, sprintf("fitting_%s_production_fit.png", hz_name)),
        p_fit,
        width = 10, height = 6, dpi = 300
    )

    # ---- 19. Compute R0 and diagnostics for this HZ ----

    tr_s <- chlaa_fit_trace(fit, burnin = 0.25, scale = "sampled")
    tp_draws <- exp(tr_s |> filter(parameter == "log_trans_prob") |> pull(value))
    fn_draws <- plogis(tr_s |> filter(parameter == "logit_frac_neff") |> pull(value))
    pars_for_r0 <- pars_warm
    pars_for_r0$trans_prob <- tp_draws
    pars_for_r0$N <- fn_draws * pop_hz
    R0_draws <- chlaa_r0(pars_for_r0)

    r0_table <- tibble::tibble(
        hz = hz_name, pop = pop_hz,
        tp_med = median(tp_draws),
        tp_lo = quantile(tp_draws, 0.025),
        tp_hi = quantile(tp_draws, 0.975),
        frac_med = median(fn_draws),
        frac_lo = quantile(fn_draws, 0.025),
        frac_hi = quantile(fn_draws, 0.975),
        R0_med = median(R0_draws),
        R0_lo = quantile(R0_draws, 0.025),
        R0_hi = quantile(R0_draws, 0.975),
        N_eff_med = median(fn_draws) * pop_hz
    )

    # Build diagnostics row
    E0_draws <- exp(tr_s |> filter(parameter == "log_E0") |> pull(value))
    obs_draws <- exp(tr_s |> filter(parameter == "log_obs_size") |> pull(value))

    diag_row <- tibble::tibble(
        hz = hz_name,
        acceptance_rate = report_prod$acceptance_rate,
        trans_prob_med = median(tp_draws),
        obs_size_med = median(obs_draws),
        E0_med = median(E0_draws),
        frac_neff_med = median(fn_draws),
        N_eff_med = median(fn_draws) * pop_hz,
        R0_med = median(R0_draws),
        n_particles_used = n_particles_used,
        loglik_var_at_n_prod = as.numeric(var_at_n_prod),
        loglik_var_target = variance_target
    )
    # Append per-parameter R-hat and ESS
    for (i in seq_len(nrow(rhat_ess))) {
        nm <- rhat_ess$variable[i]
        diag_row[[paste0("rhat_", nm)]] <- rhat_ess$rhat[i]
        diag_row[[paste0("ess_bulk_", nm)]] <- rhat_ess$ess_bulk[i]
        diag_row[[paste0("ess_tail_", nm)]] <- rhat_ess$ess_tail[i]
    }

    # --- Budget table: posterior-median deterministic run ---
    p_med <- make_start(median(tp_draws), median(obs_draws), median(E0_draws),
        frac_neff = median(fn_draws)
    )
    s_med <- chlaa_simulate(
        pars = p_med, time = hz_data_weekly$time,
        n_particles = 1, dt = 0.25, deterministic = TRUE
    )
    sat <- s_med$C / (s_med$C + p_med$contam_half_sat)
    N_eff <- median(fn_draws) * pop_hz

    budget_row <- tibble::tibble(
        hz = hz_name,
        obs_total = sum(hz_data_weekly$cases),
        model_total = sum(s_med$inc_symptoms_weekly) * RR_FIXED,
        total_ratio = model_total / obs_total,
        obs_peak = max(hz_data_weekly$cases),
        model_peak = max(s_med$inc_symptoms_weekly) * RR_FIXED,
        peak_ratio = model_peak / obs_peak,
        peak_wk_obs = which.max(hz_data_weekly$cases),
        peak_wk_mod = which.max(s_med$inc_symptoms_weekly),
        peak_wk_err = peak_wk_mod - peak_wk_obs,
        inf_needed = obs_total / RR_FIXED / 0.25,
        N_eff = N_eff,
        rounds_needed = inf_needed / N_eff,
        sat_min = min(sat), sat_max = max(sat)
    )

    # --- PPC coverage ---
    fc_ppc <- chlaa_forecast_from_fit(fit,
        time = hz_data_weekly$time,
        vars = "inc_symptoms_weekly", include_cases = TRUE,
        obs_model = "nbinom", n_draws = 200,
        burnin = 0.25, seed = seed_prod, dt = 0.25
    )
    cv <- fc_ppc |>
        filter(variable == "cases") |>
        left_join(hz_data_weekly |> select(time, cases), by = "time")
    diag_row$cover_50 <- mean(cv$cases >= cv$q0p25 & cv$cases <= cv$q0p75)
    diag_row$cover_95 <- mean(cv$cases >= cv$q0p025 & cv$cases <= cv$q0p975)

    # Save per-HZ tables (for array job mode)
    write.csv(r0_table, file.path(tables_dir, sprintf("%s_r0_table.csv", hz_name)), row.names = FALSE)
    write.csv(diag_row, file.path(tables_dir, sprintf("%s_diagnostics.csv", hz_name)), row.names = FALSE)
    write.csv(budget_row, file.path(tables_dir, sprintf("%s_budget.csv", hz_name)), row.names = FALSE)
    write.csv(variance_check_table, file.path(tables_dir, sprintf("%s_particle_variance.csv", hz_name)), row.names = FALSE)

    # ---- 20. Save fit object ----

    fit_output <- list(
        hz_name = hz_name,
        pop_hz = pop_hz,
        fit = fit,
        pars_start = fit_starts_stage2,
        pars_warm = pars_warm,
        report = report_prod,
        rhat_ess = rhat_ess,
        r0_table = r0_table,
        diagnostics = diag_row,
        budget = budget_row,
        variance_check_table = variance_check_table,
        observed = hz_data_weekly,
        outbreak_start = outbreak_start,
        outbreak_end = outbreak_end,
        total_cases = total_cases,
        n_weeks = n_weeks,
        fitted_parameters = natural_fit_names,
        vax1_schedule = vax1_schedule,
        vax1_arrays = vax1_arrays,
        vax1_total_doses = vax1_total_doses_hz,
        vax2_schedule = vax2_schedule,
        vax2_arrays = vax2_arrays,
        vax2_total_doses = vax2_total_doses_hz,
        timestamp = Sys.time()
    )

    saveRDS(fit_output, file.path(rds_dir, sprintf("%s_fit.rds", hz_name)))

    if (verbose) cat("\nFit saved to:", file.path(rds_dir, sprintf("%s_fit.rds", hz_name)), "\n")
    if (verbose) cat("Figures saved to:", fig_dir, "\n")

    return(fit_output)
}

# ---- Main Execution ----

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
            ), file.path(rds_dir, sprintf("%s_FAILED.rds", hz_to_fit)))
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

    # Combine per-HZ tables into aggregated files
    successful <- Filter(Negate(is.null), results)
    if (length(successful) > 0) {
        r0_all <- dplyr::bind_rows(lapply(successful, function(x) x$r0_table))
        diag_all <- dplyr::bind_rows(lapply(successful, function(x) x$diagnostics))
        budget_all <- dplyr::bind_rows(lapply(successful, function(x) x$budget))
        write.csv(r0_all, file.path(tables_dir, "all_hz_r0_table.csv"), row.names = FALSE)
        write.csv(diag_all, file.path(tables_dir, "all_hz_diagnostics.csv"), row.names = FALSE)
        write.csv(budget_all, file.path(tables_dir, "all_hz_budget.csv"), row.names = FALSE)
        cat("\nR0 table, diagnostics, and budget saved to:", tables_dir, "\n")
        print(r0_all, width = Inf)
    }
}

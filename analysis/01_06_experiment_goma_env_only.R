# =============================================================================
# Goma single-HZ fitting (matches 01_02 workflow for local troubleshooting)
# =============================================================================

library(chlaa)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

fig_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/figures_experiment_goma_env"
data_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/analysis/data"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

hz_name <- "goma"

# ---- 1. Load data ----

hz_params_long <- read.csv(file.path(data_dir, "hz_parameters.csv"), stringsAsFactors = FALSE)
idsr <- read.csv(file.path(data_dir, "IDSR_dataset.csv"))
vax_dates <- read.csv(file.path(data_dir, "ocv_vaccination_dates.csv"), stringsAsFactors = FALSE)

hz_rows_long <- hz_params_long %>% filter(hz == hz_name)

outbreak_start <- hz_rows_long %>%
    filter(parameter == "outbreak_start") %>% pull(value) %>% first() %>% as.Date()
outbreak_end <- hz_rows_long %>%
    filter(parameter == "outbreak_end") %>% pull(value) %>% first() %>% as.Date()

cat("Outbreak window:", as.character(outbreak_start), "to", as.character(outbreak_end), "\n")

# Burn-in: model starts this many days before first observation (day 7)
time_start <- -21L

# ---- 2. Helpers ----

safe_date_to_day <- function(date_str, origin) {
    d <- as.Date(date_str, format = "%Y-%m-%d")
    if (is.na(d)) return(0L)
    as.integer(d - origin)
}

get_param <- function(param_name) {
    val <- hz_rows_long %>% filter(parameter == param_name) %>% pull(value) %>% first()
    if (is.null(val) || length(val) == 0) return(NA)
    val
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

plot_case_fit <- function(fit, observed, title, seed, burnin = 0.25) {
    fc <- chlaa_forecast_from_fit(
        fit = fit, time = observed$time, vars = "inc_symptoms_weekly",
        include_cases = TRUE, obs_model = "nbinom",
        n_draws = 50, burnin = burnin, seed = seed, dt = 1, deterministic = FALSE
    )
    fit_cases <- fc %>%
        filter(variable == "cases") %>%
        left_join(observed %>% select(time, date), by = "time")
    ggplot() +
        geom_ribbon(data = fit_cases, aes(date, ymin = q0p025, ymax = q0p975),
            fill = "#6baed6", alpha = 0.25) +
        geom_ribbon(data = fit_cases, aes(date, ymin = q0p25, ymax = q0p75),
            fill = "#6baed6", alpha = 0.45) +
        geom_line(data = fit_cases, aes(date, q0p5), colour = "#08519c", linewidth = 0.8) +
        geom_point(data = observed, aes(date, cases), size = 1.6) +
        labs(x = NULL, y = "Weekly reported cases", title = title)
}

# ---- 3. Intervention dates ----

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

# ---- 4. Vaccination schedules ----

generate_vax_schedule <- function(total_doses, start_date, end_date, outbreak_start) {
    profile <- c(0.305, 0.377, 0.227, 0.074, 0.014, 0.003)
    n_days <- max(as.integer(end_date - start_date) + 1L, 1L)
    day_weights <- approx(seq(0, 1, length.out = length(profile)),
        profile, xout = seq(0, 1, length.out = n_days))$y
    day_weights <- day_weights / sum(day_weights)
    daily_doses <- round(total_doses * day_weights)
    daily_doses[which.max(daily_doses)] <- daily_doses[which.max(daily_doses)] +
        (total_doses - sum(daily_doses))
    start_day <- as.integer(start_date - outbreak_start)
    sched <- data.frame(time = seq(start_day, by = 1L, length.out = n_days),
        doses = as.numeric(daily_doses))
    # Ensure schedule covers time_start for interpolation
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

# Vax1
vax_hz <- vax_dates %>%
    filter(healthzone == hz_name | healthzone == gsub("_", " ", hz_name) |
        healthzone == gsub("_", "-", hz_name))

vax1_total_doses_hz <- hz_rows_long %>%
    filter(parameter == "vax1_total_doses") %>% pull(value) %>% as.numeric() %>% first()

vax1_campaigns <- vax_hz %>%
    filter(parameter %in% c("vax1_start", "vax1_end", "vax1_total_doses")) %>%
    distinct(dose_round, parameter, .keep_all = TRUE) %>%
    pivot_wider(names_from = parameter, values_from = value, id_cols = c(healthzone, dose_round))

if (!("vax1_start" %in% names(vax1_campaigns))) {
    vax1_campaigns <- vax1_campaigns[0, ]
} else {
    vax1_campaigns <- vax1_campaigns %>% filter(!is.na(vax1_start) & vax1_start != "NA")
}

if (nrow(vax1_campaigns) > 0 && !is.na(vax1_total_doses_hz) && vax1_total_doses_hz > 0) {
    vax1_camp <- vax1_campaigns[1, ]
    vax1_start_date <- as.Date(vax1_camp$vax1_start)
    vax1_end_date <- as.Date(vax1_camp$vax1_end)
    if (!is.na(vax1_start_date) && !is.na(vax1_end_date) &&
        vax1_start_date >= outbreak_start && vax1_start_date <= outbreak_end) {
        vax1_schedule <- generate_vax_schedule(vax1_total_doses_hz, vax1_start_date, vax1_end_date, outbreak_start)
        vax1_arrays <- prepare_vax_arrays(vax1_schedule)
        vax1_start_day <- min(vax1_schedule$time)
        vax1_end_day <- max(vax1_schedule$time) + 1
    } else {
        vax1_schedule <- generate_empty_vax_schedule()
        vax1_arrays <- prepare_vax_arrays(vax1_schedule)
        vax1_start_day <- 0L; vax1_end_day <- 0L
    }
} else {
    vax1_schedule <- generate_empty_vax_schedule()
    vax1_arrays <- prepare_vax_arrays(vax1_schedule)
    vax1_start_day <- 0L; vax1_end_day <- 0L
}

# Vax2
vax2_total_doses_hz <- hz_rows_long %>%
    filter(parameter == "vax2_total_doses") %>% pull(value) %>% as.numeric() %>% first()

vax2_campaigns <- vax_hz %>%
    filter(parameter %in% c("vax2_start", "vax2_end", "vax2_total_doses")) %>%
    distinct(dose_round, parameter, .keep_all = TRUE) %>%
    pivot_wider(names_from = parameter, values_from = value, id_cols = c(healthzone, dose_round))

if (!("vax2_start" %in% names(vax2_campaigns))) {
    vax2_campaigns <- vax2_campaigns[0, ]
} else {
    vax2_campaigns <- vax2_campaigns %>% filter(!is.na(vax2_start) & vax2_start != "NA")
}

if (nrow(vax2_campaigns) > 0 && !is.na(vax2_total_doses_hz) && vax2_total_doses_hz > 0) {
    vax2_camp <- vax2_campaigns[1, ]
    vax2_start_date <- as.Date(vax2_camp$vax2_start)
    vax2_end_date <- as.Date(vax2_camp$vax2_end)
    if (!is.na(vax2_start_date) && !is.na(vax2_end_date) &&
        vax2_start_date >= outbreak_start && vax2_start_date <= outbreak_end) {
        vax2_schedule <- generate_vax_schedule(vax2_total_doses_hz, vax2_start_date, vax2_end_date, outbreak_start)
        vax2_arrays <- prepare_vax_arrays(vax2_schedule)
        vax2_start_day <- min(vax2_schedule$time)
        vax2_end_day <- max(vax2_schedule$time) + 1
    } else {
        vax2_schedule <- generate_empty_vax_schedule()
        vax2_arrays <- prepare_vax_arrays(vax2_schedule)
        vax2_start_day <- 0L; vax2_end_day <- 0L
    }
} else {
    vax2_schedule <- generate_empty_vax_schedule()
    vax2_arrays <- prepare_vax_arrays(vax2_schedule)
    vax2_start_day <- 0L; vax2_end_day <- 0L
}

# ---- 5. Weekly case data ----

hz_weekly <- idsr %>%
    filter(hz == hz_name) %>%
    mutate(date = as.Date(date), deaths = replace_na(deaths, 0L)) %>%
    select(date, year, week, cases, deaths, population)

hz_outbreak <- hz_weekly %>%
    filter(date >= outbreak_start, date <= outbreak_end) %>%
    arrange(date)

pop_hz <- hz_weekly$population[1]
total_cases <- sum(hz_outbreak$cases)
cat("Population:", pop_hz, "\nOutbreak weeks:", nrow(hz_outbreak), "\nTotal cases:", total_cases, "\n")

hz_data_weekly <- hz_outbreak %>%
    mutate(time = seq_len(n()) * 7L) %>%
    select(time, date, cases, deaths) %>%
    arrange(time)

# ---- 6. E0 initialization ----

expected_reporting_rate <- 0.10
seed_date <- outbreak_start - 14
seed_row <- hz_weekly %>% filter(date <= seed_date) %>% arrange(desc(date)) %>% slice(1)

if (total_cases < 50) {
    E0_val <- max(5, ceiling(mean(hz_outbreak$cases[1:min(3, nrow(hz_outbreak))]) / expected_reporting_rate))
} else if (nrow(seed_row) > 0 && seed_row$cases > 0) {
    E0_val <- ceiling(seed_row$cases / expected_reporting_rate)
} else {
    E0_val <- ceiling(max(5, hz_outbreak$cases[1]) / expected_reporting_rate)
}
E0_val <- min(E0_val, 800)
E0_val <- max(10, E0_val)
cat("E0:", E0_val, "\n")

# ---- 7. Fitting setup (matches 01_02) ----

natural_fit_names <- c("trans_prob", "reporting_rate", "obs_size", "E0")
fit_names <- c("log_trans_prob", "logit_reporting_rate", "log_obs_size", "log_E0")
freeze_reporting_rate <- TRUE

add_transformed_values <- function(pars, freeze_reporting_rate = FALSE) {
    pars$log_trans_prob <- log(pars$trans_prob)
    if (!freeze_reporting_rate) pars$logit_reporting_rate <- qlogis(pars$reporting_rate)
    pars$log_obs_size <- log(pars$obs_size)
    pars$log_E0 <- log(pars$E0)
    pars
}

make_packer <- function(pars, freeze_reporting_rate = FALSE) {
    fit_names_use <- if (freeze_reporting_rate) {
        c("log_trans_prob", "log_obs_size", "log_E0")
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
                E0 = exp(p$log_E0)
            )
            if (!freeze_reporting_rate) out$reporting_rate <- plogis(p$logit_reporting_rate)
            out
        }
    )
}

make_start <- function(trans_prob, reporting_rate, obs_size, E0,
                       freeze_reporting_rate = FALSE) {
    pars_args <- list(
        N = pop_hz, Sev0 = 0, E0 = E0, M0 = 0,
        immunity_asym = 280, beta_p2p = 0, contam_half_sat = 1.0,
        trans_prob = trans_prob, incubation_time = 4.845, duration_sym = 14.48,
        seek_mild = 0.1, seek_severe = 0.85,
        reporting_rate = reporting_rate, fatality_treated = 0.0021, fatality_untreated = 0.5,
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
    out$vax1_schedule_time <- vax1_arrays$time
    out$vax1_schedule_doses <- vax1_arrays$doses
    out$n_vax1_schedule <- length(vax1_arrays$time)
    out$vax2_schedule_time <- vax2_arrays$time
    out$vax2_schedule_doses <- vax2_arrays$doses
    out$n_vax2_schedule <- length(vax2_arrays$time)
    add_transformed_values(out, freeze_reporting_rate = freeze_reporting_rate)
}

# ---- 8. Priors and starting points ----

fit_prior_full <- monty::monty_dsl({
    log_trans_prob ~ Normal(-7.600902, 1.2)
    logit_reporting_rate ~ Uniform(-2.944439, -1.098612)
    log_obs_size ~ Uniform(0, 5.703782)
    log_E0 ~ Uniform(2.302585, 6.684612)
}, gradient = FALSE)

fit_prior_freeze <- monty::monty_dsl({
    log_trans_prob ~ Normal(-7.600902, 1.2)
    log_obs_size ~ Uniform(0, 5.703782)
    log_E0 ~ Uniform(2.302585, 6.684612)
}, gradient = FALSE)

# trans_prob scaled by population
ref_pop <- 300000
tp_scale <- ref_pop / pop_hz
tp_starts <- c(0.003225, 1.6e-3, 4.0e-4) * tp_scale

fit_starts_full <- list(
    make_start(tp_starts[1], 0.15, 30, E0_val),
    make_start(tp_starts[2], 0.08, 20, max(10, round(E0_val * 0.5))),
    make_start(tp_starts[3], 0.20, 100, min(800, max(10, round(E0_val * 1.5))))
)

fit_starts_freeze <- list(
    make_start(tp_starts[1], 0.15, 30, E0_val, freeze_reporting_rate = TRUE),
    make_start(tp_starts[2], 0.15, 20, max(10, round(E0_val * 0.5)), freeze_reporting_rate = TRUE),
    make_start(tp_starts[3], 0.15, 100, min(800, max(10, round(E0_val * 1.5))), freeze_reporting_rate = TRUE)
)

fit_prior_stage1 <- if (freeze_reporting_rate) fit_prior_freeze else fit_prior_full
fit_starts_stage1 <- if (freeze_reporting_rate) fit_starts_freeze else fit_starts_full
fit_packer_stage1 <- make_packer(fit_starts_stage1[[1]], freeze_reporting_rate = freeze_reporting_rate)

fit_data <- data.frame(time = hz_data_weekly$time, cases = hz_data_weekly$cases,
    deaths = hz_data_weekly$deaths)

# ---- 9. Exploratory fit (freeze reporting_rate) ----

n_params <- if (freeze_reporting_rate) 3 else 4
explore_proposal <- matrix(0, n_params, n_params)
if (freeze_reporting_rate) {
    explore_proposal[1, 1] <- 0.02  # log_trans_prob
    explore_proposal[2, 2] <- 0.08  # log_obs_size
    explore_proposal[3, 3] <- 0.08  # log_E0
} else {
    explore_proposal[1, 1] <- 0.02  # log_trans_prob
    explore_proposal[2, 2] <- 0.05  # logit_reporting_rate
    explore_proposal[3, 3] <- 0.08  # log_obs_size
    explore_proposal[4, 4] <- 0.08  # log_E0
}

cat("\n=== EXPLORATORY FIT ===\n")
fit_explore <- chlaa_fit_pmcmc(
    data = fit_data,
    pars = fit_starts_stage1[[1]],
    chain_pars = fit_starts_stage1,
    n_chains = length(fit_starts_stage1),
    n_particles = 30,
    n_steps = 300,
    seed = 42,
    prior = fit_prior_stage1,
    packer = fit_packer_stage1,
    proposal_var = explore_proposal,
    obs_interval = 7,
    time_start = time_start
)

report_explore <- chlaa_fit_report(fit_explore, burnin = 0.25, thin = 2)
cat("Exploratory acceptance:", report_explore$acceptance_rate, "\n")
print(report_explore$posterior_summary)

# ---- 10. Learn covariance and prepare production ----

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

# Unfreeze reporting_rate for production
fit_starts_stage2 <- if (freeze_reporting_rate) {
    chain_median_starts(fit_explore, template_pars = fit_starts_full[[1]])
} else {
    fit_starts_full
}

fit_prior_stage2 <- fit_prior_full
fit_packer_stage2 <- make_packer(fit_starts_stage2[[1]], freeze_reporting_rate = FALSE)

prod_proposal <- tryCatch({
    cov_mat <- cov(pooled) * (2.38^2 / d)
    eig <- eigen(cov_mat, symmetric = TRUE)
    if (any(eig$values < 1e-12) || any(!is.finite(eig$values))) stop("Degenerate covariance")
    eig$values <- pmax(eig$values, 1e-12)
    eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
}, error = function(e) {
    cat("Using diagonal proposal (covariance learning failed).\n")
    diag(pmax(apply(pooled, 2, var), 1e-6)) * (2.38^2 / d)
})

# Expand proposal if freeze -> full
if (freeze_reporting_rate && ncol(prod_proposal) == 3) {
    prop4 <- matrix(0, 4, 4)
    idx_map <- c(1L, 3L, 4L)
    prop4[idx_map, idx_map] <- prod_proposal
    prop4[2, 2] <- 0.05 * (2.38^2 / 4)
    prod_proposal <- prop4
}
diag(prod_proposal) <- pmax(diag(prod_proposal), 1e-6)

rm(fit_explore, report_explore, packer, pooled, pars_mat)
gc()

# ---- 11. Production fit (all 4 params) ----

cat("\n=== PRODUCTION FIT ===\n")
fit <- chlaa_fit_pmcmc(
    data = fit_data,
    pars = pars_warm,
    chain_pars = fit_starts_stage2,
    n_chains = length(fit_starts_stage2),
    n_particles = 50,
    n_steps = 1500,
    seed = 123,
    prior = fit_prior_stage2,
    packer = fit_packer_stage2,
    proposal_var = prod_proposal,
    obs_interval = 7,
    time_start = time_start
)

report_prod <- chlaa_fit_report(fit, burnin = 0.25, thin = 2)
cat("Production acceptance:", report_prod$acceptance_rate, "\n")
print(report_prod$posterior_summary)

# ---- 12. Diagnostics ----

p_trace <- chlaa_plot_trace(fit, parameters = natural_fit_names, burnin = 0.25, scale = "natural")
ggsave(file.path(fig_dir, "goma_production_trace.png"), p_trace, width = 12, height = 8, dpi = 300)

p_lltrace <- chlaa_plot_likelihood_trace(fit, burnin = 0.25, thin = 2)
ggsave(file.path(fig_dir, "goma_production_likelihood_trace.png"), p_lltrace, width = 12, height = 8, dpi = 300)

p_pairs <- chlaa_plot_parameter_pairs(fit, parameters = natural_fit_names,
    burnin = 0.25, scale = "natural", max_points = 3000)
ggsave(file.path(fig_dir, "goma_production_pairs.png"), p_pairs, width = 10, height = 10, dpi = 300)

p_dist <- chlaa_plot_parameter_distributions(fit, parameters = natural_fit_names,
    burnin = 0.25, scale = "natural")
ggsave(file.path(fig_dir, "goma_production_distributions.png"), p_dist, width = 12, height = 8, dpi = 300)

p_fit <- plot_case_fit(fit, hz_data_weekly,
    sprintf("%s weekly reported cases", tools::toTitleCase(hz_name)), seed = 123)
ggsave(file.path(fig_dir, "goma_production_fit.png"), p_fit, width = 10, height = 6, dpi = 300)

# ---- 13. Save ----

saveRDS(
    list(fit = fit, observed = hz_data_weekly, report = report_prod,
         parameter_summary = chlaa_posterior_summary(fit, burnin = 0.25)),
    file.path(fig_dir, "goma_fit.rds")
)

cat("\nDone. Figures saved to:", fig_dir, "\n")

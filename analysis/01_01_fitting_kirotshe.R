# =============================================================================
# Kirotshe Fitting Vignette
# Reproduces: https://ojwatson.github.io/chlaa/articles/fitting.html
#
# This script fits the chlaa cholera model to Kirotshe health zone data using
# pMCMC. It progresses through: (1) synthetic data validation on raw and
# transformed scales, (2) pilot runs to learn proposal covariance, (3)
# deterministic main runs, and (4) stochastic particle-filter runs on both
# synthetic and real data.
# =============================================================================

library(chlaa)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(coda)

# Directory for saving figures
fig_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/figures"

# -----------------------------------------------------------------------------
# 1. File path helpers
# -----------------------------------------------------------------------------
# These utilities manage file paths for locating data files within the
# repository structure.

repo_file <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...))
  out <- candidates[file.exists(candidates)][1]
  if (is.na(out)) stop("Could not find file: ", file.path(...), call. = FALSE)
  out
}

repo_output_file <- function(...) {
  roots <- c(".", "..")
  root <- roots[file.exists(file.path(roots, "DESCRIPTION"))][1]
  if (is.na(root)) stop("Could not find repository root", call. = FALSE)
  file.path(root, ...)
}

extdata_file <- function(...) {
  out <- system.file("extdata", ..., package = "chlaa")
  if (nzchar(out)) {
    return(out)
  }
  repo_file("inst", "extdata", ...)
}

# -----------------------------------------------------------------------------
# 2. Load Kirotshe data
# -----------------------------------------------------------------------------
# Imports weekly case counts and intervention metadata for the Kirotshe health
# zone, establishing temporal bounds for the outbreak window.

hz_name <- "kirotshe"

data_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/analysis/data"

kirotshe_meta <- read_csv(
  file.path(data_dir, "hz_parameters.csv"),
  col_types = cols(.default = col_character()),
  na = c("NA", "")
) |>
  filter(hz == hz_name) |>
  select(-hz) |>
  pivot_wider(names_from = parameter, values_from = value)

outbreak_start_date <- as.Date(kirotshe_meta$outbreak_start)
outbreak_end_date <- as.Date(kirotshe_meta$outbreak_end)

kirotshe <- read_csv(file.path(data_dir, "IDSR_dataset.csv"), show_col_types = FALSE) |>
  filter(hz == hz_name) |>
  mutate(date = as.Date(date)) |>
  filter(date >= outbreak_start_date, date <= outbreak_end_date) |>
  mutate(time = seq_len(n()) * 7L) |>
  select(date, time, cases, deaths, population, cases_pop)
# NOTE: deaths column retains NAs; chlaa_prepare_data() converts them to
# deaths=0 + has_deaths=0 (non-informative), while actual death counts get
# has_deaths=1 (informative).

outbreak_start <- outbreak_start_date
outbreak_end <- outbreak_end_date

kirotshe |>
  select(date, time, cases, population, cases_pop) |>
  head()

# -----------------------------------------------------------------------------
# 3. Plot observed data
# -----------------------------------------------------------------------------
# Bar chart showing weekly reported case counts over time.

p_observed <- ggplot(kirotshe, aes(date, cases)) +
  geom_col(width = 6, fill = "grey65") +
  labs(
    x = NULL,
    y = "Weekly reported cases",
    title = "Kirotshe IDSR outbreak window"
  )
print(p_observed)

ggsave(
  file.path(fig_dir, "fitting_observed_cases.png"),
  plot = p_observed, width = 10, height = 6, dpi = 300
)

# -----------------------------------------------------------------------------
# 4. Convert dates to model time offsets
# -----------------------------------------------------------------------------
# Translates calendar dates to day offsets relative to outbreak start, handling
# various public health interventions (ORC, CTC, chlorination, hygiene, CATI,
# and latrines).

date_to_day <- function(x, origin) {
  d <- as.Date(x)
  if (length(d) == 0 || is.na(d)) 0L else as.integer(d - origin)
}

num_or <- function(x, default = 0) {
  y <- suppressWarnings(as.numeric(x))
  if (length(y) == 0 || is.na(y)) default else y
}

intervention_days <- list(
  orc_start = date_to_day(kirotshe_meta$orc_start, outbreak_start),
  orc_end = date_to_day(kirotshe_meta$orc_end, outbreak_start),
  ctc_start = date_to_day(kirotshe_meta$ctc_start, outbreak_start),
  ctc_end = date_to_day(kirotshe_meta$ctc_end, outbreak_start),
  chlor_start = date_to_day(kirotshe_meta$chlor_start, outbreak_start),
  chlor_end = date_to_day(kirotshe_meta$chlor_end, outbreak_start),
  chlor_effect = num_or(kirotshe_meta$chlor_effect),
  hyg_start = date_to_day(kirotshe_meta$hyg_start, outbreak_start),
  hyg_end = date_to_day(kirotshe_meta$hyg_end, outbreak_start),
  hyg_effect = num_or(kirotshe_meta$hyg_effect),
  cati_start = date_to_day(kirotshe_meta$cati_start, outbreak_start),
  cati_end = date_to_day(kirotshe_meta$cati_end, outbreak_start),
  cati_effect = num_or(kirotshe_meta$cati_effect),
  lat_start = date_to_day(kirotshe_meta$lat_start, outbreak_start),
  lat_end = date_to_day(kirotshe_meta$lat_end, outbreak_start),
  lat_effect = num_or(kirotshe_meta$lat_effect)
)

tibble(
  quantity = names(intervention_days),
  value = unlist(intervention_days)
)

# Burn-in: model starts this many days before first observation (day 7)
# to let environmental contamination build up before the particle filter
time_start <- -21L

# Reference population for scaling contam_half_sat (Kirotshe ~516k)
# This makes the environmental FOI depend on per-capita prevalence rather
# than absolute case counts, so trans_prob stays comparable across zones.
pop_ref <- 516000
contam_half_sat_scaled <- 0.01 * (kirotshe$population[1] / pop_ref)

# -----------------------------------------------------------------------------
# 5. Parameter factory function
# -----------------------------------------------------------------------------
# Creates parameter objects for simulation and fitting, accepting four free
# parameters while fixing others to epidemiologically informed defaults.

make_kirotshe_pars <- function(trans_prob,
                               reporting_rate,
                               obs_size,
                               E0) {
  out <- do.call(
    chlaa_parameters,
    c(
      list(
        N = kirotshe$population[1],
        contact_rate = 0,
        contam_half_sat = contam_half_sat_scaled,
        trans_prob = trans_prob,
        E0 = E0,
        Sev0 = 0,
        M0 = 0,
        C0 = 0,
        reporting_rate = reporting_rate,
        obs_size = obs_size,
        seek_severe = 0.85,
        fatality_treated = 0.0021,
        fatality_untreated = 0.5,
        death_reporting_rate = 0.5,
        obs_size_deaths = 1.0
      ),
      intervention_days
    )
  )
  # Empty vaccination schedule arrays required by the odin model;
  # must cover time_start for dust2 interpolation
  out$vax1_schedule_time <- c(time_start, time_start + 1L)
  out$vax1_schedule_doses <- c(0, 0)
  out$n_vax1_schedule <- 2L
  out$vax2_schedule_time <- c(time_start, time_start + 1L)
  out$vax2_schedule_doses <- c(0, 0)
  out$n_vax2_schedule <- 2L
  out
}

# -----------------------------------------------------------------------------
# 6. Generate synthetic data
# -----------------------------------------------------------------------------
# Simulates an outbreak with known parameter values, generating synthetic
# weekly case counts from a negative binomial observation model.

truth_pars <- make_kirotshe_pars(
  trans_prob = 8.5e-4,
  reporting_rate = 0.15,
  obs_size = 120,
  E0 = 80
)

truth <- chlaa_simulate(
  pars = truth_pars,
  time = kirotshe$time,
  n_particles = 1,
  seed = 1,
  dt = 1,
  deterministic = TRUE
)

set.seed(1)
synthetic_weekly <- kirotshe |>
  transmute(
    time,
    date,
    population,
    inc_symptoms_truth = truth$inc_symptoms_weekly,
    mu_cases = truth_pars$reporting_rate * truth$inc_symptoms_weekly,
    cases = rnbinom(n(), mu = pmax(mu_cases, 0.01), size = truth_pars$obs_size)
  )


natural_fit_names <- c("trans_prob", "reporting_rate", "obs_size", "E0")
truth_vec <- unlist(truth_pars[natural_fit_names])
truth_values <- tibble(parameter = names(truth_vec), truth = as.numeric(truth_vec))

truth_values


# -----------------------------------------------------------------------------
# 7. Plot synthetic data
# -----------------------------------------------------------------------------
# Overlays noisy case observations with the underlying expected incidence curve.

p_synthetic <- ggplot() +
  geom_col(data = synthetic_weekly, aes(date, cases), width = 6, fill = "grey70") +
  geom_line(data = synthetic_weekly, aes(date, mu_cases), colour = "#238b45", linewidth = 0.8) +
  labs(
    x = NULL,
    y = "Weekly reported cases",
    title = "Synthetic weekly data",
    subtitle = "Bars are noisy observations; the green line is the observation mean"
  )
print(p_synthetic)

ggsave(
  file.path(fig_dir, "fitting_synthetic_data.png"),
  plot = p_synthetic, width = 10, height = 6, dpi = 300
)

# -----------------------------------------------------------------------------
# 8. Diagnostic configuration
# -----------------------------------------------------------------------------
# Defines MCMC configuration: three parallel chains with specified iteration
# counts for pilot, deterministic, and particle-filter phases.

n_chains <- 3L
pilot_steps <- 5000L
deterministic_steps <- 20000L
particle_steps <- 5000L
particle_count <- 50L
# Synthetic data has no deaths observations; set deaths = 0 with has_deaths
# handled by chlaa_prepare_data (deaths = 0 + has_deaths = 0 → non-informative).
fit_data_synthetic <- synthetic_weekly |> select(time, cases)
fit_data_real <- kirotshe |> select(time, cases, deaths)

# -----------------------------------------------------------------------------
# 9. Trace extraction function
# -----------------------------------------------------------------------------
# Converts trace data to wide format for convenience in downstream analyses
# after removing burn-in period.

draws_wide <- function(fit, burnin = 0.25, scale = c("sampled", "natural")) {
  scale <- match.arg(scale)
  chlaa_fit_trace(fit, burnin = burnin, scale = scale) |>
    pivot_wider(names_from = parameter, values_from = value) |>
    arrange(chain, iteration)
}

# -----------------------------------------------------------------------------
# 10. Acceptance summary function
# -----------------------------------------------------------------------------
# Extracts and labels acceptance rates across chains at different fitting stages.

acceptance_summary <- function(fit, stage, burnin = 0.25) {
  chlaa_fit_report(fit, burnin = burnin)$acceptance_by_chain |>
    mutate(stage = stage, .before = 1)
}

# -----------------------------------------------------------------------------
# 11. ESS summary function
# -----------------------------------------------------------------------------
# Computes effective sample size per chain, translating autocorrelated draws
# into equivalent independent samples.

ess_summary <- function(fit, stage, burnin = 0.25, parameters = NULL) {
  dr <- draws_wide(fit, burnin = burnin, scale = "sampled")
  if (is.null(parameters)) {
    parameters <- setdiff(names(dr), c("chain", "iteration"))
  }

  dr |>
    group_split(chain) |>
    map_dfr(function(d) {
      ess <- coda::effectiveSize(as.matrix(d[, parameters, drop = FALSE]))
      tibble(chain = d$chain[[1]], parameter = names(ess), ess = as.numeric(ess))
    }) |>
    mutate(stage = stage, .before = 1)
}

# -----------------------------------------------------------------------------
# 12. R-hat summary function
# -----------------------------------------------------------------------------
# Computes Gelman-Rubin convergence diagnostics assessing whether chains from
# different starting points agree.

rhat_summary <- function(fit, stage, burnin = 0.25, parameters = NULL) {
  dr <- draws_wide(fit, burnin = burnin, scale = "sampled")
  if (is.null(parameters)) {
    parameters <- setdiff(names(dr), c("chain", "iteration"))
  }

  chains <- dr |>
    group_split(chain) |>
    map(\(d) coda::mcmc(as.matrix(d[, parameters, drop = FALSE])))

  psrf <- coda::gelman.diag(
    coda::mcmc.list(chains),
    autoburnin = FALSE,
    multivariate = FALSE
  )$psrf

  tibble(
    stage = stage,
    parameter = rownames(psrf),
    rhat = psrf[, "Point est."],
    rhat_upper = psrf[, "Upper C.I."]
  )
}

# -----------------------------------------------------------------------------
# 13. Parameter summary function
# -----------------------------------------------------------------------------
# Summarizes posterior credible intervals on natural scale, optionally checking
# whether true values are covered.

parameter_summary <- function(fit, stage, burnin = 0.25, truth = NULL) {
  out <- draws_wide(fit, burnin = burnin, scale = "natural") |>
    select(all_of(natural_fit_names)) |>
    pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
    group_by(parameter) |>
    summarise(
      q025 = quantile(value, 0.025),
      median = median(value),
      q975 = quantile(value, 0.975),
      .groups = "drop"
    ) |>
    mutate(stage = stage, .before = 1)

  if (!is.null(truth)) {
    out <- out |>
      left_join(tibble(parameter = names(truth), truth = as.numeric(truth)), by = "parameter") |>
      mutate(covers_truth = q025 <= truth & truth <= q975)
  }

  out
}

# -----------------------------------------------------------------------------
# 14. Positive definite covariance helper
# -----------------------------------------------------------------------------
# Ensures a covariance matrix is positive definite by flooring small
# eigenvalues and reconstructing.

make_pd <- function(x, min_eig = 1e-10) {
  x <- (x + t(x)) / 2
  eig <- eigen(x, symmetric = TRUE)
  eig$values <- pmax(eig$values, min_eig)
  out <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
  dimnames(out) <- dimnames(x)
  out
}

# -----------------------------------------------------------------------------
# 15. Extract proposal covariance from fit
# -----------------------------------------------------------------------------
# Learns proposal covariance from pilot chain output, scaling using the
# Roberts-Rosenthal adaptation formula.

proposal_from_fit <- function(fit, burnin = 0.25, scale = 1) {
  dr <- draws_wide(fit, burnin = burnin, scale = "sampled")
  theta <- as.matrix(dr[, setdiff(names(dr), c("chain", "iteration")), drop = FALSE])
  make_pd(cov(theta) * (2.38^2 / ncol(theta)) * scale)
}

# -----------------------------------------------------------------------------
# 16. Parameters from theta vector
# -----------------------------------------------------------------------------
# Converts parameter vector to full parameter object by unpacking and
# selectively updating.

pars_from_theta <- function(theta, template_pars, packer) {
  theta <- as.numeric(theta)
  names(theta) <- packer[["names"]]()
  unpacked <- packer[["unpack"]](theta)
  fixed_names <- names(packer[["inputs"]]()$fixed)
  update_names <- setdiff(names(unpacked), fixed_names)

  out <- template_pars
  for (nm in intersect(update_names, names(out))) {
    out[[nm]] <- unpacked[[nm]]
  }
  out
}

# -----------------------------------------------------------------------------
# 17. Chain median starting points
# -----------------------------------------------------------------------------
# Extracts posterior medians from each deterministic chain to initialize
# particle-filter chains.

chain_median_starts <- function(fit, template_pars, burnin = 0.25) {
  dr <- draws_wide(fit, burnin = burnin, scale = "sampled")
  packer <- attr(fit, "packer", exact = TRUE)
  fit_names <- packer[["names"]]()

  dr |>
    group_split(chain) |>
    map(function(d) {
      theta <- vapply(fit_names, \(nm) median(d[[nm]]), numeric(1))
      pars_from_theta(theta, template_pars, packer)
    })
}

# -----------------------------------------------------------------------------
# 18. Posterior predictive plot function
# -----------------------------------------------------------------------------
# Generates posterior predictive plots showing observation intervals and medians
# against observed data.

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

# -----------------------------------------------------------------------------
# 19. Raw-scale starting points
# -----------------------------------------------------------------------------
# Constructs three dispersed starting points on natural scale for the
# diagnostic raw-scale fit.

make_raw_start <- function(trans_prob, reporting_rate, obs_size, E0) {
  make_kirotshe_pars(
    trans_prob = trans_prob,
    reporting_rate = reporting_rate,
    obs_size = obs_size,
    E0 = E0
  )
}

raw_starts <- list(
  make_raw_start(trans_prob = 1.6e-3, reporting_rate = 0.12, obs_size = 20, E0 = 35),
  make_raw_start(trans_prob = 4.0e-4, reporting_rate = 0.25, obs_size = 200, E0 = 200),
  make_raw_start(trans_prob = 1.2e-3, reporting_rate = 0.15, obs_size = 45, E0 = 120)
)

bind_rows(
  tibble(row = "truth", parameter = names(truth_vec), value = as.numeric(truth_vec)),
  imap_dfr(raw_starts, \(pars, chain) tibble(
    row = paste0("chain_", chain),
    parameter = natural_fit_names,
    value = unlist(pars[natural_fit_names])
  ))
) |>
  pivot_wider(names_from = parameter, values_from = value)

# -----------------------------------------------------------------------------
# 20. Raw-scale prior and packer
# -----------------------------------------------------------------------------
# Defines prior distributions and proposal variances for raw-scale fitting.

raw_prior <- monty::monty_dsl(
  {
    trans_prob ~ Uniform(1e-4, 1e-2)
    reporting_rate ~ Uniform(0.06, 0.30)
    obs_size ~ Uniform(1, 300)
    E0 ~ Uniform(10, 800)
  },
  gradient = FALSE
)

raw_packer <- function(pars) {
  fixed <- pars[setdiff(names(pars), natural_fit_names)]
  monty::monty_packer(scalar = natural_fit_names, fixed = fixed)
}

raw_proposal <- c(
  trans_prob = 1e-8,
  reporting_rate = 2e-4,
  obs_size = 9,
  E0 = 400
)

# -----------------------------------------------------------------------------
# 21. Raw-scale deterministic fit
# -----------------------------------------------------------------------------
# Runs pMCMC with deterministic filter on raw parameter scale using synthetic
# data.

raw_fit <- chlaa_fit_pmcmc(
  data = fit_data_synthetic,
  pars = raw_starts[[1]],
  chain_pars = raw_starts,
  n_chains = n_chains,
  n_particles = 1,
  n_steps = pilot_steps,
  seed = 91,
  prior = raw_prior,
  packer = raw_packer(raw_starts[[1]]),
  proposal_var = raw_proposal,
  obs_interval = 7,
  time_start = time_start,
  deterministic = TRUE
)

# -----------------------------------------------------------------------------
# 22. Raw-scale diagnostics
# -----------------------------------------------------------------------------
# Trace plots showing poor mixing and low ESS values on raw scale.

acceptance_summary(raw_fit, "raw deterministic", burnin = 0.25)
ess_summary(raw_fit, "raw deterministic", burnin = 0.25)
rhat_summary(raw_fit, "raw deterministic", burnin = 0.25)

p_raw_trace <- chlaa_plot_trace(
  raw_fit,
  parameters = natural_fit_names,
  burnin = 0.25,
  scale = "natural"
)
print(p_raw_trace)

ggsave(
  file.path(fig_dir, "fitting_raw_scale_trace.png"),
  plot = p_raw_trace, width = 12, height = 8, dpi = 300
)

# -----------------------------------------------------------------------------
# 23. Raw-scale parameter pairs plot
# -----------------------------------------------------------------------------
# Scatter plot matrix of raw-scale posterior with truth overlaid.

p_raw_pairs <- chlaa_plot_parameter_pairs(
  raw_fit,
  parameters = natural_fit_names,
  burnin = 0.25,
  scale = "natural",
  truth = truth_vec,
  max_points = 2500
)
print(p_raw_pairs)

ggsave(
  file.path(fig_dir, "fitting_raw_scale_pairs.png"),
  plot = p_raw_pairs, width = 10, height = 10, dpi = 300
)

# -----------------------------------------------------------------------------
# 24. Transformed parameterization
# -----------------------------------------------------------------------------
# Implements log and logit transformations to improve parameter exploration
# geometry, with automatic back-transformation.
# log_trans_prob for a positive transmission probability;
# logit_reporting_rate for a probability bounded by zero and one;
# log_obs_size for a positive overdispersion parameter;
# log_E0 for a positive initial condition.
# The sampler moves on this transformed scale, but the plots below use the natural scale by default.

fit_names <- c(
  "log_trans_prob",
  "logit_reporting_rate",
  "log_obs_size",
  "log_E0"
)

fit_prior <- monty::monty_dsl(
  {
    log_trans_prob ~ Uniform(-9.21034, -4.60517) # log(1e-4) to log(1e-2)
    logit_reporting_rate ~ Uniform(-2.751535, -0.847298) # qlogis(c(0.06, 0.30))
    log_obs_size ~ Uniform(0, 5.703782)
    log_E0 ~ Uniform(2.302585, 6.684612)
  },
  gradient = FALSE
)

# Freeze-reporting_rate prior (3 params, for pilot/exploratory stage)
fit_prior_freeze <- monty::monty_dsl(
  {
    log_trans_prob ~ Uniform(-9.21034, -4.60517)
    log_obs_size ~ Uniform(0, 5.703782)
    log_E0 ~ Uniform(2.302585, 6.684612)
  },
  gradient = FALSE
)

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
  make_kirotshe_pars(
    trans_prob = trans_prob,
    reporting_rate = reporting_rate,
    obs_size = obs_size,
    E0 = E0
  ) |>
    add_transformed_values(freeze_reporting_rate = freeze_reporting_rate)
}

synthetic_starts <- list(
  make_start(trans_prob = 1.6e-3, reporting_rate = 0.12, obs_size = 20, E0 = 35),
  make_start(trans_prob = 4.0e-4, reporting_rate = 0.25, obs_size = 200, E0 = 200),
  make_start(trans_prob = 1.2e-3, reporting_rate = 0.15, obs_size = 45, E0 = 120)
)

# -----------------------------------------------------------------------------
# 25. Synthetic pilot run
# -----------------------------------------------------------------------------
# First transformed-scale fit with simple diagonal proposal to learn
# covariance structure. The transformed workflow begins with a deterministic pilot.
# The proposal is diagonal and deliberately simple; it only needs to move enough for us to learn the covariance structure.
# The pilot acceptance rates should not be interpreted as final diagnostics. They tell us whether the simple proposal
# was able to move enough to estimate a useful covariance.

pilot_proposal <- c(
  log_trans_prob = 0.02,
  logit_reporting_rate = 0.05,
  log_obs_size = 0.08,
  log_E0 = 0.08
)

synthetic_pilot <- chlaa_fit_pmcmc(
  data = fit_data_synthetic,
  pars = synthetic_starts[[1]],
  chain_pars = synthetic_starts,
  n_chains = n_chains,
  n_particles = 1,
  n_steps = pilot_steps,
  seed = 101,
  prior = fit_prior,
  packer = make_packer(synthetic_starts[[1]]),
  proposal_var = pilot_proposal,
  obs_interval = 7,
  time_start = time_start,
  deterministic = TRUE
)

acceptance_summary(synthetic_pilot, "synthetic pilot")

# -----------------------------------------------------------------------------
# 26. Learned proposal covariance
# -----------------------------------------------------------------------------
# Extracts correlation structure showing strong negative correlation between
# transmission and reporting parameters.
# The learned covariance is full rather than diagonal.
# Off-diagonal correlations are useful here because trans_prob, reporting_rate, and E0 can compensate for one another.

synthetic_det_proposal <- proposal_from_fit(
  synthetic_pilot,
  burnin = 0.25,
  scale = 1
)

round(cov2cor(synthetic_det_proposal), 2)

# -----------------------------------------------------------------------------
# 27. Synthetic deterministic main run
# -----------------------------------------------------------------------------
# Extended deterministic run with learned proposal, showing improved
# diagnostics and parameter recovery.
# With the learned covariance, the deterministic chains can be run longer.
# The main visual check is that chains from different starts overlap after burn-in.

synthetic_det <- chlaa_fit_pmcmc(
  data = fit_data_synthetic,
  pars = synthetic_starts[[1]],
  chain_pars = synthetic_starts,
  n_chains = n_chains,
  n_particles = 1,
  n_steps = deterministic_steps,
  seed = 201,
  prior = fit_prior,
  packer = make_packer(synthetic_starts[[1]]),
  proposal_var = synthetic_det_proposal,
  obs_interval = 7,
  time_start = time_start,
  deterministic = TRUE
)

acceptance_summary(synthetic_det, "synthetic deterministic")
ess_summary(synthetic_det, "synthetic deterministic")
rhat_summary(synthetic_det, "synthetic deterministic")

parameter_summary(
  synthetic_det,
  "synthetic deterministic",
  truth = truth_vec
)

# -----------------------------------------------------------------------------
# 28. Synthetic deterministic trace and pairs plots
# -----------------------------------------------------------------------------
# Trace, marginal distribution, and bivariate scatter plots verifying
# posterior geometry against synthetic truth.
# The deterministic fit is also useful for diagnosing posterior geometry.
# The dashed lines mark the synthetic truth.


p_synth_det_trace <- chlaa_plot_trace(
  synthetic_det,
  parameters = natural_fit_names,
  burnin = 0.25,
  scale = "natural"
)
print(p_synth_det_trace)

ggsave(
  file.path(fig_dir, "fitting_synthetic_det_trace.png"),
  plot = p_synth_det_trace, width = 12, height = 8, dpi = 300
)

p_synth_det_dist <- chlaa_plot_parameter_distributions(
  synthetic_det,
  parameters = natural_fit_names,
  burnin = 0.25,
  scale = "natural",
  truth = truth_vec
)
print(p_synth_det_dist)

ggsave(
  file.path(fig_dir, "fitting_synthetic_det_distributions.png"),
  plot = p_synth_det_dist, width = 12, height = 8, dpi = 300
)

p_synth_det_pairs <- chlaa_plot_parameter_pairs(
  synthetic_det,
  parameters = natural_fit_names,
  burnin = 0.25,
  scale = "natural",
  truth = truth_vec,
  max_points = 3000
)
print(p_synth_det_pairs)

ggsave(
  file.path(fig_dir, "fitting_synthetic_det_pairs.png"),
  plot = p_synth_det_pairs, width = 10, height = 10, dpi = 300
)

# -----------------------------------------------------------------------------
# 29. Synthetic particle filter starting points
# -----------------------------------------------------------------------------
# The particle-filter run uses the deterministic posterior medians as starting points.
# We keep the covariance direction learned by the deterministic fit and shrink the scale
# because the particle likelihood is noisy.

synthetic_particle_starts <- synthetic_starts

synthetic_particle_proposal <- proposal_from_fit(
  synthetic_det,
  burnin = 0.25,
  scale = 0.8
)

# -----------------------------------------------------------------------------
# 30. Synthetic particle filter run
# -----------------------------------------------------------------------------
# Final stochastic fit with 50 particles per likelihood evaluation across
# three chains.
# For the particle chains, acceptance is typically lower than the deterministic pilot
# because each likelihood evaluation is noisy.
# ESS and R-hat are complementary. ESS asks how many independent samples the autocorrelated chain is worth;
# R-hat asks whether the different chains agree.
# The synthetic parameters are recovered if the posterior intervals cover the true values and the posterior mass is concentrated near them.

synthetic_particle <- chlaa_fit_pmcmc(
  data = fit_data_synthetic,
  pars = synthetic_particle_starts[[1]],
  chain_pars = synthetic_particle_starts,
  n_chains = n_chains,
  n_particles = particle_count,
  n_steps = particle_steps,
  seed = 301,
  prior = fit_prior,
  packer = make_packer(synthetic_particle_starts[[1]]),
  proposal_var = synthetic_particle_proposal,
  obs_interval = 7,
  time_start = time_start,
  deterministic = FALSE
)

acceptance_summary(synthetic_particle, "synthetic particle")
ess_summary(synthetic_particle, "synthetic particle")
rhat_summary(synthetic_particle, "synthetic particle")

parameter_summary(
  synthetic_particle,
  "synthetic particle",
  truth = truth_vec
)

# -----------------------------------------------------------------------------
# 31. Synthetic particle filter diagnostics and plots
# -----------------------------------------------------------------------------
# Trace, distribution, and pair plots for particle chains; posterior
# predictive fit with 95% and 50% credible intervals.

p_synth_part_trace <- chlaa_plot_trace(
  synthetic_particle,
  parameters = natural_fit_names,
  burnin = 0.25,
  scale = "natural"
)
print(p_synth_part_trace)

ggsave(
  file.path(fig_dir, "fitting_synthetic_particle_trace.png"),
  plot = p_synth_part_trace, width = 12, height = 8, dpi = 300
)

p_synth_part_dist <- chlaa_plot_parameter_distributions(
  synthetic_particle,
  parameters = natural_fit_names,
  burnin = 0.25,
  scale = "natural",
  truth = truth_vec
)
print(p_synth_part_dist)

ggsave(
  file.path(fig_dir, "fitting_synthetic_particle_distributions.png"),
  plot = p_synth_part_dist, width = 12, height = 8, dpi = 300
)

p_synth_part_pairs <- chlaa_plot_parameter_pairs(
  synthetic_particle,
  parameters = natural_fit_names,
  burnin = 0.25,
  scale = "natural",
  truth = truth_vec,
  max_points = 3000
)
print(p_synth_part_pairs)

ggsave(
  file.path(fig_dir, "fitting_synthetic_particle_pairs.png"),
  plot = p_synth_part_pairs, width = 10, height = 10, dpi = 300
)

p_synth_part_fit <- plot_case_fit(
  synthetic_particle,
  synthetic_weekly,
  "Posterior predictive fit to synthetic weekly cases",
  seed = 401
)
print(p_synth_part_fit)

ggsave(
  file.path(fig_dir, "fitting_synthetic_posterior_predictive.png"),
  plot = p_synth_part_fit, width = 10, height = 6, dpi = 300
)

# -----------------------------------------------------------------------------
# 32. Real data starting points
# -----------------------------------------------------------------------------
# Defines three dispersed starting points for real Kirotshe data fitting.
# The pilot stage freezes reporting_rate to let other parameters settle first.

real_starts_freeze <- list(
  make_start(trans_prob = 1e-3, reporting_rate = 0.15, obs_size = 30, E0 = 80,
             freeze_reporting_rate = TRUE),
  make_start(trans_prob = 5e-4, reporting_rate = 0.15, obs_size = 20, E0 = 40,
             freeze_reporting_rate = TRUE),
  make_start(trans_prob = 5e-3, reporting_rate = 0.15, obs_size = 100, E0 = 200,
             freeze_reporting_rate = TRUE)
)

real_starts_full <- list(
  make_start(trans_prob = 1e-3, reporting_rate = 0.15, obs_size = 30, E0 = 80),
  make_start(trans_prob = 5e-4, reporting_rate = 0.08, obs_size = 20, E0 = 40),
  make_start(trans_prob = 5e-3, reporting_rate = 0.20, obs_size = 100, E0 = 200)
)

# -----------------------------------------------------------------------------
# 33. Real data pilot (freeze reporting_rate)
# -----------------------------------------------------------------------------
# Pilot run with reporting_rate frozen to let trans_prob, obs_size, and E0
# settle before unfreezing all four parameters.

freeze_pilot_proposal <- c(
  log_trans_prob = 0.02,
  log_obs_size = 0.08,
  log_E0 = 0.08
)

real_pilot <- chlaa_fit_pmcmc(
  data = fit_data_real,
  pars = real_starts_freeze[[1]],
  chain_pars = real_starts_freeze,
  n_chains = n_chains,
  n_particles = 1,
  n_steps = pilot_steps,
  seed = 501,
  prior = fit_prior_freeze,
  packer = make_packer(real_starts_freeze[[1]], freeze_reporting_rate = TRUE),
  proposal_var = freeze_pilot_proposal,
  obs_interval = 7,
  time_start = time_start,
  deterministic = TRUE
)

acceptance_summary(real_pilot, "real pilot (freeze)")

# Learn covariance from frozen pilot and expand 3x3 → 4x4
pilot_cov3 <- proposal_from_fit(real_pilot, burnin = 0.25, scale = 1)
real_det_proposal <- matrix(0, 4, 4)
idx_map <- c(1L, 3L, 4L) # log_trans_prob, log_obs_size, log_E0 → positions in full
real_det_proposal[idx_map, idx_map] <- pilot_cov3
real_det_proposal[2, 2] <- 0.05 * (2.38^2 / 4) # logit_reporting_rate default
diag(real_det_proposal) <- pmax(diag(real_det_proposal), 1e-6)

# Warm-start full starts from pilot medians
real_starts_full <- chain_median_starts(
  real_pilot,
  template_pars = real_starts_full[[1]],
  burnin = 0.25
)

# -----------------------------------------------------------------------------
# 34. Real data deterministic run (full 4-param)
# -----------------------------------------------------------------------------
# Extended deterministic run on real data with all four parameters unfrozen
# and learned proposal covariance expanded from the frozen pilot.

real_det <- chlaa_fit_pmcmc(
  data = fit_data_real,
  pars = real_starts_full[[1]],
  chain_pars = real_starts_full,
  n_chains = n_chains,
  n_particles = 1,
  n_steps = deterministic_steps,
  seed = 601,
  prior = fit_prior,
  packer = make_packer(real_starts_full[[1]]),
  proposal_var = real_det_proposal,
  obs_interval = 7,
  time_start = time_start,
  deterministic = TRUE
)

acceptance_summary(real_det, "real deterministic")
ess_summary(real_det, "real deterministic")
rhat_summary(real_det, "real deterministic")
parameter_summary(real_det, "real deterministic")

# -----------------------------------------------------------------------------
# 35. Real data deterministic plots
# -----------------------------------------------------------------------------
# Trace and pair plots for deterministic chains on real data.

p_real_det_trace <- chlaa_plot_trace(
  real_det,
  parameters = natural_fit_names,
  burnin = 0.25,
  scale = "natural"
)
print(p_real_det_trace)

ggsave(
  file.path(fig_dir, "fitting_real_det_trace.png"),
  plot = p_real_det_trace, width = 12, height = 8, dpi = 300
)

p_real_det_pairs <- chlaa_plot_parameter_pairs(
  real_det,
  parameters = natural_fit_names,
  burnin = 0.25,
  scale = "natural",
  max_points = 3000
)
print(p_real_det_pairs)

ggsave(
  file.path(fig_dir, "fitting_real_det_pairs.png"),
  plot = p_real_det_pairs, width = 10, height = 10, dpi = 300
)

# -----------------------------------------------------------------------------
# 36. Real data particle filter setup
# -----------------------------------------------------------------------------
# Extracts posterior medians and applies larger proposal shrinkage for real
# data's tighter posterior.

real_particle_starts <- chain_median_starts(
  real_det,
  template_pars = real_starts_full[[1]],
  burnin = 0.25
)

real_particle_proposal <- proposal_from_fit(
  real_det,
  burnin = 0.25,
  scale = 0.2
)

# -----------------------------------------------------------------------------
# 37. Real data particle filter run
# -----------------------------------------------------------------------------
# Final stochastic fit on real data using 50 particles.

real_particle <- chlaa_fit_pmcmc(
  data = fit_data_real,
  pars = real_particle_starts[[1]],
  chain_pars = real_particle_starts,
  n_chains = n_chains,
  n_particles = particle_count,
  n_steps = particle_steps,
  seed = 701,
  prior = fit_prior,
  packer = make_packer(real_particle_starts[[1]]),
  proposal_var = real_particle_proposal,
  obs_interval = 7,
  time_start = time_start,
  deterministic = FALSE
)

acceptance_summary(real_particle, "real particle")
ess_summary(real_particle, "real particle")
rhat_summary(real_particle, "real particle")
parameter_summary(real_particle, "real particle")

# -----------------------------------------------------------------------------
# 38. Real data particle filter diagnostics and plots
# -----------------------------------------------------------------------------
# Trace, distribution, pair plots, and posterior predictive fit for real data
# chains.

p_real_part_trace <- chlaa_plot_trace(
  real_particle,
  parameters = natural_fit_names,
  burnin = 0.25,
  scale = "natural"
)
print(p_real_part_trace)

ggsave(
  file.path(fig_dir, "fitting_real_particle_trace.png"),
  plot = p_real_part_trace, width = 12, height = 8, dpi = 300
)

p_real_part_dist <- chlaa_plot_parameter_distributions(
  real_particle,
  parameters = natural_fit_names,
  burnin = 0.25,
  scale = "natural"
)
print(p_real_part_dist)

ggsave(
  file.path(fig_dir, "fitting_real_particle_distributions.png"),
  plot = p_real_part_dist, width = 12, height = 8, dpi = 300
)

p_real_part_pairs <- chlaa_plot_parameter_pairs(
  real_particle,
  parameters = natural_fit_names,
  burnin = 0.25,
  scale = "natural",
  max_points = 3000
)
print(p_real_part_pairs)

ggsave(
  file.path(fig_dir, "fitting_real_particle_pairs.png"),
  plot = p_real_part_pairs, width = 10, height = 10, dpi = 300
)

p_real_part_fit <- plot_case_fit(
  real_particle,
  kirotshe,
  "Posterior predictive fit to Kirotshe weekly cases",
  seed = 801
)
print(p_real_part_fit)

ggsave(
  file.path(fig_dir, "fitting_real_posterior_predictive.png"),
  plot = p_real_part_fit, width = 10, height = 6, dpi = 300
)

# -----------------------------------------------------------------------------
# 39. Save fitted artifact
# -----------------------------------------------------------------------------
# Packages the fitted posterior and metadata for reuse in downstream scenario
# and health economics analyses.

fit_artifact <- list(
  fit = real_particle,
  pars = attr(real_particle, "start_pars", exact = TRUE),
  observed = kirotshe,
  interventions = kirotshe_meta,
  outbreak_start = outbreak_start,
  outbreak_end = outbreak_end,
  burnin = 0.25,
  particle_count = particle_count
)

artifact_path <- file.path(data_dir, "kirotshe_particle_fit.rds")
saveRDS(fit_artifact, artifact_path)
artifact_path

message("Fitting script complete. All figures saved to: ", fig_dir)

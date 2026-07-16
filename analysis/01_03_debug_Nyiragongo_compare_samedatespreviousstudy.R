# =============================================================================
# Nyiragongo Production Fit — comparison with previous study
# =============================================================================
#
# Fits the chlaa cholera model to Nyiragongo health zone data using pMCMC,
# matching the current pipeline model choices for comparison against a
# previously published analysis of the same outbreak.
#
# Model choices (matching 01_02_fitting_all_HZs.R):
#   - frac_neff fitted (N_eff = frac_neff * pop), reporting_rate fixed at 0.30
#   - contam_half_sat decoupled from N_eff (census-based)
#   - seed_state quasi-equilibrium initial conditions from E0
#   - immunity_asym=280, seek_mild=0.1, seek_severe=0.85
#   - R0-based starting points
#
# Steps (matching 01_02 2-stage pipeline):
#   1. Load data and set up parameters with current defaults + interventions
#   2. Exploratory fit (100 particles, 1000 steps, 3 chains)
#   3. Learn covariance, warm-start from exploratory
#   4. Production fit (200 particles, 10000 steps, 3 chains)
#   5. Diagnostics, plots, and save
# =============================================================================

library(chlaa)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(posterior)

fig_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/figures"
data_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/analysis/data"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load Nyiragongo data
# -----------------------------------------------------------------------------

hz_name <- "nyiragongo"
outbreak_start <- as.Date("2022-08-01")
outbreak_end <- as.Date("2024-09-30")

nyiragongo <- read_csv(file.path(data_dir, "IDSR_dataset.csv"), show_col_types = FALSE) |>
  filter(hz == hz_name) |>
  mutate(date = as.Date(date)) |>
  filter(date >= outbreak_start, date <= outbreak_end) |>
  arrange(date) |>
  mutate(time = seq_len(n()) * 7L) |>
  select(date, time, cases, deaths, population)

cat("Nyiragongo outbreak data:\n")
cat("  Weeks:", nrow(nyiragongo), "\n")
cat("  Total cases:", sum(nyiragongo$cases), "\n")
cat("  Population:", nyiragongo$population[1], "\n")

time_start <- -21L
pop_hz <- nyiragongo$population[1]

# Global constants (matching 01_02 pipeline)
H_REF <- 1.0
POP_REF <- 516000
RR_FIXED <- 0.30
E0_MAX <- 800
K_R0 <- POP_REF * 5.1446e-3 / H_REF

# Quasi-equilibrium seeding
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

# -----------------------------------------------------------------------------
# 2. Intervention dates (from CERF/WHO reports)
# -----------------------------------------------------------------------------

orc_start_day <- as.integer(as.Date("2023-01-16") - outbreak_start)
orc_end_day <- as.integer(as.Date("2023-07-15") - outbreak_start)
ctc_start_day <- as.integer(as.Date("2023-01-16") - outbreak_start)
ctc_end_day <- as.integer(as.Date("2023-07-15") - outbreak_start)
chlor_start_day <- as.integer(as.Date("2023-01-30") - outbreak_start)
chlor_end_day <- as.integer(as.Date("2023-07-30") - outbreak_start)
hyg_start_day <- as.integer(as.Date("2023-01-30") - outbreak_start)
hyg_end_day <- as.integer(as.Date("2023-07-30") - outbreak_start)
cati_start_day <- as.integer(as.Date("2023-01-30") - outbreak_start)
cati_end_day <- as.integer(as.Date("2023-07-30") - outbreak_start)
# No latrines in this period
lat_start_day <- 0L
lat_end_day <- 0L

# Vaccination: Loo supplementary, RDC MoH 2023; 264,824 doses, single-dose 5-day campaign
vax1_start_day <- as.integer(as.Date("2023-01-23") - outbreak_start)
vax1_end_day <- as.integer(as.Date("2023-01-27") - outbreak_start)
vax1_total_doses <- 264824L

cat(sprintf("  ORC: days %d–%d\n", orc_start_day, orc_end_day))
cat(sprintf("  CTC: days %d–%d\n", ctc_start_day, ctc_end_day))
cat(sprintf("  Chlor/Hyg/CATI: days %d–%d\n", chlor_start_day, chlor_end_day))
cat(sprintf("  Vax1: days %d–%d (%d doses)\n", vax1_start_day, vax1_end_day, vax1_total_doses))

# Build vaccination schedule arrays
vax1_n_days <- vax1_end_day - vax1_start_day + 1L
vax1_daily <- vax1_total_doses / vax1_n_days
vax1_sched_time <- c(as.integer(time_start), as.integer(vax1_start_day), as.integer(vax1_end_day + 1L))
vax1_sched_doses <- c(0, vax1_daily, 0)

# Empty vax2 schedule
vax2_sched_time <- c(as.integer(time_start), as.integer(time_start + 1L))
vax2_sched_doses <- c(0, 0)

# -----------------------------------------------------------------------------
# 3. Parameter factory (current pipeline model choices + interventions)
# -----------------------------------------------------------------------------

make_nyiragongo_pars <- function(trans_prob, obs_size, E0, frac_neff = 0.10) {
  N_eff <- frac_neff * pop_hz
  h <- H_REF * (pop_hz / POP_REF)
  out <- chlaa_parameters(
    N = N_eff,
    contact_rate = 0,
    contam_half_sat = h,
    trans_prob = trans_prob,
    E0 = E0,
    incubation_time = 4.845,
    duration_sym = 14.48,
    immunity_asym = 280,
    seek_mild = 0.1,
    seek_severe = 0.85,
    reporting_rate = RR_FIXED,
    obs_size = obs_size,
    fatality_treated = 0.0021,
    fatality_untreated = 0.5,
    death_reporting_rate = 0.5,
    obs_size_deaths = 1.0,
    # Interventions
    orc_start = orc_start_day,
    orc_end = orc_end_day,
    ctc_start = ctc_start_day,
    ctc_end = ctc_end_day,
    chlor_start = chlor_start_day,
    chlor_end = chlor_end_day,
    chlor_effect = 0.10,
    hyg_start = hyg_start_day,
    hyg_end = hyg_end_day,
    hyg_effect = 0.10,
    cati_start = cati_start_day,
    cati_end = cati_end_day,
    cati_effect = 0.10,
    lat_start = lat_start_day,
    lat_end = lat_end_day,
    lat_effect = 0,
    vax1_start = vax1_start_day,
    vax1_end = vax1_end_day,
    vax1_total_doses = vax1_total_doses
  )
  ss <- seed_state(E0, out)
  for (nm in names(ss)) out[[nm]] <- ss[[nm]]
  out$frac_neff <- frac_neff
  # Vaccination schedule arrays (required by odin model)
  out$vax1_schedule_time <- vax1_sched_time
  out$vax1_schedule_doses <- vax1_sched_doses
  out$n_vax1_schedule <- length(vax1_sched_time)
  out$vax2_schedule_time <- vax2_sched_time
  out$vax2_schedule_doses <- vax2_sched_doses
  out$n_vax2_schedule <- length(vax2_sched_time)
  out
}

# -----------------------------------------------------------------------------
# 4. Fitting configuration (matching 01_02 pipeline)
# -----------------------------------------------------------------------------

natural_fit_names <- c("trans_prob", "obs_size", "E0", "frac_neff")
fit_names <- c("log_trans_prob", "log_obs_size", "log_E0", "logit_frac_neff")

n_explore <- 100L
n_explore_steps <- 1000L
n_prod <- 200L
n_prod_steps <- 10000L
seed_explore <- 42L
seed_prod <- 123L

fit_data <- nyiragongo |> select(time, cases, deaths)

# Dynamic E0 initialization
expected_reporting_rate <- 0.10
E0_val <- ceiling(max(5, mean(nyiragongo$cases[1:min(3, nrow(nyiragongo))])) / expected_reporting_rate)
E0_val <- min(E0_val, 0.9 * E0_MAX)
E0_val <- max(10, E0_val)
cat(sprintf("  Initial E0: %d\n", E0_val))

# -----------------------------------------------------------------------------
# 5. Helper functions
# -----------------------------------------------------------------------------

draws_wide <- function(fit, burnin = 0.25, scale = c("sampled", "natural")) {
  scale <- match.arg(scale)
  chlaa_fit_trace(fit, burnin = burnin, scale = scale) |>
    pivot_wider(names_from = parameter, values_from = value) |>
    arrange(chain, iteration)
}

pars_from_theta <- function(theta, template_pars, packer) {
  theta <- as.numeric(theta)
  names(theta) <- packer[["names"]]()
  unpacked <- packer[["unpack"]](theta)
  fixed_names <- names(packer[["inputs"]]()$fixed)
  update_names <- setdiff(names(unpacked), fixed_names)
  out <- template_pars
  for (nm in intersect(update_names, names(out))) out[[nm]] <- unpacked[[nm]]
  out
}

chain_median_starts <- function(fit, template_pars, burnin = 0.25) {
  dr <- draws_wide(fit, burnin = burnin, scale = "sampled")
  packer <- attr(fit, "packer", exact = TRUE)
  fit_names_local <- packer[["names"]]()
  dr |>
    group_split(chain) |>
    map(function(d) {
      theta <- vapply(fit_names_local, \(nm) median(d[[nm]]), numeric(1))
      pars_from_theta(theta, template_pars, packer)
    })
}

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
    geom_ribbon(data = fit_cases, aes(date, ymin = q0p025, ymax = q0p975),
                fill = "#911e12", alpha = 0.25) +
    geom_ribbon(data = fit_cases, aes(date, ymin = q0p25, ymax = q0p75),
                fill = "#e2b19b", alpha = 0.45) +
    geom_line(data = fit_cases, aes(date, q0p5), colour = "#2f6a4e", linewidth = 0.8) +
    geom_line(data = observed, aes(date, cases), colour = "#494949", linewidth = 0.5) +
    labs(x = NULL, y = "Weekly reported cases", title = title) +
    theme_minimal()
}

# -----------------------------------------------------------------------------
# 6. Priors and parameterization (current pipeline choices)
# -----------------------------------------------------------------------------

fit_prior <- monty::monty_dsl({
  log_trans_prob ~ Uniform(-9.21034, -2.995732)   # log(c(1e-4, 5e-2))
  log_obs_size ~ Uniform(0, 5.703782)             # log(c(1, 300))
  log_E0 ~ Uniform(2.302585, 6.684612)            # log(c(10, 800))
  logit_frac_neff ~ Uniform(-4.6, 2.944439)       # ~qlogis(c(0.01, 0.95))
}, gradient = FALSE)

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
  make_nyiragongo_pars(
    trans_prob = trans_prob,
    obs_size = obs_size,
    E0 = E0,
    frac_neff = frac_neff
  ) |>
    add_transformed_values()
}

# -----------------------------------------------------------------------------
# 7. R0-based starting points (all 4 params)
# -----------------------------------------------------------------------------

r0_targets  <- c(1.5, 2.5, 4.0)
frac_starts <- c(0.10, 0.05, 0.20)
tp_starts   <- r0_targets / (frac_starts * K_R0)

starts <- list(
  make_start(tp_starts[1], 30, E0_val, frac_neff = frac_starts[1]),
  make_start(tp_starts[2], 20, max(10, round(E0_val * 0.5)), frac_neff = frac_starts[2]),
  make_start(tp_starts[3], 100, min(0.9 * E0_MAX, max(10, round(E0_val * 1.5))), frac_neff = frac_starts[3])
)

# =============================================================================
# FITTING PIPELINE (2-stage, matching 01_02_fitting_all_HZs.R)
# =============================================================================

# -----------------------------------------------------------------------------
# 8. Exploratory fit (100 particles, 1000 steps, 3 chains)
# -----------------------------------------------------------------------------

cat("\n=== EXPLORATORY FIT ===\n")

explore_proposal <- matrix(0, 4, 4)
explore_proposal[1, 1] <- 0.02  # log_trans_prob
explore_proposal[2, 2] <- 0.08  # log_obs_size
explore_proposal[3, 3] <- 0.08  # log_E0
explore_proposal[4, 4] <- 0.10  # logit_frac_neff

fit_packer_stage1 <- make_packer(starts[[1]])

fit_explore <- chlaa_fit_pmcmc(
  data = fit_data,
  pars = starts[[1]],
  chain_pars = starts,
  n_chains = length(starts),
  n_particles = n_explore,
  n_steps = n_explore_steps,
  seed = seed_explore,
  prior = fit_prior,
  packer = fit_packer_stage1,
  proposal_var = explore_proposal,
  obs_interval = 7,
  time_start = time_start
)

report_explore <- chlaa_fit_report(fit_explore, burnin = 0.25, thin = 2)
cat("Exploratory acceptance rate:", report_explore$acceptance_rate, "\n")
print(report_explore$posterior_summary)

# -----------------------------------------------------------------------------
# 9. Learn covariance and prepare production
# -----------------------------------------------------------------------------

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
pars_warm <- starts[[1]]
for (nm in names(warm_vec)) pars_warm[[nm]] <- warm_vec[[nm]]

fit_starts_stage2 <- chain_median_starts(fit_explore, template_pars = starts[[1]])
fit_packer_stage2 <- make_packer(fit_starts_stage2[[1]])

prod_proposal <- tryCatch(
  {
    cov_mat <- cov(pooled) * (2.38^2 / d)
    eig <- eigen(cov_mat, symmetric = TRUE)
    if (any(eig$values < 1e-12) || any(!is.finite(eig$values))) {
      cat("Warning: Degenerate covariance matrix. Using diagonal fallback.\n")
      stop("Degenerate covariance")
    }
    eig$values <- pmax(eig$values, 1e-12)
    eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  },
  error = function(e) {
    cat("Using diagonal proposal (covariance learning failed).\n")
    diag(pmax(apply(pooled, 2, var), 1e-6)) * (2.38^2 / d)
  }
)
diag(prod_proposal) <- pmax(diag(prod_proposal), 1e-6)

cat("\nWarm-start parameter values:\n")
for (nm in packer$names()) {
  cat(sprintf("  %s = %.6f\n", nm, pars_warm[[nm]]))
}

rm(fit_explore, report_explore, packer, pooled)
gc()

# -----------------------------------------------------------------------------
# 10. Production fit (200 particles, 10000 steps, 3 chains)
# -----------------------------------------------------------------------------

cat("\n=== PRODUCTION FIT ===\n")

fit <- chlaa_fit_pmcmc(
  data = fit_data,
  pars = pars_warm,
  chain_pars = fit_starts_stage2,
  n_chains = length(fit_starts_stage2),
  n_particles = n_prod,
  n_steps = n_prod_steps,
  seed = seed_prod,
  prior = fit_prior,
  packer = fit_packer_stage2,
  proposal_var = prod_proposal,
  obs_interval = 7,
  time_start = time_start
)

report_prod <- chlaa_fit_report(fit, burnin = 0.25, thin = 2)
cat("Production acceptance rate:", report_prod$acceptance_rate, "\n")
print(report_prod$posterior_summary)

if (report_prod$acceptance_rate < 0.10) {
  cat("WARNING: Low acceptance rate may indicate identifiability issues.\n")
}

# R-hat and ESS
pars_arr <- fit$pars
dimnames(pars_arr) <- list(attr(fit, "packer")$names(), NULL, NULL)
draws_diag <- posterior::as_draws_array(aperm(pars_arr, c(2, 3, 1)))
rhat_ess <- posterior::summarise_draws(draws_diag, "rhat", "ess_bulk", "ess_tail")
cat("\n=== R-hat and ESS ===\n")
print(rhat_ess)
rm(pars_arr, draws_diag)

# -----------------------------------------------------------------------------
# 11. Diagnostic plots
# -----------------------------------------------------------------------------

p_trace <- chlaa_plot_trace(fit, parameters = natural_fit_names, burnin = 0.25, scale = "natural")
ggsave(file.path(fig_dir, "fitting_nyiragongo_comparative_trace.png"),
       p_trace, width = 12, height = 8, dpi = 300)

p_dist <- chlaa_plot_parameter_distributions(
  fit, parameters = natural_fit_names, burnin = 0.25, scale = "natural"
)
ggsave(file.path(fig_dir, "fitting_nyiragongo_comparative_distributions.png"),
       p_dist, width = 12, height = 8, dpi = 300)

p_fit <- plot_case_fit(
  fit, nyiragongo,
  "Nyiragongo weekly cases (08-2022 to 09-2024)",
  seed = seed_prod
)
print(p_fit)
ggsave(file.path(fig_dir, "fitting_nyiragongo_comparative_fit.png"),
       plot = p_fit, width = 10, height = 6, dpi = 300)

cat("\nPlots saved to:", fig_dir, "\n")

# -----------------------------------------------------------------------------
# 12. Save fitted artifact
# -----------------------------------------------------------------------------

rds_dir <- file.path(fig_dir, ".rds files")
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

fit_artifact <- list(
  hz_name = hz_name,
  pop_hz = pop_hz,
  fit = fit,
  pars_start = fit_starts_stage2,
  pars_warm = pars_warm,
  report = report_prod,
  rhat_ess = rhat_ess,
  observed = nyiragongo,
  outbreak_start = outbreak_start,
  outbreak_end = outbreak_end,
  total_cases = sum(nyiragongo$cases),
  n_weeks = nrow(nyiragongo),
  fitted_parameters = natural_fit_names,
  notes = "Current pipeline choices: frac_neff fitted, RR=0.30 fixed, seed_state, immunity_asym=280, seek_mild=0.1, seek_severe=0.85, contam_half_sat census-based",
  timestamp = Sys.time()
)

saveRDS(fit_artifact, file.path(rds_dir, "nyiragongo_comparative_fit.rds"))

cat("Fit artifact saved to:", file.path(rds_dir, "nyiragongo_comparative_fit.rds"), "\n")
message("Nyiragongo comparative fitting complete.")

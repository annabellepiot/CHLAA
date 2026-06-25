# =============================================================================
# Scenario Workflow Vignette
# Reproduces: https://ojwatson.github.io/chlaa/articles/scenario_workflow.html
#
# This script performs scenario analysis starting from a fitted pMCMC model.
# It defines no-intervention, anticipatory-action, and anticipatory-action +
# vaccination scenarios, then produces posterior scenario forecasts and
# decision summary tables.
# =============================================================================

library(chlaa)
library(ggplot2)
library(dplyr)

# Directory for saving figures
fig_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/CHLAA/figures"

# -----------------------------------------------------------------------------
# 1. Load the fitted model
# -----------------------------------------------------------------------------
# Scenario analysis starts from a model that can reproduce the outbreak we are
# using as the baseline. This reads the pre-fitted stochastic pMCMC model for
# Kirotshe from package data rather than refitting.

data_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/analysis/data"
fit_obj <- readRDS(file.path(data_dir, "kirotshe_particle_fit.rds"))

fit <- fit_obj$fit
base_pars <- fit_obj$pars
observed <- fit_obj$observed
burnin <- fit_obj$burnin

pars_fit <- chlaa_update_from_fit(
  fit = fit,
  pars = base_pars,
  draw = "median",
  burnin = burnin
)

# -----------------------------------------------------------------------------
# 2. Posterior predictive check
# -----------------------------------------------------------------------------
# The posterior predictive check shows the fitted distribution for weekly
# reported cases. The points are the observed Kirotshe data; the ribbons and
# line come from posterior draws.

fit_fc <- chlaa_forecast_from_fit(
  fit = fit,
  pars = base_pars,
  time = observed$time,
  vars = "inc_symptoms_weekly",
  include_cases = TRUE,
  obs_interval = 7,
  obs_model = "nbinom",
  n_draws = 80,
  burnin = burnin,
  seed = 11,
  dt = 1
)

p_posterior_check <- chlaa_plot_forecast(
  fit_fc,
  var = "cases",
  data = observed,
  data_y = "cases"
)
print(p_posterior_check)

ggsave(
  file.path(fig_dir, "scenario_posterior_predictive_check.png"),
  plot = p_posterior_check, width = 10, height = 6, dpi = 300
)

# -----------------------------------------------------------------------------
# 3. Define scenario parameters
# -----------------------------------------------------------------------------
# The baseline scenario is the fitted outbreak with the recorded intervention
# timings. We compare it with one no-intervention counterfactual and two
# anticipatory-action examples.

horizon <- max(observed$time) + 182
scenario_time <- seq(7, horizon, by = 7)
trigger_time <- min(observed$time[observed$cases >= 50])
response_end <- horizon + 1
campaign_days <- 28
vaccine_doses <- floor(0.20 * pars_fit$N)

# No intervention counterfactual: all interventions turned off
no_intervention <- list(
  chlor_start = 0, chlor_end = 0, chlor_effect = 0,
  hyg_start = 0, hyg_end = 0, hyg_effect = 0,
  lat_start = 0, lat_end = 0, lat_effect = 0,
  cati_start = 0, cati_end = 0, cati_effect = 0,
  orc_start = 0, orc_end = 0, orc_capacity = 0,
  ctc_start = 0, ctc_end = 0, ctc_capacity = 0,
  vax1_start = 0, vax1_end = 0, vax1_doses_per_day = 0, vax1_total_doses = 0,
  vax2_start = 0, vax2_end = 0, vax2_doses_per_day = 0, vax2_total_doses = 0
)

# Anticipatory-action response: WASH + treatment triggered at 50-case threshold
aa_response <- list(
  chlor_start = trigger_time, chlor_end = response_end, chlor_effect = pars_fit$chlor_effect,
  hyg_start = trigger_time, hyg_end = response_end, hyg_effect = pars_fit$hyg_effect,
  lat_start = trigger_time, lat_end = response_end, lat_effect = max(pars_fit$lat_effect, 0.10),
  cati_start = trigger_time, cati_end = response_end, cati_effect = pars_fit$cati_effect,
  orc_start = trigger_time, orc_end = response_end, orc_capacity = pars_fit$orc_capacity,
  ctc_start = trigger_time, ctc_end = response_end, ctc_capacity = pars_fit$ctc_capacity,
  vax1_start = 0, vax1_end = 0, vax1_doses_per_day = 0, vax1_total_doses = 0,
  vax2_start = 0, vax2_end = 0, vax2_doses_per_day = 0, vax2_total_doses = 0
)

# Anticipatory-action response + vaccination campaign
aa_vaccination <- modifyList(aa_response, list(
  vax1_start = trigger_time + 14,
  vax1_end = trigger_time + 14 + campaign_days,
  vax1_total_doses = vaccine_doses,
  vax1_doses_per_day = vaccine_doses / campaign_days
))

scenarios <- list(
  chlaa_scenario("no_interventions", no_intervention),
  chlaa_scenario("aa_response", aa_response),
  chlaa_scenario("aa_response_plus_vaccine", aa_vaccination)
)

vapply(scenarios, `[[`, character(1), "name")

# -----------------------------------------------------------------------------
# 4. Posterior scenario forecasts
# -----------------------------------------------------------------------------
# chlaa_forecast_scenarios_from_fit() uses the same posterior draws for each
# scenario, then reports absolute forecasts and paired differences against the
# baseline. That pairing matters because it removes some Monte Carlo noise
# from scenario differences.

scenario_fc <- chlaa_forecast_scenarios_from_fit(
  fit = fit,
  pars = base_pars,
  scenarios = scenarios,
  baseline_name = "fitted_response",
  time = scenario_time,
  vars = c("inc_symptoms_weekly", "cum_symptoms", "cum_deaths"),
  include_cases = TRUE,
  obs_interval = 7,
  n_draws = 60,
  burnin = burnin,
  seed = 12,
  dt = 1
)

# Absolute scenario forecasts: weekly cases across all intervention scenarios
p_scenario_absolute <- chlaa_plot_scenario_forecasts(
  scenario_fc,
  var = "cases",
  type = "absolute",
  data = observed,
  data_y = "cases"
)
print(p_scenario_absolute)

ggsave(
  file.path(fig_dir, "scenario_absolute_forecasts_cases.png"),
  plot = p_scenario_absolute, width = 12, height = 7, dpi = 300
)

# Difference plot: cumulative deaths relative to baseline
p_scenario_diff_deaths <- chlaa_plot_scenario_forecasts(
  scenario_fc,
  var = "cum_deaths",
  type = "difference",
  include_baseline = FALSE
)
print(p_scenario_diff_deaths)

ggsave(
  file.path(fig_dir, "scenario_difference_cumulative_deaths.png"),
  plot = p_scenario_diff_deaths, width = 12, height = 7, dpi = 300
)

# -----------------------------------------------------------------------------
# 5. Decision summary
# -----------------------------------------------------------------------------
# For compact decision tables we simulate the same scenarios using the
# posterior median parameter set. This is faster and easier to inspect, while
# the forecast plots above show posterior uncertainty.

scenario_runs <- chlaa_run_scenarios(
  pars = pars_fit,
  scenarios = c(list(chlaa_scenario("fitted_response", list())), scenarios),
  time = scenario_time,
  n_particles = 50,
  dt = 1,
  seed = 13
)

# Summary table with total cases, deaths, cases/deaths averted, peak
# incidence, time to peak, and time to control
scenario_summary <- chlaa_scenario_summary(
  scenario_runs,
  baseline = "fitted_response",
  incidence_var = "inc_symptoms_weekly"
)
print(scenario_summary)

# Comprehensive comparison including health economic metrics (costs, DALYs,
# ICERs, net monetary benefit)
scenario_comparison <- chlaa_compare_scenarios(
  scenario_runs,
  baseline = "fitted_response",
  include_econ = TRUE,
  wtp = 1500
)
print(scenario_comparison)

message("Scenario workflow script complete. All figures saved to: ", fig_dir)

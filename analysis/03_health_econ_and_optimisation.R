# =============================================================================
# Health Economics and Optimisation Vignette
# Reproduces: https://ojwatson.github.io/chlaa/articles/health_econ_and_optimisation.html
#
# This script performs health-economic analysis using the saved Kirotshe pMCMC
# fit. It covers: posterior predictive checks, intervention scenario simulation,
# cost-effectiveness analysis, cost-effectiveness acceptability curves (CEAC),
# and budget optimisation across vaccination, WASH, and care allocations.
# =============================================================================

library(chlaa)
library(ggplot2)
library(dplyr)

# Directory for saving figures
fig_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/CHLAA/figures"

# -----------------------------------------------------------------------------
# 1. Load fitted model
# -----------------------------------------------------------------------------
# The health-economic workflow uses the same saved Kirotshe pMCMC fit as the
# scenario vignette. We first check the fitted weekly case curve, then use the
# posterior-median parameter set for scenario simulation and costing.

fit_obj <- readRDS(system.file(
  "extdata", "kirotshe_particle_fit.rds",
  package = "chlaa",
  mustWork = TRUE
))

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
# 2. Forecast and visualise fitted case curve
# -----------------------------------------------------------------------------
# Generates forecast from the fitted model and plots the weekly case curve.

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
  seed = 21,
  dt = 1
)

p_fitted_cases <- chlaa_plot_forecast(
  fit_fc,
  var = "cases",
  data = observed,
  data_y = "cases"
)
print(p_fitted_cases)

ggsave(
  file.path(fig_dir, "health_econ_fitted_case_curve.png"),
  plot = p_fitted_cases, width = 10, height = 6, dpi = 300
)

# -----------------------------------------------------------------------------
# 3. Set up intervention scenarios
# -----------------------------------------------------------------------------
# The economic comparison uses one fitted-response baseline, a no-intervention
# counterfactual, and two anticipatory-action examples. Costs and DALYs are
# then computed from the simulated cumulative cases, deaths, treatments, WASH
# windows, and vaccine doses.

horizon <- max(observed$time) + 182
scenario_time <- seq(7, horizon, by = 7)
trigger_time <- min(observed$time[observed$cases >= 50])
response_end <- horizon + 1
campaign_days <- 28
vaccine_doses <- floor(0.20 * pars_fit$N)

# No intervention counterfactual
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

# Anticipatory-action response
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
  chlaa_scenario("fitted_response", list()),
  chlaa_scenario("no_interventions", no_intervention),
  chlaa_scenario("aa_response", aa_response),
  chlaa_scenario("aa_response_plus_vaccine", aa_vaccination)
)

# Run all scenarios
runs <- chlaa_run_scenarios(
  pars = pars_fit,
  scenarios = scenarios,
  time = scenario_time,
  n_particles = 50,
  dt = 1,
  seed = 22
)

# Overlay of weekly incidence across scenarios
p_scenario_overlay <- chlaa_plot_scenario_overlay(
  runs,
  var = "inc_symptoms_weekly"
)
print(p_scenario_overlay)

ggsave(
  file.path(fig_dir, "health_econ_scenario_overlay.png"),
  plot = p_scenario_overlay, width = 12, height = 7, dpi = 300
)

# -----------------------------------------------------------------------------
# 4. Economic analysis and cost-effectiveness
# -----------------------------------------------------------------------------
# The default economics inputs are deliberately transparent package data. They
# are useful for workflow development and should be replaced with local values
# for decision-making.

econ <- chlaa_econ_defaults()
head(chlaa_econ_sources(), 8)

cmp <- chlaa_compare_scenarios(
  runs,
  baseline = "fitted_response",
  include_econ = TRUE,
  econ = econ,
  wtp = 1500
)

print(cmp)

# Cost-effectiveness plane
p_ce_plane <- chlaa_plot_ce_plane(cmp)
print(p_ce_plane)

ggsave(
  file.path(fig_dir, "health_econ_cost_effectiveness_plane.png"),
  plot = p_ce_plane, width = 10, height = 8, dpi = 300
)

# Deaths averted by scenario
p_deaths_averted <- chlaa_plot_scenarios(cmp, metric = "deaths_averted")
print(p_deaths_averted)

ggsave(
  file.path(fig_dir, "health_econ_deaths_averted.png"),
  plot = p_deaths_averted, width = 10, height = 6, dpi = 300
)

# -----------------------------------------------------------------------------
# 5. Cost-effectiveness acceptability curves (CEAC)
# -----------------------------------------------------------------------------
# The CEAC estimates which scenario has the highest net monetary benefit across
# particles at each willingness-to-pay threshold.

ceac_tbl <- chlaa_ceac(
  runs,
  baseline = "fitted_response",
  wtp = seq(0, 3000, by = 250),
  econ = econ
)

head(ceac_tbl)

p_ceac <- ggplot(ceac_tbl, aes(x = wtp, y = prob_best, colour = scenario)) +
  geom_line(linewidth = 0.8) +
  labs(
    x = "WTP (USD per DALY averted)",
    y = "Probability cost-effective",
    colour = "Scenario",
    title = "Cost-effectiveness acceptability curves"
  ) +
  theme_minimal()
print(p_ceac)

ggsave(
  file.path(fig_dir, "health_econ_ceac.png"),
  plot = p_ceac, width = 10, height = 6, dpi = 300
)

# -----------------------------------------------------------------------------
# 6. Budget optimisation
# -----------------------------------------------------------------------------
# The optimiser is a simple search over vaccination, WASH, and care
# allocations. It is meant as a transparent starting point: the key object is
# the evaluations table, which shows how simulated outcomes change across
# feasible allocations.

opt <- chlaa_optimise_budget(
  pars = pars_fit,
  budget = 5e5,
  time = seq(trigger_time, horizon, by = 7),
  n_particles = 20,
  dt = 1,
  grid_size = 8,
  min_fraction = list(vax = 0.1),
  max_fraction = list(wash = 0.6),
  max_vax_doses_per_day = 5000,
  max_total_doses = vaccine_doses
)

print(opt$best)

p_budget_surface <- ggplot(opt$evaluations, aes(x = budget_vax, y = budget_wash, colour = deaths)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_viridis_c() +
  labs(
    x = "Budget allocated to vaccination",
    y = "Budget allocated to WASH",
    colour = "Expected deaths",
    title = "Budget allocation search surface"
  ) +
  theme_minimal()
print(p_budget_surface)

ggsave(
  file.path(fig_dir, "health_econ_budget_allocation_surface.png"),
  plot = p_budget_surface, width = 10, height = 8, dpi = 300
)

message("Health economics script complete. All figures saved to: ", fig_dir)

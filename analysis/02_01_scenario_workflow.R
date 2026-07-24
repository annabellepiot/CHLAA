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
fig_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

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

# Helper to generate a uniform vaccination schedule for scenario use.
# The odin model requires interpolated schedule arrays, not just
# vax1_start/end/total_doses.
make_scenario_vax_schedule <- function(total_doses, start_day, end_day) {
  n_days <- max(end_day - start_day, 1L)
  daily_doses <- total_doses / n_days
  sched_time <- as.integer(c(start_day, end_day))
  sched_doses <- c(daily_doses, 0)
  # Ensure schedule covers early times for interpolation
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

# Empty vaccination schedule arrays (for scenarios with no vaccination)
empty_vax <- list(
  vax1_schedule_time = c(0L, 1L),
  vax1_schedule_doses = c(0, 0),
  n_vax1_schedule = 2L,
  vax2_schedule_time = c(0L, 1L),
  vax2_schedule_doses = c(0, 0),
  n_vax2_schedule = 2L
)

# No intervention counterfactual: all interventions turned off
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

# Anticipatory-action response: WASH + treatment triggered at 50-case threshold
aa_response <- c(list(
  chlor_start = trigger_time, chlor_end = response_end, chlor_effect = pars_fit$chlor_effect,
  hyg_start = trigger_time, hyg_end = response_end, hyg_effect = pars_fit$hyg_effect,
  lat_start = trigger_time, lat_end = response_end, lat_effect = max(pars_fit$lat_effect, 0.10),
  cati_start = trigger_time, cati_end = response_end, cati_effect = pars_fit$cati_effect,
  orc_start = trigger_time, orc_end = response_end, orc_capacity = pars_fit$orc_capacity,
  ctc_start = trigger_time, ctc_end = response_end, ctc_capacity = pars_fit$ctc_capacity,
  vax1_start = 0, vax1_end = 0, vax1_total_doses = 0,
  vax2_start = 0, vax2_end = 0, vax2_total_doses = 0
), empty_vax)

# Anticipatory-action response + vaccination campaign
vax1_start_day <- trigger_time + 14
vax1_end_day <- vax1_start_day + campaign_days
vax_sched <- make_scenario_vax_schedule(vaccine_doses, vax1_start_day, vax1_end_day)
aa_vaccination <- modifyList(aa_response, c(list(
  vax1_start = vax1_start_day,
  vax1_end = vax1_end_day,
  vax1_total_doses = vaccine_doses
), vax_sched))

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

# Cumulative excess cases & deaths at time snapshots
# Two separate 2x2 faceted plots: one for cases, one for deaths

# Find nearest available weekly time points to 100, 200, 300, 400 days
target_days <- c(100, 200, 300, 400)
snap_times <- vapply(target_days, function(d) {
  scenario_time[which.min(abs(scenario_time - d))]
}, numeric(1))

# Scenario colours and ordering
scenario_colours <- c(
  "aa_response"              = "#cb86ff",
  "aa_response_plus_vaccine" = "#88b517",
  "no_interventions"         = "#f7776d"
)
baseline_colour <- "#0abfc6"
scenario_order <- c("aa_response", "aa_response_plus_vaccine", "no_interventions")
scenario_labels <- c(
  "aa_response"               = "AA response",
  "aa_response_plus_vaccine"  = "AA response\n+ vaccination",
  "no_interventions"          = "No\ninterventions"
)

# Legend labels (single-line for legend)
scenario_legend_labels <- c(
  "aa_response"               = "AA response",
  "aa_response_plus_vaccine"  = "AA response + vaccination",
  "no_interventions"          = "No interventions"
)

# Helper: build one excess plot (cases or deaths)
build_excess_plot <- function(var_name, y_label) {
  dat <- scenario_fc %>%
    filter(type == "difference", variable == var_name,
           time %in% snap_times, scenario != "fitted_response") %>%
    mutate(
      facet = factor(paste0(target_days[match(time, snap_times)], " days"),
                     levels = paste0(target_days, " days")),
      scenario = factor(scenario, levels = scenario_order)
    )

  # Numeric labels: "median (95% UI lower to upper)"
  dat <- dat %>%
    mutate(num_label = paste0(
      round(q0p5), "\n(",
      round(q0p025), " to ",
      round(q0p975), ")"
    ))

  # Position labels above/below the whisker depending on sign
  dat <- dat %>%
    mutate(label_y = ifelse(q0p5 >= 0, q0p975, q0p025),
           label_vjust = ifelse(q0p5 >= 0, -0.15, 1.15))

  # "Fitted response" label: placed at right side of 200 days facet
  label_df <- data.frame(
    facet = factor("200 days", levels = paste0(target_days, " days")),
    x = "no_interventions", y = 0, label = "Fitted response "
  )
  label_df$x <- factor(label_df$x, levels = scenario_order)

  # Compute a tighter y-axis range based on data
  y_lo <- min(dat$q0p025, na.rm = TRUE)
  y_hi <- max(dat$q0p975, na.rm = TRUE)
  y_pad <- (y_hi - y_lo) * 0.25
  y_limits <- c(y_lo - y_pad, y_hi + y_pad)

  ggplot(dat, aes(x = scenario)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = baseline_colour,
               linewidth = 0.8) +
    geom_crossbar(aes(y = q0p5, ymin = q0p25, ymax = q0p75, fill = scenario),
                  width = 0.55, alpha = 0.8, colour = "grey30",
                  linewidth = 0.3, middle.linewidth = 0.6) +
    geom_errorbar(aes(ymin = q0p025, ymax = q0p975),
                  width = 0.28, linewidth = 0.4, colour = "grey30") +
    geom_text(aes(y = label_y, label = num_label, vjust = label_vjust),
              size = 3.5, family = "Helvetica", lineheight = 0.85) +
    geom_text(data = label_df, aes(x = x, y = y, label = label),
              hjust = 1.1, vjust = -0.6, size = 3.5, colour = baseline_colour,
              fontface = "bold", family = "Helvetica") +
    facet_wrap(~ facet, nrow = 2) +
    scale_x_discrete(labels = NULL, drop = FALSE) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = 0.02)) +
    scale_fill_manual(
      values = scenario_colours,
      labels = scenario_legend_labels
    ) +
    labs(
      x = NULL, y = y_label,
      title = "Kirotshe",
      subtitle = "Days since start of modelled outbreak",
      fill = NULL
    ) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.9))) +
    theme_minimal(base_family = "Helvetica", base_size = 13) +
    theme(
      panel.grid         = element_blank(),
      panel.background   = element_rect(fill = "white", colour = "grey70"),
      panel.border       = element_rect(fill = NA, colour = "grey70",
                                        linewidth = 0.5),
      plot.background    = element_rect(fill = "white", colour = NA),
      strip.text         = element_text(face = "bold", size = 12),
      axis.line          = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.ticks.y       = element_line(colour = "grey40"),
      axis.text.x        = element_blank(),
      legend.position    = "right",
      legend.text        = element_text(size = 10),
      plot.title         = element_text(face = "bold", size = 14),
      plot.subtitle      = element_text(size = 11, colour = "grey40"),
      panel.spacing      = unit(1.2, "lines")
    )
}

# Cases plot
p_excess_cases <- build_excess_plot(
  "cum_symptoms",
  "Excess cumulative cases"
)
print(p_excess_cases)

ggsave(file.path(fig_dir, "scenario_excess_cases_kirotshe.png"),
       plot = p_excess_cases, width = 10, height = 9, dpi = 300)
ggsave(file.path(fig_dir, "scenario_excess_cases_kirotshe.pdf"),
       plot = p_excess_cases, width = 10, height = 9)

# Deaths plot
p_excess_deaths <- build_excess_plot(
  "cum_deaths",
  "Excess cumulative deaths"
)
print(p_excess_deaths)

ggsave(file.path(fig_dir, "scenario_excess_deaths_kirotshe.png"),
       plot = p_excess_deaths, width = 10, height = 9, dpi = 300)
ggsave(file.path(fig_dir, "scenario_excess_deaths_kirotshe.pdf"),
       plot = p_excess_deaths, width = 10, height = 9)

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



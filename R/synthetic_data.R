# Synthetic example data helpers

.chlaa_case_study_day <- function(date, start_date) {
  as.numeric(as.Date(date) - as.Date(start_date))
}

.chlaa_case_study_base_parameters <- function(obs_size) {
  # Paper-aligned defaults (population and broad epidemiological scale):
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC12477517/
  chlaa_parameters(
    N = 540000,
    Sev0 = 2,
    contact_rate = 0.05,
    trans_prob = 0.055,
    reporting_rate = 0.12,
    obs_size = obs_size,
    chlor_start = 0, chlor_end = 0, chlor_effect = 0,
    hyg_start = 0, hyg_end = 0, hyg_effect = 0,
    lat_start = 0, lat_end = 0, lat_effect = 0,
    cati_start = 0, cati_end = 0, cati_effect = 0,
    orc_start = 0, orc_end = 0,
    ctc_start = 0, ctc_end = 0,
    vax1_start = 0, vax1_end = 0, vax1_total_doses = 0,
    vax2_start = 0, vax2_end = 0, vax2_total_doses = 0
  )
}

.chlaa_case_study_scenarios <- function(pars,
                                        start_date,
                                        trigger_date,
                                        declaration_date,
                                        late_vax_start_date,
                                        vax_total_doses,
                                        campaign_days) {
  trigger_day <- .chlaa_case_study_day(trigger_date, start_date)
  declaration_day <- .chlaa_case_study_day(declaration_date, start_date)
  late_vax_start_day <- .chlaa_case_study_day(late_vax_start_date, start_date)

  response_duration <- 240
  response_end <- declaration_day + response_duration

  response_effects <- list(
    chlor_end = response_end, chlor_effect = 0.45,
    hyg_end = response_end, hyg_effect = 0.45,
    lat_end = response_end, lat_effect = 0.25,
    cati_end = response_end, cati_effect = 0.25,
    orc_end = response_end,
    ctc_end = response_end
  )

  response_from <- function(start_day) {
    c(
      list(
        chlor_start = start_day,
        hyg_start = start_day,
        lat_start = start_day,
        cati_start = start_day,
        orc_start = start_day,
        ctc_start = start_day
      ),
      response_effects
    )
  }

  late_campaign <- list(
    vax1_start = late_vax_start_day,
    vax1_end = late_vax_start_day + campaign_days,
    vax1_total_doses = vax_total_doses,
    vax2_start = 0, vax2_end = 0,
    vax2_total_doses = 0
  )

  early_one_dose <- list(
    vax1_start = trigger_day + 30,
    vax1_end = trigger_day + 30 + campaign_days,
    vax1_total_doses = vax_total_doses,
    vax2_start = 0, vax2_end = 0,
    vax2_total_doses = 0
  )

  # Keep the second-dose campaign conservative to avoid infeasible
  # administration trajectories in stochastic particles when many first-dose
  # recipients leave eligibility compartments before dose 2.
  vax2_total_stable <- min(vax_total_doses / 2, 40000)

  early_two_dose <- list(
    vax1_start = trigger_day + 30,
    vax1_end = trigger_day + 30 + campaign_days,
    vax1_total_doses = vax_total_doses / 2,
    vax2_start = trigger_day + 30 + campaign_days + 14,
    vax2_end = trigger_day + 30 + campaign_days + 14 + campaign_days,
    vax2_total_doses = vax2_total_stable
  )

  list(
    chlaa_scenario("scenario_1_baseline", c(response_from(declaration_day), late_campaign)),
    chlaa_scenario("scenario_2_anticipatory_action", c(response_from(trigger_day), late_campaign)),
    chlaa_scenario("scenario_3_anticipatory_action_plus_one_vaccine_dose", c(response_from(trigger_day), early_one_dose)),
    chlaa_scenario("scenario_4_anticipatory_action_plus_two_vaccine_doses", c(response_from(trigger_day), early_two_dose))
  )
}

#' Generate A Synthetic Cholera Outbreak Time Series
#'
#' Generates a paper-aligned synthetic outbreak curve from the package model and
#' returns the observed case data frame. This is intended as a small example data
#' generator rather than a shared vignette setup object.
#'
#' @param time Numeric vector of simulation times (days).
#' @param start_date Start date corresponding to `time = 0`.
#' @param trigger_date Anticipatory action trigger date.
#' @param declaration_date Outbreak declaration date.
#' @param late_vax_start_date Planned campaign start date for the generated
#'   response scenario.
#' @param seed Integer random seed.
#' @param n_particles Number of particles used to generate latent incidence.
#' @param dt Model time step.
#' @param obs_size Observation over-dispersion (Negative Binomial size).
#' @param vax_total_doses Total vaccine doses used in the generated response
#'   scenario.
#' @param campaign_days Vaccination campaign duration (days).
#'
#' @return A data.frame with columns `date`, `time`, `cases`, `mu_cases`,
#'   `inc_symptoms_truth`, and `inc_infections_truth`.
#'   The generating parameter set is attached as `attr(x, "truth_parameters")`.
#' @export
chlaa_generate_example_outbreak_data <- function(time = 0:915,
                                                 start_date = as.Date("2022-07-01"),
                                                 trigger_date = as.Date("2022-10-25"),
                                                 declaration_date = as.Date("2022-12-14"),
                                                 late_vax_start_date = as.Date("2023-01-20"),
                                                 seed = 42,
                                                 n_particles = 20,
                                                 dt = 1,
                                                 obs_size = 18,
                                                 vax_total_doses = 280000,
                                                 campaign_days = 150) {
  if (!is.numeric(time) || length(time) < 2) {
    stop("time must be a numeric vector with length >= 2", call. = FALSE)
  }
  time <- sort(as.numeric(time))

  date_args <- list(
    start_date = start_date,
    trigger_date = trigger_date,
    declaration_date = declaration_date,
    late_vax_start_date = late_vax_start_date
  )
  for (nm in names(date_args)) {
    if (!inherits(date_args[[nm]], "Date") || length(date_args[[nm]]) != 1) {
      stop(nm, " must be a single Date", call. = FALSE)
    }
  }
  if (!(trigger_date < declaration_date && declaration_date <= late_vax_start_date)) {
    stop("Expected trigger_date < declaration_date <= late_vax_start_date", call. = FALSE)
  }

  pars <- .chlaa_case_study_base_parameters(obs_size = obs_size)
  scenarios <- .chlaa_case_study_scenarios(
    pars = pars,
    start_date = start_date,
    trigger_date = trigger_date,
    declaration_date = declaration_date,
    late_vax_start_date = late_vax_start_date,
    vax_total_doses = vax_total_doses,
    campaign_days = campaign_days
  )

  truth <- chlaa_run_scenarios(
    pars = pars,
    scenarios = scenarios[1],
    time = time,
    n_particles = n_particles,
    dt = dt,
    seed = seed
  )

  truth_sym <- stats::aggregate(truth$inc_symptoms, by = list(time = truth$time), FUN = mean)
  truth_inf <- stats::aggregate(truth$inc_infections, by = list(time = truth$time), FUN = mean)

  mu_cases <- pmax(pars$reporting_rate * truth_sym$x, 0.01)

  set.seed(seed)
  cases <- stats::rnbinom(length(mu_cases), mu = mu_cases, size = obs_size)

  out <- data.frame(
    date = start_date + truth_sym$time,
    time = truth_sym$time,
    cases = as.integer(cases),
    mu_cases = as.numeric(mu_cases),
    inc_symptoms_truth = as.numeric(truth_sym$x),
    inc_infections_truth = as.numeric(truth_inf$x),
    stringsAsFactors = FALSE
  )

  attr(out, "truth_parameters") <- pars
  attr(out, "paper_reference") <- "https://pmc.ncbi.nlm.nih.gov/articles/PMC12477517/"
  out
}

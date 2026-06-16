# Scenario builder utilities (anticipatory action)

#' Derive an anticipatory action trigger time from a simulation
#'
#' Finds the first time where an aggregated incidence stream reaches or exceeds a threshold.
#'
#' @param sim Output of `chlaa_simulate()` (or a single-scenario subset of `chlaa_run_scenarios()`).
#' @param threshold Numeric threshold.
#' @param var Variable name in `sim` to threshold on. Defaults to `inc_symptoms`.
#' @param fun Aggregation function across particles at each time (defaults to `mean`).
#'
#' @return A single numeric trigger time (same units as `sim$time`), or `NA_real_` if never reached.
#' @export
chlaa_trigger_time_from_sim <- function(sim, threshold, var = "inc_symptoms", fun = mean) {
  if (!is.data.frame(sim)) stop("sim must be a data.frame", call. = FALSE)
  if (!all(c("time", "particle", var) %in% names(sim))) {
    stop("sim must contain columns: time, particle, ", var, call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1) stop("threshold must be a single number", call. = FALSE)

  agg <- tapply(sim[[var]], sim$time, fun)
  tt <- as.numeric(names(agg))
  idx <- which(agg >= threshold)
  if (length(idx) == 0) return(NA_real_)
  min(tt[idx], na.rm = TRUE)
}

.chlaa_window <- function(start, duration, horizon = NULL) {
  if (is.null(horizon)) {
    return(list(start = start, end = start + duration))
  }
  list(start = start, end = min(start + duration, horizon))
}

.normalise_regimen <- function(regimen) {
  r <- tolower(as.character(regimen[1]))
  if (r %in% c("1", "1dose", "one", "single")) return("1dose")
  if (r %in% c("2", "2dose", "two", "double")) return("2dose")
  stop("Invalid regimen: ", regimen, " (use 1/2 or '1dose'/'2dose')", call. = FALSE)
}

.chlaa_make_vax_plan <- function(N,
                                  trigger_time,
                                  total_doses,
                                  regimen = c("1dose", "2dose"),
                                  delay = 0,
                                  campaign_days = 6,
                                  dose_interval = 14,
                                  horizon = NULL) {
  regimen <- .normalise_regimen(regimen)

  if (!is.numeric(total_doses) || length(total_doses) != 1 || total_doses < 0) {
    stop("total_doses must be a single non-negative number", call. = FALSE)
  }
  if (!is.numeric(campaign_days) || length(campaign_days) != 1 || campaign_days <= 0) {
    stop("campaign_days must be > 0", call. = FALSE)
  }

  start1 <- trigger_time + delay
  end1 <- start1 + campaign_days
  if (!is.null(horizon)) end1 <- min(end1, horizon)

  if (regimen == "1dose") {
    doses1 <- total_doses
    rate1 <- doses1 / campaign_days
    return(list(
      vax1_start = start1, vax1_end = end1,
      vax1_total_doses = doses1, vax1_doses_per_day = rate1,
      vax2_start = 0, vax2_end = 0,
      vax2_total_doses = 0, vax2_doses_per_day = 0
    ))
  }

  people <- floor(total_doses / 2)
  doses1 <- people
  doses2 <- people

  rate1 <- doses1 / campaign_days

  start2 <- end1 + dose_interval
  end2 <- start2 + campaign_days
  if (!is.null(horizon)) {
    start2 <- min(start2, horizon)
    end2 <- min(end2, horizon)
  }
  rate2 <- doses2 / campaign_days

  list(
    vax1_start = start1, vax1_end = end1,
    vax1_total_doses = doses1, vax1_doses_per_day = rate1,
    vax2_start = start2, vax2_end = end2,
    vax2_total_doses = doses2, vax2_doses_per_day = rate2
  )
}

#' Build a standard set of anticipatory action scenarios
#'
#' @param pars Baseline parameter list.
#' @param trigger_time Numeric day the AA package starts. If `NULL`, derived from `baseline_sim`
#'   and `trigger_threshold`.
#' @param baseline_sim Baseline simulation output used to derive trigger time when `trigger_time` is `NULL`.
#' @param trigger_threshold Threshold for `inc_symptoms` used to derive trigger time.
#' @param horizon Optional cap on intervention end times (typically max simulation time + 1).
#' @param include_no_intervention Include a "no_intervention" scenario that zeros response levers.
#' @param include_aa_no_vax Include an "aa_no_vax" scenario (WASH + care, no vaccination).
#' @param wash_duration Duration in days for WASH components.
#' @param care_duration Duration in days for ORC/CTC availability.
#' @param chlor_effect,hyg_effect,lat_effect,cati_effect Fractional effects (0-1).
#' @param orc_capacity,ctc_capacity Capacities (persons/day) when active.
#' @param vax_total_doses Total doses available.
#' @param vax_delay Days from trigger to start vaccination.
#' @param vax_campaign_days Campaign duration for each round.
#' @param vax_dose_interval Days between rounds (2-dose scenario).
#' @param ve_1,ve_2 Vaccine efficacies.
#' @param vax_immunity_1,vax_immunity_2 Vaccine immunity durations in days.
#' @param baseline_name Name for baseline scenario.
#'
#' @return A list of `chlaa_scenario` objects.
#' @export
chlaa_make_aa_scenarios <- function(pars,
                                      trigger_time = NULL,
                                      baseline_sim = NULL,
                                      trigger_threshold = NULL,
                                      horizon = NULL,
                                      include_no_intervention = TRUE,
                                      include_aa_no_vax = TRUE,
                                      wash_duration = 120,
                                      care_duration = 120,
                                      chlor_effect = 0.20,
                                      hyg_effect = 0.20,
                                      lat_effect = 0.10,
                                      cati_effect = 0.10,
                                      orc_capacity = NULL,
                                      ctc_capacity = NULL,
                                      vax_total_doses = 0,
                                      vax_delay = 0,
                                      vax_campaign_days = 6,
                                      vax_dose_interval = 14,
                                      ve_1 = NULL,
                                      ve_2 = NULL,
                                      vax_immunity_1 = NULL,
                                      vax_immunity_2 = NULL,
                                      baseline_name = "baseline") {
  .check_named_list(pars, "pars")

  if (is.null(trigger_time)) {
    if (is.null(baseline_sim) || is.null(trigger_threshold)) {
      stop("If trigger_time is NULL you must provide baseline_sim and trigger_threshold", call. = FALSE)
    }
    trigger_time <- chlaa_trigger_time_from_sim(baseline_sim, threshold = trigger_threshold)
    if (is.na(trigger_time)) stop("Trigger threshold was never reached in baseline_sim", call. = FALSE)
  }

  scenarios <- list(chlaa_scenario(baseline_name, list()))

  if (isTRUE(include_no_intervention)) {
    scenarios <- c(scenarios, list(
      chlaa_scenario("no_intervention", list(
        chlor_start = 0, chlor_end = 0, chlor_effect = 0,
        hyg_start = 0, hyg_end = 0, hyg_effect = 0,
        lat_start = 0, lat_end = 0, lat_effect = 0,
        cati_start = 0, cati_end = 0, cati_effect = 0,
        orc_start = 0, orc_end = 0,
        ctc_start = 0, ctc_end = 0,
        orc_capacity = 0,
        ctc_capacity = 0,
        vax1_start = 0, vax1_end = 0, vax1_doses_per_day = 0, vax1_total_doses = 0,
        vax2_start = 0, vax2_end = 0, vax2_doses_per_day = 0, vax2_total_doses = 0
      ))
    ))
  }

  w_wash <- .chlaa_window(trigger_time, wash_duration, horizon = horizon)
  w_care <- .chlaa_window(trigger_time, care_duration, horizon = horizon)

  orc_cap_use <- if (is.null(orc_capacity)) pars$orc_capacity else orc_capacity
  ctc_cap_use <- if (is.null(ctc_capacity)) pars$ctc_capacity else ctc_capacity

  ve1_use <- if (is.null(ve_1)) pars$ve_1 else ve_1
  ve2_use <- if (is.null(ve_2)) pars$ve_2 else ve_2
  vim1_use <- if (is.null(vax_immunity_1)) pars$vax_immunity_1 else vax_immunity_1
  vim2_use <- if (is.null(vax_immunity_2)) pars$vax_immunity_2 else vax_immunity_2

  aa_common <- list(
    chlor_start = w_wash$start, chlor_end = w_wash$end, chlor_effect = chlor_effect,
    hyg_start = w_wash$start, hyg_end = w_wash$end, hyg_effect = hyg_effect,
    lat_start = w_wash$start, lat_end = w_wash$end, lat_effect = lat_effect,
    cati_start = w_wash$start, cati_end = w_wash$end, cati_effect = cati_effect,
    orc_start = w_care$start, orc_end = w_care$end, orc_capacity = orc_cap_use,
    ctc_start = w_care$start, ctc_end = w_care$end, ctc_capacity = ctc_cap_use
  )

  if (isTRUE(include_aa_no_vax)) {
    scenarios <- c(scenarios, list(
      chlaa_scenario("aa_no_vax", c(
        aa_common,
        list(
          vax1_start = 0, vax1_end = 0, vax1_doses_per_day = 0, vax1_total_doses = 0,
          vax2_start = 0, vax2_end = 0, vax2_doses_per_day = 0, vax2_total_doses = 0,
          ve_1 = ve1_use, ve_2 = ve2_use,
          vax_immunity_1 = vim1_use, vax_immunity_2 = vim2_use
        )
      ))
    ))
  }

  if (vax_total_doses > 0) {
    N <- pars$N

    plan1 <- .chlaa_make_vax_plan(
      N = N, trigger_time = trigger_time, total_doses = vax_total_doses,
      regimen = "1dose", delay = vax_delay, campaign_days = vax_campaign_days,
      dose_interval = vax_dose_interval, horizon = horizon
    )
    scenarios <- c(scenarios, list(
      chlaa_scenario("aa_vax_1dose", c(
        aa_common,
        plan1,
        list(ve_1 = ve1_use, ve_2 = ve2_use,
             vax_immunity_1 = vim1_use, vax_immunity_2 = vim2_use)
      ))
    ))

    plan2 <- .chlaa_make_vax_plan(
      N = N, trigger_time = trigger_time, total_doses = vax_total_doses,
      regimen = "2dose", delay = vax_delay, campaign_days = vax_campaign_days,
      dose_interval = vax_dose_interval, horizon = horizon
    )
    scenarios <- c(scenarios, list(
      chlaa_scenario("aa_vax_2dose", c(
        aa_common,
        plan2,
        list(ve_1 = ve1_use, ve_2 = ve2_use,
             vax_immunity_1 = vim1_use, vax_immunity_2 = vim2_use)
      ))
    ))
  }

  scenarios
}

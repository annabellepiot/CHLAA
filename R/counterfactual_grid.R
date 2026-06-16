# Counterfactual grid utilities

.chlaa_num_tag <- function(x, digits = 3) {
  if (is.integer(x)) return(as.character(x))
  if (!is.numeric(x)) return(as.character(x))
  s <- format(round(x, digits = digits), trim = TRUE, scientific = FALSE)
  s <- gsub("-", "m", s)
  s <- gsub("\\.", "p", s)
  s
}

.chlaa_make_name <- function(prefix, row) {
  parts <- c(prefix)
  for (nm in names(row)) {
    parts <- c(parts, paste0(nm, .chlaa_num_tag(row[[nm]])))
  }
  nm <- paste(parts, collapse = "_")
  gsub("[^A-Za-z0-9_]+", "", nm)
}

#' Build a counterfactual grid of AA-like scenarios
#'
#' This creates a grid of scenarios by crossing selected policy knobs.
#'
#' @param trigger_time Numeric trigger day.
#' @param horizon End cap for intervention windows (typically max(time)+1).
#' @param aa_start_offset Vector of offsets (days) added to trigger_time (negative = earlier).
#' @param wash_duration Vector of durations (days).
#' @param care_duration Vector of durations (days).
#' @param chlor_effect,hyg_effect,lat_effect,cati_effect Vectors of fractional effects.
#' @param orc_capacity,ctc_capacity Vectors of capacities.
#' @param vax_regimen Vector of regimens: "none", "1dose", "2dose".
#' @param vax_total_doses Vector of total doses.
#' @param vax_delay Vector of delays from AA start (days).
#' @param vax_campaign_days Campaign duration (days) per round.
#' @param vax_dose_interval Interval between rounds (2-dose) in days.
#' @param prefix Scenario name prefix.
#'
#' @return A data.frame with columns: scenario, modify, and knob columns.
#' @export
chlaa_counterfactual_grid <- function(trigger_time,
                                       horizon,
                                       aa_start_offset = 0,
                                       wash_duration = 120,
                                       care_duration = 120,
                                       chlor_effect = 0.2,
                                       hyg_effect = 0.2,
                                       lat_effect = 0.1,
                                       cati_effect = 0.1,
                                       orc_capacity = 500,
                                       ctc_capacity = 100,
                                       vax_regimen = c("none", "1dose", "2dose"),
                                       vax_total_doses = c(0, 280000),
                                       vax_delay = c(0, 14),
                                       vax_campaign_days = 6,
                                       vax_dose_interval = 14,
                                       prefix = "cf") {
  if (!is.numeric(trigger_time) || length(trigger_time) != 1) stop("trigger_time must be a single number", call. = FALSE)
  if (!is.numeric(horizon) || length(horizon) != 1) stop("horizon must be a single number", call. = FALSE)

  vax_regimen <- as.character(vax_regimen)
  ok <- vax_regimen %in% c("none", "1dose", "2dose")
  if (!all(ok)) stop("vax_regimen must be among: none, 1dose, 2dose", call. = FALSE)

  grid <- expand.grid(
    aa_start_offset = aa_start_offset,
    wash_duration = wash_duration,
    care_duration = care_duration,
    chlor_effect = chlor_effect,
    hyg_effect = hyg_effect,
    lat_effect = lat_effect,
    cati_effect = cati_effect,
    orc_capacity = orc_capacity,
    ctc_capacity = ctc_capacity,
    vax_regimen = vax_regimen,
    vax_total_doses = vax_total_doses,
    vax_delay = vax_delay,
    stringsAsFactors = FALSE
  )

  modify <- vector("list", nrow(grid))
  scenario <- character(nrow(grid))

  for (i in seq_len(nrow(grid))) {
    r <- as.list(grid[i, , drop = FALSE])

    aa_start <- trigger_time + as.numeric(r$aa_start_offset)

    w_wash <- .chlaa_window(aa_start, as.numeric(r$wash_duration), horizon = horizon)
    w_care <- .chlaa_window(aa_start, as.numeric(r$care_duration), horizon = horizon)

    m <- list(
      chlor_start = w_wash$start, chlor_end = w_wash$end, chlor_effect = as.numeric(r$chlor_effect),
      hyg_start = w_wash$start, hyg_end = w_wash$end, hyg_effect = as.numeric(r$hyg_effect),
      lat_start = w_wash$start, lat_end = w_wash$end, lat_effect = as.numeric(r$lat_effect),
      cati_start = w_wash$start, cati_end = w_wash$end, cati_effect = as.numeric(r$cati_effect),

      orc_start = w_care$start, orc_end = w_care$end, orc_capacity = as.numeric(r$orc_capacity),
      ctc_start = w_care$start, ctc_end = w_care$end, ctc_capacity = as.numeric(r$ctc_capacity)
    )

    reg <- as.character(r$vax_regimen)
    doses <- as.numeric(r$vax_total_doses)
    delay <- as.numeric(r$vax_delay)

    if (reg == "none" || doses <= 0) {
      m <- c(m, list(
        vax1_start = 0, vax1_end = 0, vax1_total_doses = 0, vax1_doses_per_day = 0,
        vax2_start = 0, vax2_end = 0, vax2_total_doses = 0, vax2_doses_per_day = 0
      ))
    } else if (reg == "1dose") {
      plan <- .chlaa_make_vax_plan(
        N = 1, trigger_time = aa_start, total_doses = doses,
        regimen = "1dose", delay = delay, campaign_days = vax_campaign_days,
        dose_interval = vax_dose_interval, horizon = horizon
      )
      m <- c(m, plan)
    } else if (reg == "2dose") {
      plan <- .chlaa_make_vax_plan(
        N = 1, trigger_time = aa_start, total_doses = doses,
        regimen = "2dose", delay = delay, campaign_days = vax_campaign_days,
        dose_interval = vax_dose_interval, horizon = horizon
      )
      m <- c(m, plan)
    }

    modify[[i]] <- m

    name_fields <- list(
      off = r$aa_start_offset,
      v = r$vax_regimen,
      d = r$vax_total_doses,
      del = r$vax_delay,
      ce = r$chlor_effect,
      orc = r$orc_capacity,
      ctc = r$ctc_capacity
    )
    scenario[[i]] <- .chlaa_make_name(prefix, name_fields)
  }

  res <- grid
  res$scenario <- scenario
  res$modify <- I(modify)

  if (requireNamespace("tibble", quietly = TRUE)) {
    res <- tibble::as_tibble(res)
  }
  res
}

#' Convert a counterfactual grid to scenario objects
#'
#' @param grid Output from `chlaa_counterfactual_grid()`.
#'
#' @return A list of `chlaa_scenario` objects.
#' @export
chlaa_scenarios_from_grid <- function(grid) {
  if (!is.data.frame(grid)) stop("grid must be a data.frame", call. = FALSE)
  if (!all(c("scenario", "modify") %in% names(grid))) {
    stop("grid must contain columns scenario and modify", call. = FALSE)
  }

  out <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    out[[i]] <- chlaa_scenario(grid$scenario[[i]], grid$modify[[i]])
  }
  out
}

#' Run a counterfactual grid (convenience)
#'
#' @param pars Baseline parameter list.
#' @param grid Counterfactual grid from `chlaa_counterfactual_grid()`.
#' @param time Simulation times.
#' @param baseline_name Name of baseline scenario to include.
#' @param n_particles Particles per scenario.
#' @param dt Time step.
#' @param seed Seed.
#'
#' @return Output of `chlaa_run_scenarios()` (baseline + all grid scenarios).
#' @export
chlaa_run_counterfactual_grid <- function(pars,
                                           grid,
                                           time,
                                           baseline_name = "baseline",
                                           n_particles = 200,
                                           dt = 0.25,
                                           seed = 1) {
  .check_named_list(pars, "pars")
  sc_base <- chlaa_scenario(baseline_name, list())
  sc_grid <- chlaa_scenarios_from_grid(grid)
  chlaa_run_scenarios(
    pars = pars,
    scenarios = c(list(sc_base), sc_grid),
    time = time,
    n_particles = n_particles,
    dt = dt,
    seed = seed
  )
}

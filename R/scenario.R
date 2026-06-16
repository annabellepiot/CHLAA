# Scenario running and comparison

#' Create a counterfactual scenario
#'
#' @param name Scenario name.
#' @param modify Named list of parameter modifications.
#'
#' @return An object of class `chlaa_scenario`.
#' @export
chlaa_scenario <- function(name, modify) {
  if (!is.character(name) || length(name) != 1) stop("name must be a single string", call. = FALSE)
  if (!is.list(modify)) stop("modify must be a list", call. = FALSE)
  if (length(modify) > 0) .check_named_list(modify, "modify")
  structure(list(name = name, modify = modify), class = "chlaa_scenario")
}

#' Run multiple scenarios
#'
#' @param pars Baseline parameters.
#' @param scenarios List of `chlaa_scenario` objects.
#' @param time Simulation times.
#' @param n_particles Particles per scenario.
#' @param dt Time step.
#' @param seed Base seed (scenario i uses seed + i).
#'
#' @return A data.frame of simulation outputs with a `scenario` column. The returned data frame
#'   has an attribute `scenario_parameters`, a named list of the parameter lists used per scenario.
#' @export
chlaa_run_scenarios <- function(pars,
                                 scenarios,
                                 time,
                                 n_particles = 200,
                                 dt = 0.25,
                                 seed = 1) {
  if (!is.list(scenarios) || length(scenarios) == 0) {
    stop("scenarios must be a non-empty list of chlaa_scenario objects", call. = FALSE)
  }

  pars_used <- list()

  out <- lapply(seq_along(scenarios), function(i) {
    sc <- scenarios[[i]]
    if (!inherits(sc, "chlaa_scenario")) stop("All scenarios must be chlaa_scenario objects", call. = FALSE)

    p <- utils::modifyList(pars, sc$modify)
    chlaa_parameters_validate(p)
    pars_used[[sc$name]] <<- p

    sim <- chlaa_simulate(p, time = time, n_particles = n_particles, dt = dt, seed = seed + i)
    sim$scenario <- sc$name
    sim
  })

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    res <- do.call(rbind, out)
    attr(res, "scenario_parameters") <- pars_used
    return(res)
  }

  res <- dplyr::bind_rows(out) |>
    dplyr::relocate("scenario", "time", "particle") |>
    dplyr::arrange(.data$scenario, .data$time, .data$particle)

  attr(res, "scenario_parameters") <- pars_used
  res
}

.chlaa_end_of_horizon <- function(scenario_runs) {
  req <- c(
    "scenario","time","particle",
    "cum_infections","cum_symptoms","cum_deaths",
    "cum_vax1","cum_vax2",
    "cum_orc_treated","cum_ctc_treated"
  )
  missing <- setdiff(req, names(scenario_runs))
  if (length(missing) > 0) {
    stop("scenario_runs missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    # base fallback
    end_time <- max(scenario_runs$time)
    end <- scenario_runs[scenario_runs$time == end_time, , drop = FALSE]
    end$doses_total <- end$cum_vax1 + end$cum_vax2
    return(end)
  }

  scenario_runs |>
    dplyr::group_by(.data$scenario, .data$particle) |>
    dplyr::filter(.data$time == max(.data$time)) |>
    dplyr::ungroup() |>
    dplyr::mutate(doses_total = .data$cum_vax1 + .data$cum_vax2)
}

#' Compare scenario outcomes against a baseline
#'
#' Produces a scenario comparison table including infections, symptomatic cases, deaths, doses,
#' treatment counts, and (optionally) costs and DALYs with ICERs vs baseline.
#'
#' @param scenario_runs Output from `chlaa_run_scenarios()`.
#' @param baseline Baseline scenario name.
#' @param include_econ If TRUE, include health economic outputs (requires `chlaa_health_econ()`).
#' @param econ Optional named list of economic parameter overrides (passed to `chlaa_health_econ()`).
#' @param wtp Optional willingness-to-pay per DALY averted for NMB/INMB outputs.
#'
#' @return A data.frame of scenario-level means and deltas vs baseline.
#' @export
chlaa_compare_scenarios <- function(scenario_runs,
                                      baseline,
                                      include_econ = TRUE,
                                      econ = NULL,
                                      wtp = NULL) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr is required for chlaa_compare_scenarios()", call. = FALSE)
  }

  end <- .chlaa_end_of_horizon(scenario_runs)

  summ <- end |>
    dplyr::group_by(.data$scenario) |>
    dplyr::summarise(
      infections = mean(.data$cum_infections),
      cases_symptomatic = mean(.data$cum_symptoms),
      deaths = mean(.data$cum_deaths),
      doses = mean(.data$doses_total),
      orc_treated = mean(.data$cum_orc_treated),
      ctc_treated = mean(.data$cum_ctc_treated),
      .groups = "drop"
    )

  if (!baseline %in% summ$scenario) stop("baseline scenario not found", call. = FALSE)
  base <- summ[summ$scenario == baseline, , drop = FALSE]

  summ <- summ |>
    dplyr::mutate(
      infections_averted = base$infections - .data$infections,
      cases_averted = base$cases_symptomatic - .data$cases_symptomatic,
      deaths_averted = base$deaths - .data$deaths,
      infections_diff = .data$infections - base$infections,
      cases_diff = .data$cases_symptomatic - base$cases_symptomatic,
      deaths_diff = .data$deaths - base$deaths
    )

  if (!isTRUE(include_econ)) return(summ)

  he <- chlaa_health_econ(scenario_runs, econ = econ)
  out <- dplyr::left_join(summ, he, by = "scenario")

  base2 <- out[out$scenario == baseline, , drop = FALSE]

  if (!is.null(wtp)) {
    if (!is.numeric(wtp) || length(wtp) != 1 || wtp < 0) {
      stop("wtp must be a single non-negative number", call. = FALSE)
    }
    out <- out |>
      dplyr::mutate(
        nmb = wtp * (base2$mean_dalys - .data$mean_dalys) - (.data$mean_cost - base2$mean_cost),
        inmb = .data$nmb
      )
  }
  res <- out |>
    dplyr::mutate(
      cost = .data$mean_cost,
      dalys = .data$mean_dalys,
      cost_diff = .data$cost - base2$mean_cost,
      dalys_averted = base2$mean_dalys - .data$dalys,
      icer_cost_per_daly_averted = dplyr::if_else(
        .data$dalys_averted > 0,
        .data$cost_diff / .data$dalys_averted,
        NA_real_
      ),
      icer_cost_per_death_averted = dplyr::if_else(
        .data$deaths_averted > 0,
        .data$cost_diff / .data$deaths_averted,
        NA_real_
      )
    ) |>
    dplyr::select(
      "scenario",
      "infections", "cases_symptomatic", "deaths",
      "doses", "orc_treated", "ctc_treated",
      "infections_averted", "cases_averted", "deaths_averted",
      "cost", "dalys",
      "cost_diff", "dalys_averted",
      "icer_cost_per_daly_averted", "icer_cost_per_death_averted",
      "mean_cost_vax", "mean_cost_care", "mean_cost_wash"
    )

  if (!is.null(wtp)) {
    res <- dplyr::mutate(res, nmb = out$nmb, inmb = out$inmb)
  }
  res
}

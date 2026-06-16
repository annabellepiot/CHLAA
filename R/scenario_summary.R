# Scenario templates and reporting helpers

#' Build standard operational scenario templates
#'
#' @param pars Baseline parameters.
#' @param trigger_time Trigger time for anticipatory interventions.
#' @param horizon End of horizon (used to cap intervention windows).
#' @param vax_total_doses Total vaccine doses for vaccination scenarios.
#' @param early_offset Days relative to trigger for early vaccination start.
#' @param late_offset Days relative to trigger for late vaccination start.
#' @param campaign_days Vaccination campaign duration.
#' @param wash_duration WASH intervention duration.
#' @param care_duration Care intervention duration.
#'
#' @return A list of `chlaa_scenario` objects.
#' @export
chlaa_standard_scenarios <- function(pars,
                                       trigger_time,
                                       horizon,
                                       vax_total_doses = 280000,
                                       early_offset = -14,
                                       late_offset = 14,
                                       campaign_days = 6,
                                       wash_duration = 120,
                                       care_duration = 120) {
  .check_named_list(pars, "pars")
  if (!is.numeric(trigger_time) || length(trigger_time) != 1) stop("trigger_time must be a single number", call. = FALSE)
  if (!is.numeric(horizon) || length(horizon) != 1 || horizon <= trigger_time) {
    stop("horizon must be > trigger_time", call. = FALSE)
  }

  sc <- list(
    chlaa_scenario("baseline", list())
  )

  no_vax <- chlaa_make_aa_scenarios(
    pars = pars,
    trigger_time = trigger_time,
    horizon = horizon,
    include_no_intervention = FALSE,
    include_aa_no_vax = TRUE,
    vax_total_doses = 0,
    wash_duration = wash_duration,
    care_duration = care_duration
  )
  no_vax <- no_vax[vapply(no_vax, function(x) x$name == "aa_no_vax", logical(1))]
  if (length(no_vax) == 1) {
    no_vax[[1]]$name <- "no_vaccination"
    sc <- c(sc, no_vax)
  }

  plan_early <- .chlaa_make_vax_plan(
    N = pars$N,
    trigger_time = trigger_time,
    total_doses = vax_total_doses,
    regimen = "1dose",
    delay = early_offset,
    campaign_days = campaign_days,
    horizon = horizon
  )
  plan_late <- .chlaa_make_vax_plan(
    N = pars$N,
    trigger_time = trigger_time,
    total_doses = vax_total_doses,
    regimen = "1dose",
    delay = late_offset,
    campaign_days = campaign_days,
    horizon = horizon
  )

  w_wash <- .chlaa_window(trigger_time, wash_duration, horizon = horizon)
  w_care <- .chlaa_window(trigger_time, care_duration, horizon = horizon)

  wash_mod <- list(
    chlor_start = w_wash$start, chlor_end = w_wash$end, chlor_effect = 0.2,
    hyg_start = w_wash$start, hyg_end = w_wash$end, hyg_effect = 0.2,
    lat_start = w_wash$start, lat_end = w_wash$end, lat_effect = 0.1,
    cati_start = w_wash$start, cati_end = w_wash$end, cati_effect = 0.1,
    orc_start = 0, orc_end = 0,
    ctc_start = 0, ctc_end = 0
  )

  care_mod <- list(
    orc_start = w_care$start, orc_end = w_care$end, orc_capacity = pars$orc_capacity,
    ctc_start = w_care$start, ctc_end = w_care$end, ctc_capacity = pars$ctc_capacity,
    chlor_start = 0, chlor_end = 0, chlor_effect = 0,
    hyg_start = 0, hyg_end = 0, hyg_effect = 0,
    lat_start = 0, lat_end = 0, lat_effect = 0,
    cati_start = 0, cati_end = 0, cati_effect = 0
  )

  combined_mod <- utils::modifyList(wash_mod, care_mod)
  combined_mod <- utils::modifyList(combined_mod, plan_early)

  early_mod <- utils::modifyList(combined_mod, plan_early)
  late_mod <- utils::modifyList(combined_mod, plan_late)

  sc <- c(
    sc,
    list(
      chlaa_scenario("early_vaccination", early_mod),
      chlaa_scenario("late_vaccination", late_mod),
      chlaa_scenario("wash_only", wash_mod),
      chlaa_scenario("cm_only", care_mod),
      chlaa_scenario("combined_package", combined_mod)
    )
  )

  sc
}

#' Summarise scenario outcomes with uncertainty and baseline deltas
#'
#' @param scenario_runs Output from `chlaa_run_scenarios()`.
#' @param baseline Baseline scenario name.
#' @param incidence_var Incidence variable used for peak and control metrics.
#' @param control_threshold Threshold used for time-to-control metrics.
#'
#' @return A data.frame with scenario summaries.
#' @export
chlaa_scenario_summary <- function(scenario_runs,
                                     baseline = "baseline",
                                     incidence_var = "inc_symptoms",
                                     control_threshold = 1) {
  .require_suggested("dplyr")
  req <- c("scenario", "time", "particle", "cum_symptoms", "cum_deaths", incidence_var)
  miss <- setdiff(req, names(scenario_runs))
  if (length(miss) > 0) stop("scenario_runs missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

  end <- scenario_runs |>
    dplyr::group_by(.data$scenario, .data$particle) |>
    dplyr::filter(.data$time == max(.data$time)) |>
    dplyr::ungroup()

  peak <- scenario_runs |>
    dplyr::group_by(.data$scenario, .data$particle) |>
    dplyr::summarise(
      peak_incidence = max(.data[[incidence_var]], na.rm = TRUE),
      time_peak = .data$time[which.max(.data[[incidence_var]])][1],
      .groups = "drop"
    )

  ttc <- scenario_runs |>
    dplyr::group_by(.data$scenario, .data$particle) |>
    dplyr::summarise(
      time_to_control = {
        idx <- which(.data[[incidence_var]] <= control_threshold)
        if (length(idx) == 0) NA_real_ else .data$time[min(idx)]
      },
      .groups = "drop"
    )

  by_particle <- end |>
    dplyr::select("scenario", "particle", "cum_symptoms", "cum_deaths") |>
    dplyr::left_join(peak, by = c("scenario", "particle")) |>
    dplyr::left_join(ttc, by = c("scenario", "particle"))

  base <- by_particle[by_particle$scenario == baseline, , drop = FALSE]
  if (nrow(base) == 0) stop("baseline scenario not found", call. = FALSE)
  base <- base[order(base$particle), , drop = FALSE]

  merged <- by_particle |>
    dplyr::left_join(
      base |>
        dplyr::rename(
          base_cases = "cum_symptoms",
          base_deaths = "cum_deaths"
        ) |>
        dplyr::select("particle", "base_cases", "base_deaths"),
      by = "particle"
    ) |>
    dplyr::mutate(
      cases_averted = .data$base_cases - .data$cum_symptoms,
      deaths_averted = .data$base_deaths - .data$cum_deaths
    )

  merged |>
    dplyr::group_by(.data$scenario) |>
    dplyr::summarise(
      total_cases = mean(.data$cum_symptoms),
      total_deaths = mean(.data$cum_deaths),
      cases_averted = mean(.data$cases_averted),
      deaths_averted = mean(.data$deaths_averted),
      peak_incidence = mean(.data$peak_incidence),
      time_peak = mean(.data$time_peak),
      time_to_control = mean(.data$time_to_control, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Produce a compact scenario analysis report object
#'
#' @param scenario_runs Output from `chlaa_run_scenarios()`.
#' @param baseline Baseline scenario name.
#' @param include_econ Include economic comparison table.
#' @param econ Economic parameter overrides passed to economics functions.
#' @param wtp Optional willingness-to-pay per DALY for NMB outputs.
#'
#' @return A list with summary table and key plots.
#' @export
chlaa_scenario_report <- function(scenario_runs,
                                    baseline = "baseline",
                                    include_econ = TRUE,
                                    econ = NULL,
                                    wtp = NULL) {
  summary_tbl <- chlaa_scenario_summary(scenario_runs, baseline = baseline)
  cmp <- chlaa_compare_scenarios(
    scenario_runs = scenario_runs,
    baseline = baseline,
    include_econ = include_econ,
    econ = econ,
    wtp = wtp
  )

  plots <- list(
    incidence = chlaa_plot_incidence(scenario_runs, var = "inc_symptoms"),
    cumulative_deaths = chlaa_plot_incidence(scenario_runs, var = "cum_deaths"),
    scenario_comparison = chlaa_plot_scenarios(cmp, metric = "deaths")
  )

  if (isTRUE(include_econ)) {
    plots$ce_plane <- chlaa_plot_ce_plane(cmp)
  }

  list(
    summary = summary_tbl,
    comparison = cmp,
    plots = plots
  )
}

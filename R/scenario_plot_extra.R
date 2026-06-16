# Additional scenario plotting helpers

#' Plot Scenario Trajectories On A Shared Axis
#'
#' Summarises each scenario over particles and overlays median and uncertainty
#' ribbons in one panel to make cross-scenario differences easier to inspect.
#'
#' @param scenario_runs Output from `chlaa_run_scenarios()` or
#'   `chlaa_run_scenarios_from_snapshot()`.
#' @param var Variable to summarise (default `inc_symptoms`).
#' @param probs Two-element interval for ribbon bounds.
#'
#' @return A ggplot object.
#' @export
chlaa_plot_scenario_overlay <- function(scenario_runs,
                                        var = "inc_symptoms",
                                        probs = c(0.1, 0.9)) {
  .require_suggested("ggplot2")
  .require_suggested("dplyr")

  if (!is.data.frame(scenario_runs)) stop("scenario_runs must be a data.frame", call. = FALSE)
  req <- c("scenario", "time", "particle", var)
  missing <- setdiff(req, names(scenario_runs))
  if (length(missing) > 0) stop("scenario_runs missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!is.numeric(probs) || length(probs) != 2 || any(probs <= 0) || any(probs >= 1) || probs[1] >= probs[2]) {
    stop("probs must be length-2 with 0 < probs[1] < probs[2] < 1", call. = FALSE)
  }

  lo <- probs[1]
  hi <- probs[2]

  summ <- scenario_runs |>
    dplyr::group_by(.data$scenario, .data$time) |>
    dplyr::summarise(
      q_lo = stats::quantile(.data[[var]], lo, na.rm = TRUE),
      q50 = stats::quantile(.data[[var]], 0.5, na.rm = TRUE),
      q_hi = stats::quantile(.data[[var]], hi, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot2::ggplot(summ, ggplot2::aes(x = .data$time, colour = .data$scenario, fill = .data$scenario)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$q_lo, ymax = .data$q_hi), alpha = 0.15, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(y = .data$q50), linewidth = 0.8) +
    ggplot2::labs(
      x = "Time (days)",
      y = var,
      title = paste0("Scenario overlay: ", var),
      subtitle = paste0("Ribbon = ", round(100 * lo), "-", round(100 * hi), "% interval")
    ) +
    ggplot2::theme_minimal()
}

#' Plot Scenario Differences Versus Baseline Through Time
#'
#' Computes `scenario - baseline` differences at each time and shows uncertainty
#' ribbons by scenario. Optionally plots cumulative differences.
#'
#' @param scenario_runs Output from `chlaa_run_scenarios()` or
#'   `chlaa_run_scenarios_from_snapshot()`.
#' @param baseline Baseline scenario name.
#' @param var Variable to compare (default `inc_symptoms`).
#' @param cumulative If TRUE, show cumulative differences over time.
#' @param probs Two-element interval for ribbon bounds.
#'
#' @return A ggplot object.
#' @export
chlaa_plot_difference_vs_baseline <- function(scenario_runs,
                                              baseline = "baseline",
                                              var = "inc_symptoms",
                                              cumulative = FALSE,
                                              probs = c(0.1, 0.9)) {
  .require_suggested("ggplot2")
  .require_suggested("dplyr")

  if (!is.data.frame(scenario_runs)) stop("scenario_runs must be a data.frame", call. = FALSE)
  req <- c("scenario", "time", "particle", var)
  missing <- setdiff(req, names(scenario_runs))
  if (length(missing) > 0) stop("scenario_runs missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!baseline %in% scenario_runs$scenario) stop("baseline scenario not found", call. = FALSE)

  lo <- probs[1]
  hi <- probs[2]

  base <- scenario_runs |>
    dplyr::filter(.data$scenario == baseline) |>
    dplyr::select("time", "particle", base_val = dplyr::all_of(var))

  cmp <- scenario_runs |>
    dplyr::filter(.data$scenario != baseline) |>
    dplyr::left_join(base, by = c("time", "particle")) |>
    dplyr::mutate(diff = .data[[var]] - .data$base_val)

  if (isTRUE(cumulative)) {
    cmp <- cmp |>
      dplyr::arrange(.data$scenario, .data$particle, .data$time) |>
      dplyr::group_by(.data$scenario, .data$particle) |>
      dplyr::mutate(diff = cumsum(.data$diff)) |>
      dplyr::ungroup()
  }

  summ <- cmp |>
    dplyr::group_by(.data$scenario, .data$time) |>
    dplyr::summarise(
      q_lo = stats::quantile(.data$diff, lo, na.rm = TRUE),
      q50 = stats::quantile(.data$diff, 0.5, na.rm = TRUE),
      q_hi = stats::quantile(.data$diff, hi, na.rm = TRUE),
      .groups = "drop"
    )

  ylab <- if (isTRUE(cumulative)) {
    paste0("Cumulative ", var, " difference vs ", baseline)
  } else {
    paste0(var, " difference vs ", baseline)
  }

  ggplot2::ggplot(summ, ggplot2::aes(x = .data$time, colour = .data$scenario, fill = .data$scenario)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$q_lo, ymax = .data$q_hi), alpha = 0.15, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(y = .data$q50), linewidth = 0.8) +
    ggplot2::labs(
      x = "Time (days)",
      y = ylab,
      title = paste0("Scenario differences vs baseline: ", var)
    ) +
    ggplot2::theme_minimal()
}

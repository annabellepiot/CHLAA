# Plotting helpers

#' Plot trajectories of selected model variables
#'
#' @param sim Output from `chlaa_simulate()` or `chlaa_run_scenarios()`.
#' @param vars Variables to plot.
#' @param facet If TRUE, facet by variable.
#'
#' @return A ggplot object.
#' @export
chlaa_plot_trajectories <- function(sim, vars = c("inc_symptoms", "cum_deaths"), facet = TRUE) {
  .require_suggested("ggplot2")
  .require_suggested("tidyr")
  .require_suggested("dplyr")

  if (!is.data.frame(sim) || !all(c("time", "particle") %in% names(sim))) {
    stop("sim must be a data.frame with columns time and particle", call. = FALSE)
  }

  missing <- setdiff(vars, names(sim))
  if (length(missing) > 0) stop("Missing variables: ", paste(missing, collapse = ", "), call. = FALSE)

  df <- sim |>
    dplyr::select(dplyr::any_of(c("scenario", "time", "particle", vars))) |>
    tidyr::pivot_longer(cols = dplyr::all_of(vars), names_to = "variable", values_to = "value")

  if (!"scenario" %in% names(df)) df$scenario <- "simulation"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$value, group = .data$particle)) +
    ggplot2::geom_line(alpha = 0.15) +
    ggplot2::facet_wrap(~ .data$scenario, scales = "free_y") +
    ggplot2::labs(x = "Time (days)", y = NULL, title = "Model trajectories") +
    ggplot2::theme_minimal()

  if (isTRUE(facet)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$value, group = .data$particle)) +
      ggplot2::geom_line(alpha = 0.15) +
      ggplot2::facet_grid(.data$variable ~ .data$scenario, scales = "free_y") +
      ggplot2::labs(x = "Time (days)", y = NULL, title = "Model trajectories") +
      ggplot2::theme_minimal()
  }

  p
}

#' Plot daily incidence summary for a variable
#'
#' @param sim Simulation output.
#' @param var Incidence variable (default inc_symptoms).
#'
#' @return A ggplot object.
#' @export
chlaa_plot_incidence <- function(sim, var = "inc_symptoms") {
  .require_suggested("ggplot2")
  .require_suggested("dplyr")

  if (!is.data.frame(sim) || !all(c("time", "particle", var) %in% names(sim))) {
    stop("sim must be a data.frame with columns time, particle, and ", var, call. = FALSE)
  }

  df <- sim
  if (!"scenario" %in% names(df)) df$scenario <- "simulation"

  summ <- df |>
    dplyr::group_by(.data$scenario, .data$time) |>
    dplyr::summarise(
      mean = mean(.data[[var]]),
      q025 = stats::quantile(.data[[var]], 0.025),
      q975 = stats::quantile(.data[[var]], 0.975),
      .groups = "drop"
    )

  ggplot2::ggplot(summ, ggplot2::aes(x = .data$time, y = .data$mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$q025, ymax = .data$q975), alpha = 0.2) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ .data$scenario, scales = "free_y") +
    ggplot2::labs(x = "Time (days)", y = var, title = paste("Incidence:", var)) +
    ggplot2::theme_minimal()
}

#' Plot scenario comparison bars for a chosen metric
#'
#' @param cmp Output from `chlaa_compare_scenarios()`.
#' @param metric Column in `cmp` to plot.
#'
#' @return A ggplot object.
#' @export
chlaa_plot_scenarios <- function(cmp, metric = "deaths") {
  .require_suggested("ggplot2")
  if (!is.data.frame(cmp) || !"scenario" %in% names(cmp) || !metric %in% names(cmp)) {
    stop("cmp must be a data.frame with columns scenario and ", metric, call. = FALSE)
  }

  ggplot2::ggplot(cmp, ggplot2::aes(x = .data$scenario, y = .data[[metric]])) +
    ggplot2::geom_col() +
    ggplot2::labs(x = NULL, y = metric, title = paste("Scenario comparison:", metric)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot a cost-effectiveness plane from scenario comparisons
#'
#' Expects output from `chlaa_compare_scenarios(..., include_econ = TRUE)`.
#'
#' @param cmp Comparison table.
#' @return A ggplot object.
#' @export
chlaa_plot_ce_plane <- function(cmp) {
  .require_suggested("ggplot2")
  req <- c("scenario", "cost_diff", "dalys_averted")
  if (!all(req %in% names(cmp))) {
    stop("cmp must contain columns: ", paste(req, collapse = ", "), call. = FALSE)
  }

  ggplot2::ggplot(cmp, ggplot2::aes(x = .data$dalys_averted, y = .data$cost_diff, label = .data$scenario)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_point() +
    ggplot2::geom_text(vjust = -0.4, check_overlap = TRUE) +
    ggplot2::labs(
      x = "DALYs averted vs baseline",
      y = "Incremental cost vs baseline",
      title = "Cost-effectiveness plane"
    ) +
    ggplot2::theme_minimal()
}

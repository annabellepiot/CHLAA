# Plotting for scenario forecast summaries

.parse_qcols <- function(df) {
  qcols <- grep("^q", names(df), value = TRUE)
  if (length(qcols) == 0) return(list(qcols = character(0), probs = numeric(0)))

  probs <- suppressWarnings(as.numeric(gsub("p", ".", sub("^q", "", qcols))))
  ok <- is.finite(probs)
  list(qcols = qcols[ok], probs = probs[ok])
}

.closest_qcol <- function(qcols, probs, target) {
  if (length(qcols) == 0) return(NULL)
  qcols[which.min(abs(probs - target))]
}

#' Plot scenario forecast summaries
#'
#' Works with output from `chlaa_forecast_scenarios_from_fit()`.
#'
#' @param forecast Scenario forecast table.
#' @param var Variable to plot (e.g. "cases", "inc_symptoms", "cum_deaths").
#' @param type "absolute" for scenario trajectories, or "difference" for scenario minus baseline.
#' @param scenarios Optional character vector of scenario names to include.
#' @param include_baseline Include baseline scenario in the plot.
#' @param facet If TRUE, facet by scenario; if FALSE, overlay scenarios.
#' @param show_mean If TRUE, add a dashed mean line.
#' @param data Optional observed data to overlay (only meaningful for type = "absolute").
#' @param data_time Time column name in `data`.
#' @param data_y Y column name in `data` (defaults to `var`).
#'
#' @return A ggplot object.
#' @export
chlaa_plot_scenario_forecasts <- function(forecast,
                                            var,
                                            type = c("absolute", "difference"),
                                            scenarios = NULL,
                                            include_baseline = NULL,
                                            facet = FALSE,
                                            show_mean = FALSE,
                                            data = NULL,
                                            data_time = "time",
                                            data_y = NULL) {
  .require_suggested("ggplot2")

  type <- match.arg(type)

  if (!is.data.frame(forecast)) stop("forecast must be a data.frame", call. = FALSE)
  req <- c("scenario", "type", "time", "variable", "mean")
  if (!all(req %in% names(forecast))) {
    stop("forecast must contain columns: ", paste(req, collapse = ", "), call. = FALSE)
  }
  if (!is.character(var) || length(var) != 1) stop("var must be a single string", call. = FALSE)

  baseline <- attr(forecast, "baseline_name")
  if (is.null(baseline) || !is.character(baseline) || length(baseline) != 1) baseline <- "baseline"

  if (is.null(include_baseline)) include_baseline <- (type == "absolute")

  df <- forecast[forecast$variable == var & forecast$type == type, , drop = FALSE]
  if (nrow(df) == 0) stop("No rows for var = ", var, " and type = ", type, call. = FALSE)

  if (!is.null(scenarios)) {
    df <- df[df$scenario %in% scenarios, , drop = FALSE]
  }
  if (!isTRUE(include_baseline)) {
    df <- df[df$scenario != baseline, , drop = FALSE]
  }

  qp <- .parse_qcols(df)
  qcols <- qp$qcols
  probs <- qp$probs
  if (length(qcols) == 0) stop("No quantile columns found (expected columns starting with 'q')", call. = FALSE)

  lo_outer <- qcols[which.min(probs)]
  hi_outer <- qcols[which.max(probs)]

  ord <- order(probs)
  qcols_ord <- qcols[ord]
  probs_ord <- probs[ord]

  lo_inner <- .closest_qcol(qcols, probs, 0.25)
  hi_inner <- .closest_qcol(qcols, probs, 0.75)
  if (is.null(lo_inner) || is.null(hi_inner)) {
    if (length(qcols_ord) >= 4) {
      lo_inner <- qcols_ord[2]
      hi_inner <- qcols_ord[length(qcols_ord) - 1]
    } else {
      lo_inner <- lo_outer
      hi_inner <- hi_outer
    }
  }

  med_col <- .closest_qcol(qcols, probs, 0.5)
  use_median <- !is.null(med_col)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time))

  if (isTRUE(facet)) {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[lo_outer]], ymax = .data[[hi_outer]]), alpha = 0.2) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[lo_inner]], ymax = .data[[hi_inner]]), alpha = 0.3)

    if (use_median) {
      p <- p + ggplot2::geom_line(ggplot2::aes(y = .data[[med_col]]))
    } else {
      p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$mean))
    }

    if (isTRUE(show_mean) && use_median) {
      p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$mean), linetype = 2)
    }

    p <- p + ggplot2::facet_wrap(~ .data$scenario, scales = "free_y")
  } else {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data[[lo_outer]], ymax = .data[[hi_outer]], fill = .data$scenario),
        alpha = 0.12
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data[[lo_inner]], ymax = .data[[hi_inner]], fill = .data$scenario),
        alpha = 0.22
      )

    if (use_median) {
      p <- p + ggplot2::geom_line(ggplot2::aes(y = .data[[med_col]], colour = .data$scenario))
    } else {
      p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$mean, colour = .data$scenario))
    }

    if (isTRUE(show_mean) && use_median) {
      p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$mean, colour = .data$scenario), linetype = 2)
    }
  }

  if (type == "difference") {
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = 2)
  }

  if (!is.null(data) && type == "absolute") {
    if (!is.data.frame(data)) stop("data must be a data.frame", call. = FALSE)
    if (!data_time %in% names(data)) stop("data_time column not found in data", call. = FALSE)
    if (is.null(data_y)) data_y <- var
    if (!data_y %in% names(data)) stop("data_y column not found in data", call. = FALSE)

    p <- p + ggplot2::geom_point(
      data = data,
      ggplot2::aes(x = .data[[data_time]], y = .data[[data_y]]),
      inherit.aes = FALSE
    )
  }

  title_txt <- if (type == "absolute") paste("Scenario forecasts:", var) else paste("Scenario difference vs baseline:", var)

  p +
    ggplot2::labs(x = "Time (days)", y = var, title = title_txt) +
    ggplot2::theme_minimal()
}

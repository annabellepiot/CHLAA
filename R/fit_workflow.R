# Fit workflow helpers

#' Prepare and validate incidence data for model fitting
#'
#' @param data Input data.frame.
#' @param time_col Name of time column.
#' @param cases_col Name of observed case-count column.
#' @param expected_step Optional expected time step. If provided, times must follow this step.
#' @param fill_missing If TRUE and `expected_step` is provided, fill missing time rows with zero cases.
#'
#' @return A data.frame with columns `time` and `cases`.
#' @export
chlaa_prepare_data <- function(data,
                                 time_col = "time",
                                 cases_col = "cases",
                                 expected_step = NULL,
                                 fill_missing = FALSE) {
  if (!is.data.frame(data)) stop("data must be a data.frame", call. = FALSE)
  if (!time_col %in% names(data)) stop("time_col not found in data", call. = FALSE)
  if (!cases_col %in% names(data)) stop("cases_col not found in data", call. = FALSE)

  out <- data.frame(
    time = as.numeric(data[[time_col]]),
    cases = data[[cases_col]],
    stringsAsFactors = FALSE
  )

  # Pass through deaths column if present
  if ("deaths" %in% names(data)) {
    out$deaths <- data[["deaths"]]
    out$deaths[is.na(out$deaths)] <- 0L
    out$deaths <- as.integer(round(out$deaths))
  }

  if (any(!is.finite(out$time))) stop("time values must be finite numeric", call. = FALSE)
  if (any(is.na(out$cases))) stop("cases contains missing values", call. = FALSE)
  if (any(out$cases < 0)) stop("cases must be non-negative", call. = FALSE)
  if (any(abs(out$cases - round(out$cases)) > sqrt(.Machine$double.eps))) {
    stop("cases must be integer-like counts", call. = FALSE)
  }
  out$cases <- as.integer(round(out$cases))

  o <- order(out$time)
  out <- out[o, , drop = FALSE]

  d <- diff(out$time)
  if (any(d <= 0)) stop("time values must be strictly increasing", call. = FALSE)

  if (!is.null(expected_step)) {
    if (!is.numeric(expected_step) || length(expected_step) != 1 || expected_step <= 0) {
      stop("expected_step must be a single positive number", call. = FALSE)
    }

    ok_step <- abs(d - expected_step) <= sqrt(.Machine$double.eps) * pmax(1, abs(expected_step))
    if (all(ok_step)) return(out)

    if (!isTRUE(fill_missing)) {
      stop("time step is not consistent with expected_step; set fill_missing = TRUE to pad missing rows", call. = FALSE)
    }

    full_time <- seq(min(out$time), max(out$time), by = expected_step)
    merged <- merge(
      data.frame(time = full_time, stringsAsFactors = FALSE),
      out,
      by = "time",
      all.x = TRUE,
      sort = TRUE
    )
    merged$cases[is.na(merged$cases)] <- 0L
    out <- merged
  }

  out
}

#' Summarise fit diagnostics and posterior summaries
#'
#' @param fit A `chlaa_fit` object.
#' @param burnin Burn-in proportion or count.
#' @param thin Thinning interval.
#' @param probs Quantiles for posterior summary.
#'
#' @return A list with acceptance rate, iterations retained, posterior summary, and trace data.
#' @export
chlaa_fit_report <- function(fit, burnin = 0.5, thin = 1, probs = c(0.025, 0.5, 0.975)) {
  fit <- chlaa_as_fit(fit)
  chain_draws <- .chlaa_fit_chain_draws(fit, burnin = burnin, thin = thin, scale = "sampled")
  param_cols <- setdiff(names(chain_draws), c("chain", "iteration"))

  acceptance_by_chain <- do.call(rbind, lapply(split(chain_draws, chain_draws$chain), function(d) {
    if (nrow(d) < 2) stop("Need at least 2 retained iterations per chain for diagnostics", call. = FALSE)
    step_changed <- apply(abs(diff(as.matrix(d[, param_cols, drop = FALSE]))) > 0, 1, any)
    data.frame(
      chain = d$chain[[1]],
      acceptance_rate = mean(step_changed),
      n_iterations = nrow(d),
      stringsAsFactors = FALSE
    )
  }))
  rownames(acceptance_by_chain) <- NULL
  if (requireNamespace("tibble", quietly = TRUE)) acceptance_by_chain <- tibble::as_tibble(acceptance_by_chain)

  acceptance_rate <- mean(acceptance_by_chain$acceptance_rate)

  trace <- chlaa_fit_trace(fit, burnin = burnin, thin = thin)
  summ <- chlaa_posterior_summary(fit, burnin = burnin, thin = thin, probs = probs)

  list(
    acceptance_rate = acceptance_rate,
    acceptance_by_chain = acceptance_by_chain,
    n_iterations = min(acceptance_by_chain$n_iterations),
    n_draws = nrow(chain_draws),
    n_parameters = length(param_cols),
    posterior_summary = summ,
    trace = trace
  )
}

#' Extract a long trace data frame from posterior draws
#'
#' @param fit A `chlaa_fit` object.
#' @param burnin Burn-in proportion or count.
#' @param thin Thinning interval.
#' @param parameters Optional vector of parameter names to keep.
#' @param scale Plot or return sampled MCMC coordinates (`"sampled"`) or
#'   unpacked model parameters (`"natural"`).
#'
#' @return A long data.frame with columns `chain`, `iteration`, `parameter`, `value`.
#' @export
chlaa_fit_trace <- function(fit,
                            burnin = 0.0,
                            thin = 1,
                            parameters = NULL,
                            scale = c("sampled", "natural")) {
  fit <- chlaa_as_fit(fit)
  scale <- match.arg(scale)
  draws <- .chlaa_fit_chain_draws(fit, burnin = burnin, thin = thin, scale = scale)
  param_cols <- setdiff(names(draws), c("chain", "iteration"))

  if (!is.null(parameters)) {
    keep <- intersect(parameters, param_cols)
    if (length(keep) == 0) stop("No requested parameters found in draws", call. = FALSE)
    param_cols <- keep
  }

  out <- do.call(rbind, lapply(param_cols, function(p) {
    data.frame(
      chain = draws$chain,
      iteration = draws$iteration,
      parameter = p,
      value = draws[[p]],
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL

  if (requireNamespace("tibble", quietly = TRUE)) out <- tibble::as_tibble(out)
  out
}

#' Plot fit traces for selected parameters
#'
#' @param fit A `chlaa_fit` object.
#' @param parameters Optional subset of parameter names.
#' @param burnin Burn-in proportion or count.
#' @param thin Thinning interval.
#' @param scale Plot sampled MCMC coordinates (`"sampled"`) or unpacked model
#'   parameters (`"natural"`).
#'
#' @return A ggplot object.
#' @export
chlaa_plot_trace <- function(fit,
                             parameters = NULL,
                             burnin = 0.0,
                             thin = 1,
                             scale = c("sampled", "natural")) {
  .require_suggested("ggplot2")
  scale <- match.arg(scale)
  tr <- chlaa_fit_trace(fit, burnin = burnin, thin = thin, parameters = parameters, scale = scale)

  ggplot2::ggplot(tr, ggplot2::aes(x = .data$iteration, y = .data$value, colour = .data$chain)) +
    ggplot2::geom_line(alpha = 0.7, linewidth = 0.3) +
    ggplot2::facet_wrap(~ .data$parameter, scales = "free_y") +
    ggplot2::labs(x = "Iteration", y = "Value", colour = "Chain", title = "pMCMC traces") +
    ggplot2::theme_minimal()
}

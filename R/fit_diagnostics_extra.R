# Extra diagnostics and plotting helpers for fit objects

.chlaa_fit_density_vector <- function(fit) {
  fit <- chlaa_as_fit(fit)
  d <- fit$density
  if (is.null(d)) stop("fit does not contain density values", call. = FALSE)

  if (is.vector(d)) return(as.numeric(d))
  if (is.matrix(d)) return(as.numeric(d))

  stop("Unsupported density structure in fit$density", call. = FALSE)
}

.chlaa_fit_density_matrix <- function(fit) {
  fit <- chlaa_as_fit(fit)
  d <- fit$density
  if (is.null(d)) stop("fit does not contain density values", call. = FALSE)

  if (is.vector(d)) {
    d <- matrix(as.numeric(d), ncol = 1)
  } else if (is.matrix(d)) {
    storage.mode(d) <- "double"
  } else {
    stop("Unsupported density structure in fit$density", call. = FALSE)
  }

  if (is.null(colnames(d)) || any(colnames(d) == "")) {
    colnames(d) <- .chlaa_fit_chain_names(ncol(d))
  }
  d
}

#' Extract Likelihood (Log-Density) Trace From A Fit
#'
#' @param fit Fit object returned by `chlaa_fit_pmcmc()`.
#' @param burnin Burn-in proportion in (0,1) or iteration count.
#' @param thin Thinning interval.
#'
#' @return A data.frame with columns `iteration` and `log_density`.
#' @export
chlaa_fit_density_trace <- function(fit, burnin = 0, thin = 1) {
  dens <- .chlaa_fit_density_matrix(fit)
  idx <- .chlaa_iteration_index(nrow(dens), burnin = burnin, thin = thin)

  rows <- lapply(seq_len(ncol(dens)), function(k) {
    data.frame(
      chain = colnames(dens)[[k]],
      iteration = idx,
      log_density = dens[idx, k],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (requireNamespace("tibble", quietly = TRUE)) out <- tibble::as_tibble(out)
  out
}

#' Plot The Likelihood Trace
#'
#' @param fit Fit object returned by `chlaa_fit_pmcmc()`.
#' @param burnin Burn-in proportion in (0,1) or iteration count.
#' @param thin Thinning interval.
#'
#' @return A ggplot object.
#' @export
chlaa_plot_likelihood_trace <- function(fit, burnin = 0, thin = 1) {
  .require_suggested("ggplot2")
  df <- chlaa_fit_density_trace(fit, burnin = burnin, thin = thin)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$log_density, colour = .data$chain)) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::labs(
      x = "Iteration",
      y = "Log posterior density",
      colour = "Chain",
      title = "Likelihood / posterior density trace"
    ) +
    ggplot2::theme_minimal()
}

#' Plot Distribution Of Likelihood Values
#'
#' @param fit Fit object returned by `chlaa_fit_pmcmc()`.
#' @param burnin Burn-in proportion in (0,1) or iteration count.
#' @param thin Thinning interval.
#'
#' @return A ggplot object.
#' @export
chlaa_plot_likelihood_density <- function(fit, burnin = 0, thin = 1) {
  .require_suggested("ggplot2")
  df <- chlaa_fit_density_trace(fit, burnin = burnin, thin = thin)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$log_density, fill = .data$chain)) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.35, position = "identity") +
    ggplot2::geom_density() +
    ggplot2::labs(
      x = "Log posterior density",
      y = "Count / density",
      fill = "Chain",
      title = "Distribution of sampled likelihood values"
    ) +
    ggplot2::theme_minimal()
}

#' Plot Pairwise Posterior Parameter Densities
#'
#' @param fit Fit object returned by `chlaa_fit_pmcmc()`.
#' @param parameters Optional parameter subset. Defaults to the first five.
#' @param burnin Burn-in proportion in (0,1) or iteration count.
#' @param thin Thinning interval.
#' @param max_points Maximum sampled posterior rows to plot.
#' @param scale Plot sampled MCMC coordinates (`"sampled"`) or unpacked model
#'   parameters (`"natural"`).
#' @param truth Optional named numeric vector of true/reference values.
#'
#' @return A ggplot object.
#' @export
chlaa_plot_parameter_pairs <- function(fit,
                                       parameters = NULL,
                                       burnin = 0.5,
                                       thin = 1,
                                       max_points = 2000,
                                       scale = c("natural", "sampled"),
                                       truth = NULL) {
  .require_suggested("ggplot2")
  scale <- match.arg(scale)

  draws <- .chlaa_fit_chain_draws(chlaa_as_fit(fit), burnin = burnin, thin = thin, scale = scale)
  param_cols <- setdiff(names(draws), c("chain", "iteration"))

  if (is.null(parameters)) {
    parameters <- param_cols[seq_len(min(5, length(param_cols)))]
  } else {
    parameters <- intersect(parameters, param_cols)
  }
  if (length(parameters) < 2) {
    stop("Need at least two parameters for pairwise plotting", call. = FALSE)
  }

  d <- as.data.frame(draws[, parameters, drop = FALSE], stringsAsFactors = FALSE)
  if (nrow(d) > max_points) {
    set.seed(1)
    d <- d[sample.int(nrow(d), max_points), , drop = FALSE]
  }
  truth <- .chlaa_truth_for_parameters(truth, parameters)

  if (requireNamespace("GGally", quietly = TRUE)) {
    p <- GGally::ggpairs(
      d,
      lower = list(continuous = GGally::wrap(.chlaa_ggpairs_lower, truth = truth)),
      diag = list(continuous = GGally::wrap(.chlaa_ggpairs_diag, truth = truth)),
      upper = list(continuous = GGally::wrap("blankDiag"))
    )
    return(
      p +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Posterior parameter corner plot")
    )
  }

  .chlaa_plot_parameter_pairs_fallback(d, parameters, truth)
}

#' Plot Marginal Posterior Parameter Distributions
#'
#' @param fit Fit object returned by `chlaa_fit_pmcmc()`.
#' @param parameters Optional parameter subset. Defaults to the first six.
#' @param burnin Burn-in proportion in (0,1) or iteration count.
#' @param thin Thinning interval.
#' @param scale Plot sampled MCMC coordinates (`"sampled"`) or unpacked model
#'   parameters (`"natural"`).
#' @param truth Optional named numeric vector of true/reference values.
#'
#' @return A ggplot object.
#' @export
chlaa_plot_parameter_distributions <- function(fit,
                                               parameters = NULL,
                                               burnin = 0.5,
                                               thin = 1,
                                               scale = c("natural", "sampled"),
                                               truth = NULL) {
  .require_suggested("ggplot2")
  scale <- match.arg(scale)

  draws <- .chlaa_fit_chain_draws(chlaa_as_fit(fit), burnin = burnin, thin = thin, scale = scale)
  param_cols <- setdiff(names(draws), c("chain", "iteration"))
  if (is.null(parameters)) {
    parameters <- param_cols[seq_len(min(6, length(param_cols)))]
  } else {
    parameters <- intersect(parameters, param_cols)
  }
  if (length(parameters) < 1) stop("No requested parameters found in draws", call. = FALSE)

  long <- do.call(rbind, lapply(parameters, function(p) {
    data.frame(
      chain = draws$chain,
      parameter = p,
      value = draws[[p]],
      stringsAsFactors = FALSE
    )
  }))
  rownames(long) <- NULL
  truth <- .chlaa_truth_for_parameters(truth, parameters)
  truth_df <- if (is.null(truth)) NULL else data.frame(parameter = names(truth), value = as.numeric(truth))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value, colour = .data$chain, fill = .data$chain)) +
    ggplot2::geom_density(alpha = 0.2) +
    ggplot2::facet_wrap(~ .data$parameter, scales = "free", ncol = 2) +
    ggplot2::labs(
      x = "Parameter value",
      y = "Posterior density",
      colour = "Chain",
      fill = "Chain",
      title = "Marginal posterior parameter distributions"
    ) +
    ggplot2::theme_minimal()

  if (!is.null(truth_df)) {
    p <- p + ggplot2::geom_vline(
      data = truth_df,
      ggplot2::aes(xintercept = .data$value),
      inherit.aes = FALSE,
      linetype = 2,
      linewidth = 0.5
    )
  }

  p
}

.chlaa_truth_for_parameters <- function(truth, parameters) {
  if (is.null(truth)) return(NULL)
  if (!is.numeric(truth)) stop("truth must be a named numeric vector", call. = FALSE)
  if (is.null(names(truth)) || any(names(truth) == "")) {
    if (length(truth) != length(parameters)) {
      stop("truth must be named, or have the same length as parameters", call. = FALSE)
    }
    names(truth) <- parameters
  }

  out <- truth[intersect(parameters, names(truth))]
  if (length(out) == 0) return(NULL)
  out
}

.chlaa_ggpairs_var <- function(mapping, which) {
  vars <- all.vars(mapping[[which]])
  if (length(vars) == 0) NA_character_ else vars[[1]]
}

.chlaa_ggpairs_lower <- function(data, mapping, truth = NULL, ...) {
  xvar <- .chlaa_ggpairs_var(mapping, "x")
  yvar <- .chlaa_ggpairs_var(mapping, "y")
  p <- ggplot2::ggplot(data = data, mapping = mapping) +
    ggplot2::geom_point(alpha = 0.15, size = 0.35)

  if (requireNamespace("MASS", quietly = TRUE)) {
    p <- p + ggplot2::geom_density_2d(linewidth = 0.25, alpha = 0.7)
  }

  if (!is.null(truth)) {
    if (!is.na(xvar) && xvar %in% names(truth)) {
      p <- p + ggplot2::geom_vline(xintercept = truth[[xvar]], linetype = 2, linewidth = 0.35)
    }
    if (!is.na(yvar) && yvar %in% names(truth)) {
      p <- p + ggplot2::geom_hline(yintercept = truth[[yvar]], linetype = 2, linewidth = 0.35)
    }
    if (!is.na(xvar) && !is.na(yvar) && all(c(xvar, yvar) %in% names(truth))) {
      p <- p + ggplot2::annotate("point", x = truth[[xvar]], y = truth[[yvar]], size = 1.2)
    }
  }

  p
}

.chlaa_ggpairs_diag <- function(data, mapping, truth = NULL, ...) {
  xvar <- .chlaa_ggpairs_var(mapping, "x")
  p <- ggplot2::ggplot(data = data, mapping = mapping) +
    ggplot2::geom_density(alpha = 0.6, fill = "#dbe9f6", colour = "#3a6ea5")

  if (!is.null(truth) && !is.na(xvar) && xvar %in% names(truth)) {
    p <- p + ggplot2::geom_vline(xintercept = truth[[xvar]], linetype = 2, linewidth = 0.35)
  }

  p
}

.chlaa_plot_parameter_pairs_fallback <- function(d, parameters, truth = NULL) {
  pairs <- utils::combn(parameters, 2, simplify = FALSE)
  long <- do.call(rbind, lapply(pairs, function(pp) {
    data.frame(
      param_x = pp[1],
      param_y = pp[2],
      x = d[[pp[1]]],
      y = d[[pp[2]]],
      stringsAsFactors = FALSE
    )
  }))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_bin2d(bins = 35) +
    ggplot2::facet_grid(.data$param_y ~ .data$param_x, scales = "free") +
    ggplot2::scale_fill_viridis_c(option = "C", guide = "none") +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = "Pairwise posterior parameter densities (fallback)"
    ) +
    ggplot2::theme_minimal()

  if (!is.null(truth)) {
    truth_pairs <- do.call(rbind, lapply(pairs, function(pp) {
      if (!all(pp %in% names(truth))) return(NULL)
      data.frame(param_x = pp[1], param_y = pp[2], x = truth[[pp[1]]], y = truth[[pp[2]]])
    }))
    if (!is.null(truth_pairs) && nrow(truth_pairs) > 0) {
      p <- p +
        ggplot2::geom_vline(
          data = truth_pairs,
          ggplot2::aes(xintercept = .data$x),
          inherit.aes = FALSE,
          linetype = 2
        ) +
        ggplot2::geom_hline(
          data = truth_pairs,
          ggplot2::aes(yintercept = .data$y),
          inherit.aes = FALSE,
          linetype = 2
        ) +
        ggplot2::geom_point(
          data = truth_pairs,
          ggplot2::aes(x = .data$x, y = .data$y),
          inherit.aes = FALSE,
          size = 1.2
        )
    }
  }

  p
}

#' Plot Parameters Against Sampled Likelihood
#'
#' @param fit Fit object returned by `chlaa_fit_pmcmc()`.
#' @param parameters Optional parameter subset. Defaults to the first six.
#' @param burnin Burn-in proportion in (0,1) or iteration count.
#' @param thin Thinning interval.
#' @param max_points Maximum sampled posterior rows to plot.
#'
#' @return A ggplot object.
#' @export
chlaa_plot_parameter_vs_likelihood <- function(fit,
                                               parameters = NULL,
                                               burnin = 0.5,
                                               thin = 1,
                                               max_points = 2000) {
  .require_suggested("ggplot2")

  fit <- chlaa_as_fit(fit)
  draws_df <- .chlaa_fit_chain_draws(fit, burnin = burnin, thin = thin, scale = "sampled")
  dens_df <- chlaa_fit_density_trace(fit, burnin = burnin, thin = thin)
  if (nrow(dens_df) != nrow(draws_df)) {
    stop("fit density does not match number of posterior draws", call. = FALSE)
  }
  draws <- as.matrix(draws_df[, setdiff(names(draws_df), c("chain", "iteration")), drop = FALSE])
  dens <- dens_df$log_density

  if (is.null(parameters)) {
    parameters <- colnames(draws)[seq_len(min(6, ncol(draws)))]
  } else {
    parameters <- intersect(parameters, colnames(draws))
  }
  if (length(parameters) < 1) stop("No requested parameters found in draws", call. = FALSE)

  d <- as.data.frame(draws[, parameters, drop = FALSE], stringsAsFactors = FALSE)
  d$log_density <- dens
  if (nrow(d) > max_points) {
    set.seed(1)
    d <- d[sample.int(nrow(d), max_points), , drop = FALSE]
  }

  long <- do.call(rbind, lapply(parameters, function(p) {
    data.frame(parameter = p, value = d[[p]], log_density = d$log_density, stringsAsFactors = FALSE)
  }))

  ggplot2::ggplot(long, ggplot2::aes(x = .data$value, y = .data$log_density)) +
    ggplot2::geom_point(alpha = 0.2, size = 0.4) +
    ggplot2::geom_smooth(se = FALSE, linewidth = 0.6) +
    ggplot2::facet_wrap(~ .data$parameter, scales = "free_x") +
    ggplot2::labs(
      x = "Parameter value",
      y = "Log posterior density",
      title = "Parameter-likelihood relationship"
    ) +
    ggplot2::theme_minimal()
}

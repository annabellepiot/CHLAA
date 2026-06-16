# S3 wrapper utilities for chlaa_fit objects

#' Coerce an object into class `chlaa_fit`
#'
#' @param x A fit object (typically from `chlaa_fit_pmcmc()`).
#'
#' @return The same object with class `chlaa_fit` prepended.
#' @export
chlaa_as_fit <- function(x) {
  if (inherits(x, "chlaa_fit")) return(x)
  class(x) <- unique(c("chlaa_fit", class(x)))
  x
}

#' Extract fit metadata (packer, prior, start parameters, data)
#'
#' @param fit A `chlaa_fit` object.
#'
#' @return A named list.
#' @export
chlaa_fit_metadata <- function(fit) {
  fit <- chlaa_as_fit(fit)
  list(
    packer = attr(fit, "packer", exact = TRUE),
    prior = attr(fit, "prior", exact = TRUE),
    start_pars = attr(fit, "start_pars", exact = TRUE),
    chain_pars = attr(fit, "chain_pars", exact = TRUE),
    data = attr(fit, "data", exact = TRUE)
  )
}

#' Posterior summary table for a fit
#'
#' @param fit A `chlaa_fit` object.
#' @param burnin Burn-in proportion (0-1) or integer iterations.
#' @param thin Thinning interval.
#' @param probs Quantiles to compute.
#'
#' @return A data.frame with columns: parameter, mean, sd, and quantiles.
#' @export
chlaa_posterior_summary <- function(fit, burnin = 0.5, thin = 1, probs = c(0.025, 0.5, 0.975)) {
  fit <- chlaa_as_fit(fit)
  draws <- .chlaa_fit_selected_draws_matrix(fit, burnin = burnin, thin = thin)

  q <- t(apply(draws, 2, stats::quantile, probs = probs, names = TRUE))
  mu <- colMeans(draws)
  sd <- apply(draws, 2, stats::sd)

  out <- data.frame(
    parameter = colnames(draws),
    mean = as.numeric(mu),
    sd = as.numeric(sd),
    stringsAsFactors = FALSE
  )

  qdf <- as.data.frame(q, stringsAsFactors = FALSE)
  colnames(qdf) <- paste0("q", gsub("\\.", "p", colnames(qdf)))
  out <- cbind(out, qdf)

  if (requireNamespace("tibble", quietly = TRUE)) out <- tibble::as_tibble(out)
  out
}

#' @export
print.chlaa_fit <- function(x, ...) {
  x <- chlaa_as_fit(x)

  cat("<chlaa_fit>\n")
  md <- chlaa_fit_metadata(x)

  n_iter <- NA_integer_
  n_chains <- NA_integer_
  n_par <- NA_integer_
  dr <- try(.chlaa_fit_draws_array(x), silent = TRUE)
  if (!inherits(dr, "try-error")) {
    n_iter <- dim(dr)[2]
    n_chains <- dim(dr)[3]
    n_par <- dim(dr)[1]
  }
  cat("Posterior draws: ", n_iter, " iterations x ", n_chains, " chains; ", n_par, " parameters\n", sep = "")

  if (is.data.frame(md$data) && all(c("time", "cases") %in% names(md$data))) {
    cat("Data: ", nrow(md$data), " observations; time range [",
        min(md$data$time), ", ", max(md$data$time), "]\n", sep = "")
  } else {
    cat("Data: <not attached>\n")
  }

  invisible(x)
}

#' @export
summary.chlaa_fit <- function(object, burnin = 0.5, thin = 1, probs = c(0.025, 0.5, 0.975), ...) {
  chlaa_posterior_summary(object, burnin = burnin, thin = thin, probs = probs)
}

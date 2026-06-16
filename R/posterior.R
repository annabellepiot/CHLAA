# Posterior utilities for chlaa

#' Extract posterior draws from a fit object
#'
#' Supports common structures including monty output where `fit$pars` is a 3D array:
#' (parameter, sample, chain).
#'
#' If column names are missing, it will try to use `attr(fit, "packer")$names`.
#'
#' @param fit Fit object returned by `chlaa_fit_pmcmc()` (or compatible).
#'
#' @return A numeric matrix with column names (iterations x parameters).
#' @export
chlaa_fit_draws <- function(fit) {
  a <- .chlaa_fit_draws_array(fit)
  pnames <- dimnames(a)[[1]]
  mats <- lapply(seq_len(dim(a)[3]), function(k) {
    t(a[, , k, drop = TRUE])
  })
  draws <- do.call(rbind, mats)
  colnames(draws) <- pnames
  draws
}

.chlaa_fit_draws_array <- function(fit) {
  draws <- NULL

  if (is.matrix(fit)) {
    draws <- fit
  } else if (is.data.frame(fit)) {
    draws <- as.matrix(fit)
  } else if (is.list(fit)) {
    if (!is.null(fit$pars) && is.array(fit$pars)) {
      a <- fit$pars
      d <- dim(a)

      if (length(d) == 3) {
        pnames <- dimnames(a)[[1]]
        draws <- a
        if (!is.null(pnames)) dimnames(draws)[[1]] <- pnames
        if (is.null(dimnames(draws)[[3]])) {
          dimnames(draws)[[3]] <- .chlaa_fit_chain_names(d[3])
        }
      } else if (length(d) == 2) {
        dn <- dimnames(a)
        rn <- dn[[1]]
        cn <- dn[[2]]

        if (!is.null(rn) && any(nzchar(rn))) {
          # monty-style: (parameter, sample)
          draws <- t(a)
          colnames(draws) <- rn
        } else if (!is.null(cn) && any(nzchar(cn))) {
          # already in (sample, parameter) orientation
          draws <- a
          colnames(draws) <- cn
        } else {
          draws <- a
        }
      }
    }

    if (is.null(draws)) {
      cand_names <- c("samples", "draws", "theta")
      for (nm in cand_names) {
        if (!is.null(fit[[nm]]) && (is.matrix(fit[[nm]]) || is.data.frame(fit[[nm]]))) {
          draws <- as.matrix(fit[[nm]])
          break
        }
      }
    }

    if (is.null(draws)) {
      idx <- which(vapply(fit, function(x) is.matrix(x) || is.data.frame(x), logical(1)))
      if (length(idx) > 0) {
        draws <- as.matrix(fit[[idx[1]]])
      }
    }
  }

  if (is.null(draws)) stop("Could not extract posterior draws from `fit`.", call. = FALSE)
  if (!is.numeric(draws)) storage.mode(draws) <- "double"

  if (is.array(draws) && length(dim(draws)) == 3) {
    if (is.null(dimnames(draws)[[1]]) || any(dimnames(draws)[[1]] == "")) {
      packer <- attr(fit, "packer", exact = TRUE)
      if (!is.null(packer) && !is.null(packer[["names"]])) {
        pnames <- packer[["names"]]()
        if (!is.null(pnames) && length(pnames) == dim(draws)[1]) {
          dimnames(draws)[[1]] <- pnames
        }
      }
    }
    if (is.null(dimnames(draws)[[1]]) || any(dimnames(draws)[[1]] == "")) {
      stop("Posterior draws are missing parameter names; cannot map draws to parameters.", call. = FALSE)
    }
    return(draws)
  }

  if (is.null(colnames(draws)) || any(colnames(draws) == "")) {
    packer <- attr(fit, "packer", exact = TRUE)
    if (!is.null(packer) && !is.null(packer[["names"]])) {
      pnames <- packer[["names"]]()
      if (!is.null(pnames) && length(pnames) == ncol(draws)) {
        colnames(draws) <- pnames
      }
    }
  }
  if (is.null(colnames(draws)) || any(colnames(draws) == "")) {
    stop("Posterior draws are missing column names; cannot map draws to parameters.", call. = FALSE)
  }

  a <- array(
    t(draws),
    dim = c(ncol(draws), nrow(draws), 1L),
    dimnames = list(colnames(draws), NULL, "chain_1")
  )
  a
}

.chlaa_fit_chain_names <- function(n_chains) {
  paste0("chain_", seq_len(n_chains))
}

.chlaa_iteration_index <- function(n, burnin = 0.5, thin = 1) {
  if (!is.numeric(thin) || length(thin) != 1 || thin < 1) stop("thin must be >= 1", call. = FALSE)
  if (n < 1) stop("draws has no rows", call. = FALSE)

  b <- burnin
  if (!is.numeric(b) || length(b) != 1 || b < 0) stop("burnin must be a non-negative number", call. = FALSE)

  if (b > 0 && b < 1) {
    start <- floor(b * n) + 1
  } else {
    start <- as.integer(b) + 1
  }
  start <- min(max(1, start), n)

  seq(from = start, to = n, by = as.integer(thin))
}

#' Select iterations from a posterior draws matrix
#'
#' @param draws Matrix of posterior draws (iterations x parameters).
#' @param burnin Burn-in, either proportion in (0,1) or an integer count.
#' @param thin Thinning interval (integer >= 1).
#'
#' @return A matrix subset of draws.
#' @export
chlaa_fit_select_iterations <- function(draws, burnin = 0.5, thin = 1) {
  if (!is.matrix(draws)) stop("draws must be a matrix", call. = FALSE)
  idx <- .chlaa_iteration_index(nrow(draws), burnin = burnin, thin = thin)
  draws[idx, , drop = FALSE]
}

.chlaa_fit_chain_draws <- function(fit,
                                   burnin = 0.5,
                                   thin = 1,
                                   scale = c("sampled", "natural")) {
  scale <- match.arg(scale)
  fit <- chlaa_as_fit(fit)
  a <- .chlaa_fit_draws_array(fit)

  idx <- .chlaa_iteration_index(dim(a)[2], burnin = burnin, thin = thin)
  sampled_names <- dimnames(a)[[1]]
  chain_names <- dimnames(a)[[3]]
  if (is.null(chain_names) || any(chain_names == "")) {
    chain_names <- .chlaa_fit_chain_names(dim(a)[3])
  }

  rows <- vector("list", length = dim(a)[3])
  for (k in seq_len(dim(a)[3])) {
    m <- t(a[, idx, k, drop = TRUE])
    colnames(m) <- sampled_names
    d <- as.data.frame(m, stringsAsFactors = FALSE)
    if (scale == "natural") {
      d <- .chlaa_natural_draws_from_sampled(d, fit)
    }
    d$chain <- chain_names[[k]]
    d$iteration <- idx
    rows[[k]] <- d[, c("chain", "iteration", setdiff(names(d), c("chain", "iteration"))), drop = FALSE]
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (requireNamespace("tibble", quietly = TRUE)) out <- tibble::as_tibble(out)
  out
}

.chlaa_natural_draws_from_sampled <- function(draws, fit) {
  packer <- attr(fit, "packer", exact = TRUE)
  if (is.null(packer) || is.null(packer[["unpack"]])) return(draws)

  sampled_names <- colnames(draws)
  first <- packer[["unpack"]](as.numeric(draws[1, sampled_names, drop = TRUE]))
  fixed_names <- character()
  if (!is.null(packer[["inputs"]])) {
    inputs <- packer[["inputs"]]()
    fixed_names <- names(inputs$fixed)
  }
  update_names <- setdiff(names(first), fixed_names)
  process_names <- setdiff(update_names, sampled_names)
  natural_names <- if (length(process_names) > 0) process_names else update_names

  rows <- lapply(seq_len(nrow(draws)), function(i) {
    theta <- as.numeric(draws[i, sampled_names, drop = TRUE])
    names(theta) <- sampled_names
    unpacked <- packer[["unpack"]](theta)
    unlist(unpacked[natural_names], use.names = TRUE)
  })

  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(out) <- natural_names
  out
}

.chlaa_fit_selected_draws_matrix <- function(fit,
                                             burnin = 0.5,
                                             thin = 1,
                                             scale = c("sampled", "natural")) {
  scale <- match.arg(scale)
  draws <- .chlaa_fit_chain_draws(fit, burnin = burnin, thin = thin, scale = scale)
  param_cols <- setdiff(names(draws), c("chain", "iteration"))
  out <- as.matrix(draws[, param_cols, drop = FALSE])
  colnames(out) <- param_cols
  out
}

.chlaa_update_pars_from_theta <- function(theta, pars, fit = NULL) {
  p <- pars
  packer <- if (!is.null(fit)) attr(fit, "packer", exact = TRUE) else NULL

  if (!is.null(packer) && !is.null(packer[["unpack"]])) {
    unpacked <- packer[["unpack"]](theta)
    inputs <- packer[["inputs"]]()
    fixed_names <- names(inputs$fixed)
    update_names <- setdiff(names(unpacked), fixed_names)

    for (nm in intersect(update_names, names(p))) {
      p[[nm]] <- unpacked[[nm]]
    }
  } else {
    common <- intersect(names(theta), names(p))
    if (length(common) > 0) p[common] <- as.list(as.numeric(theta[common]))
  }

  p
}

#' Update a parameter list using posterior information from a fit
#'
#' @param fit Fit object returned by `chlaa_fit_pmcmc()`.
#' @param pars Baseline parameter list.
#' @param draw Which posterior summary to use: "mean", "median", "sample", or "index".
#' @param burnin Burn-in, proportion or integer count.
#' @param thin Thinning interval.
#' @param index If `draw = "index"`, the 1-based index into the retained iterations.
#' @param seed Seed used when `draw = "sample"`.
#' @param validate Validate the resulting parameter list.
#'
#' @return A named list of parameters.
#' @export
chlaa_update_from_fit <- function(fit,
                                   pars,
                                   draw = c("median", "mean", "sample", "index"),
                                   burnin = 0.5,
                                   thin = 1,
                                   index = NULL,
                                   seed = 1,
                                   validate = TRUE) {
  draw <- match.arg(draw)
  .check_named_list(pars, "pars")

  draws2 <- .chlaa_fit_selected_draws_matrix(fit, burnin = burnin, thin = thin)
  if (nrow(draws2) < 1) stop("No posterior iterations remaining after burn-in/thinning.", call. = FALSE)

  theta <- switch(
    draw,
    mean = colMeans(draws2),
    median = apply(draws2, 2, stats::median),
    sample = {
      set.seed(seed)
      draws2[sample.int(nrow(draws2), 1), ]
    },
    index = {
      if (is.null(index) || !is.numeric(index) || length(index) != 1) {
        stop("When draw = 'index', provide a single numeric `index`.", call. = FALSE)
      }
      ii <- as.integer(index)
      if (ii < 1 || ii > nrow(draws2)) stop("index out of range for retained iterations.", call. = FALSE)
      draws2[ii, ]
    }
  )

  theta <- as.numeric(theta)
  names(theta) <- colnames(draws2)

  out <- .chlaa_update_pars_from_theta(theta, pars, fit)

  if (isTRUE(validate)) chlaa_parameters_validate(out)
  out
}

#' Create AA scenario set using a fitted posterior baseline
#'
#' Convenience wrapper: updates parameters from posterior draws and calls `chlaa_make_aa_scenarios()`.
#'
#' @param fit Fit object.
#' @param pars Baseline parameter list.
#' @param draw Posterior draw selection passed to `chlaa_update_from_fit()`.
#' @param burnin,thin,seed Passed to `chlaa_update_from_fit()`.
#' @param trigger_time Optional explicit trigger time.
#' @param trigger_threshold If trigger_time is NULL, derive from baseline simulation using this threshold.
#' @param trigger_time_var Variable used for thresholding (default `inc_symptoms`).
#' @param trigger_sim_time Time vector used for baseline simulation when deriving trigger.
#' @param trigger_sim_particles Particles used for baseline simulation when deriving trigger.
#' @param dt Model time step.
#' @param horizon Optional cap on intervention end times.
#' @param ... Other arguments forwarded to `chlaa_make_aa_scenarios()` (e.g. vax_total_doses).
#'
#' @return A list of scenarios with attribute `baseline_pars`.
#' @export
chlaa_make_aa_scenarios_from_fit <- function(fit,
                                               pars,
                                               draw = c("median", "mean", "sample", "index"),
                                               burnin = 0.5,
                                               thin = 1,
                                               seed = 1,
                                               trigger_time = NULL,
                                               trigger_threshold = 15,
                                               trigger_time_var = "inc_symptoms",
                                               trigger_sim_time = 0:365,
                                               trigger_sim_particles = 50,
                                               dt = 0.25,
                                               horizon = NULL,
                                               ...) {
  draw <- match.arg(draw)

  pars_fit <- chlaa_update_from_fit(
    fit = fit, pars = pars, draw = draw,
    burnin = burnin, thin = thin, seed = seed,
    validate = TRUE
  )

  if (is.null(horizon)) horizon <- max(trigger_sim_time) + 1

  if (is.null(trigger_time)) {
    base_sim <- chlaa_simulate(
      pars_fit, time = trigger_sim_time,
      n_particles = trigger_sim_particles, dt = dt, seed = seed
    )
    trigger_time <- chlaa_trigger_time_from_sim(
      base_sim, threshold = trigger_threshold, var = trigger_time_var
    )
    if (is.na(trigger_time)) stop("Trigger threshold was never reached in the baseline simulation.", call. = FALSE)
  }

  sc <- chlaa_make_aa_scenarios(
    pars = pars_fit,
    trigger_time = trigger_time,
    horizon = horizon,
    ...
  )
  attr(sc, "baseline_pars") <- pars_fit
  sc
}

#' Posterior predictive simulation
#'
#' Runs the simulation model for multiple posterior draws and returns a long data frame.
#'
#' @param fit Fit object.
#' @param pars Baseline parameter list.
#' @param time Simulation times.
#' @param n_draws Number of posterior draws to simulate.
#' @param burnin,thin Passed to draw selection.
#' @param seed Seed.
#' @param dt Model time step.
#' @param n_particles Particles per posterior draw.
#'
#' @return A data.frame with columns: draw, time, particle, plus model variables.
#' @export
chlaa_simulate_posterior <- function(fit,
                                      pars,
                                      time,
                                      n_draws = 50,
                                      burnin = 0.5,
                                      thin = 1,
                                      seed = 1,
                                      dt = 0.25,
                                      n_particles = 1) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required", call. = FALSE)

  draws2 <- .chlaa_fit_selected_draws_matrix(fit, burnin = burnin, thin = thin)
  if (nrow(draws2) < 1) stop("No posterior iterations remaining after burn-in/thinning.", call. = FALSE)

  set.seed(seed)
  idx <- sample.int(nrow(draws2), size = n_draws, replace = n_draws > nrow(draws2))

  out <- vector("list", length(idx))
  for (i in seq_along(idx)) {
    theta <- draws2[idx[i], , drop = TRUE]
    p <- .chlaa_update_pars_from_theta(theta, pars, fit)
    chlaa_parameters_validate(p)

    sim <- chlaa_simulate(p, time = time, n_particles = n_particles, dt = dt, seed = seed + i)
    sim$draw <- i
    out[[i]] <- sim
  }

  dplyr::bind_rows(out) |>
    dplyr::relocate(.data$draw, .data$time, .data$particle)
}

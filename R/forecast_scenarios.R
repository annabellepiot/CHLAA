# Scenario forecasts from a fitted posterior

.chlaa_normalise_scenarios_input <- function(scenarios) {
  if (is.null(scenarios)) return(list())

  if (is.data.frame(scenarios)) {
    if (!all(c("scenario", "modify") %in% names(scenarios))) {
      stop("If scenarios is a data.frame it must have columns: scenario, modify", call. = FALSE)
    }
    out <- vector("list", nrow(scenarios))
    for (i in seq_len(nrow(scenarios))) {
      out[[i]] <- chlaa_scenario(scenarios$scenario[[i]], scenarios$modify[[i]])
    }
    return(out)
  }

  if (is.list(scenarios) && length(scenarios) > 0) {
    if (all(vapply(scenarios, function(x) inherits(x, "chlaa_scenario"), logical(1)))) {
      return(scenarios)
    }

    if (!is.null(names(scenarios)) && all(names(scenarios) != "")) {
      out <- lapply(names(scenarios), function(nm) chlaa_scenario(nm, scenarios[[nm]]))
      return(out)
    }
  }

  stop(
    "scenarios must be one of: (1) list of chlaa_scenario, (2) named list of modify lists, ",
    "(3) data.frame with columns scenario and modify",
    call. = FALSE
  )
}

.chlaa_quantile_colnames <- function(q) {
  paste0("q", gsub("\\.", "p", as.character(q)))
}

.chlaa_simulate_posterior_matrix <- function(draws,
                                               idx,
                                               fit,
                                               base_pars,
                                               modify = NULL,
                                               time,
                                               vars_use,
                                               include_cases,
                                               obs_model,
                                               obs_interval,
                                               dt,
                                               seed,
                                               n_particles,
                                               n_threads,
                                               deterministic) {
  if (!is.null(modify) && length(modify) > 0) .check_named_list(modify, "modify")
  .check_named_list(base_pars, "base_pars")

  T <- length(time)
  n_draws <- length(idx)
  ns <- n_draws * n_particles

  mats <- lapply(vars_use, function(v) matrix(NA_real_, nrow = ns, ncol = T))
  names(mats) <- vars_use

  row0 <- 0L
  for (i in seq_len(n_draws)) {
    theta <- draws[idx[i], , drop = TRUE]
    p <- .chlaa_update_pars_from_theta(theta, base_pars, fit)
    if (!is.null(modify) && length(modify) > 0) {
      p <- utils::modifyList(p, modify)
    }
    chlaa_parameters_validate(p)

    sim <- chlaa_simulate(
      pars = p,
      time = time,
      n_particles = n_particles,
      dt = dt,
      seed = seed + i,
      n_threads = n_threads,
      deterministic = deterministic
    )

    for (v in setdiff(vars_use, "cases")) {
      if (!v %in% names(sim)) stop("Variable not in simulation output: ", v, call. = FALSE)
      m <- matrix(sim[[v]], nrow = n_particles, ncol = T)
      mats[[v]][(row0 + 1L):(row0 + n_particles), ] <- m
    }

    if (isTRUE(include_cases)) {
      obs_incidence_var <- .chlaa_obs_incidence_var(obs_interval)
      if (!(obs_incidence_var %in% names(sim))) {
        stop(obs_incidence_var, " required to generate cases", call. = FALSE)
      }
      if (!all(c("reporting_rate", "obs_size") %in% names(p))) {
        stop("reporting_rate and obs_size must be present in parameters", call. = FALSE)
      }

      mu <- pmax(0, p$reporting_rate * sim[[obs_incidence_var]])
      cases_vec <- if (obs_model == "mean") {
        mu
      } else {
        stats::rnbinom(n = length(mu), mu = mu, size = p$obs_size)
      }

      mats[["cases"]][(row0 + 1L):(row0 + n_particles), ] <- matrix(cases_vec, nrow = n_particles, ncol = T)
    }

    row0 <- row0 + n_particles
  }

  mats
}

.chlaa_summarise_mats <- function(mats, time, scenario, type, quantiles) {
  qcols <- .chlaa_quantile_colnames(quantiles)
  out_list <- vector("list", length(mats))

  for (j in seq_along(mats)) {
    v <- names(mats)[j]
    m <- mats[[j]]

    q <- apply(m, 2, stats::quantile, probs = quantiles, names = FALSE)
    if (!is.matrix(q)) q <- matrix(q, nrow = length(quantiles), ncol = length(time))

    df <- data.frame(
      scenario = scenario,
      type = type,
      time = time,
      variable = v,
      mean = colMeans(m),
      stringsAsFactors = FALSE
    )
    for (k in seq_along(quantiles)) {
      df[[qcols[k]]] <- as.numeric(q[k, ])
    }
    df$n_samples <- nrow(m)
    out_list[[j]] <- df
  }

  do.call(rbind, out_list)
}

#' Forecast multiple scenarios from a fitted posterior
#'
#' Selects posterior iterations once (shared across scenarios), simulates baseline and each scenario
#' using the same draw indices and RNG seed streams, and returns both absolute forecasts and paired
#' differences vs baseline.
#'
#' @param fit A `chlaa_fit` object.
#' @param pars Baseline parameter list. If NULL uses `attr(fit, "start_pars")` else `chlaa_parameters()`.
#' @param scenarios Scenarios to run (list of `chlaa_scenario`, named list of modify lists, or a grid data.frame).
#' @param baseline_name Baseline scenario name (modify list may be empty).
#' @param time Simulation times. If NULL uses `fit` data times.
#' @param vars Model variables to summarise.
#' @param include_cases Include predicted observed cases variable "cases".
#' @param obs_model One of "nbinom" or "mean".
#' @param obs_interval Observation interval in days for generated observed
#'   cases. If NULL, uses `attr(fit, "obs_interval")`, falling back to 1.
#' @param quantiles Quantiles to compute.
#' @param n_draws Number of posterior draws to use.
#' @param burnin Burn-in proportion or integer.
#' @param thin Thinning interval.
#' @param seed Seed.
#' @param dt Model time step.
#' @param n_particles Particles per draw.
#' @param n_threads Threads for dust2.
#' @param deterministic Deterministic process model toggle (if supported).
#' @param include_baseline_in_scenarios If TRUE, ensures baseline is included even if not passed.
#'
#' @return A tidy data.frame with columns: scenario, type, time, variable, mean, quantiles, n_samples.
#' @export
chlaa_forecast_scenarios_from_fit <- function(fit,
                                                pars = NULL,
                                                scenarios = NULL,
                                                baseline_name = "baseline",
                                                time = NULL,
                                                vars = c("inc_symptoms", "cum_symptoms", "cum_deaths"),
                                                include_cases = TRUE,
                                                obs_model = c("nbinom", "mean"),
                                                obs_interval = NULL,
                                                quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                                n_draws = 100,
                                                burnin = 0.5,
                                                thin = 1,
                                                seed = 1,
                                                dt = 0.25,
                                                n_particles = 1,
                                                n_threads = 1,
                                                deterministic = FALSE,
                                                include_baseline_in_scenarios = TRUE) {
  obs_model <- match.arg(obs_model)
  fit <- chlaa_as_fit(fit)

  if (is.null(time)) {
    dat <- attr(fit, "data", exact = TRUE)
    if (is.data.frame(dat) && "time" %in% names(dat)) {
      time <- sort(unique(dat$time))
    } else {
      stop("time must be provided (or fit must have data with a time column)", call. = FALSE)
    }
  }
  time <- sort(as.numeric(time))

  if (is.null(pars)) {
    pars <- attr(fit, "start_pars", exact = TRUE)
    if (is.null(pars)) pars <- chlaa_parameters()
  }
  .check_named_list(pars, "pars")
  chlaa_parameters_validate(pars)
  obs_interval <- .chlaa_forecast_obs_interval(obs_interval, fit)

  draws <- .chlaa_fit_selected_draws_matrix(fit, burnin = burnin, thin = thin)
  if (nrow(draws) < 1) stop("No posterior iterations remain after burn-in/thinning", call. = FALSE)

  set.seed(seed)
  idx <- sample.int(nrow(draws), size = n_draws, replace = n_draws > nrow(draws))

  vars_use <- vars
  if (isTRUE(include_cases)) vars_use <- unique(c(vars_use, "cases"))

  sc_list <- .chlaa_normalise_scenarios_input(scenarios)

  has_baseline <- any(vapply(sc_list, function(s) identical(s$name, baseline_name), logical(1)))
  if (isTRUE(include_baseline_in_scenarios) && !has_baseline) {
    sc_list <- c(list(chlaa_scenario(baseline_name, list())), sc_list)
  }

  base_modify <- list()
  for (s in sc_list) {
    if (identical(s$name, baseline_name)) {
      base_modify <- s$modify
      break
    }
  }
  if (length(base_modify) > 0) .check_named_list(base_modify, "baseline modify")

  mats_base <- .chlaa_simulate_posterior_matrix(
    draws = draws,
    idx = idx,
    fit = fit,
    base_pars = pars,
    modify = base_modify,
    time = time,
    vars_use = vars_use,
    include_cases = include_cases,
    obs_model = obs_model,
    obs_interval = obs_interval,
    dt = dt,
    seed = seed,
    n_particles = n_particles,
    n_threads = n_threads,
    deterministic = deterministic
  )

  out <- list()

  out[[length(out) + 1L]] <- .chlaa_summarise_mats(
    mats = mats_base, time = time, scenario = baseline_name, type = "absolute", quantiles = quantiles
  )

  mats0 <- lapply(mats_base, function(m) m - m)
  out[[length(out) + 1L]] <- .chlaa_summarise_mats(
    mats = mats0, time = time, scenario = baseline_name, type = "difference", quantiles = quantiles
  )

  for (s in sc_list) {
    if (identical(s$name, baseline_name)) next

    mats_s <- .chlaa_simulate_posterior_matrix(
      draws = draws,
      idx = idx,
      fit = fit,
      base_pars = pars,
      modify = s$modify,
      time = time,
      vars_use = vars_use,
      include_cases = include_cases,
      obs_model = obs_model,
      obs_interval = obs_interval,
      dt = dt,
      seed = seed,
      n_particles = n_particles,
      n_threads = n_threads,
      deterministic = deterministic
    )

    out[[length(out) + 1L]] <- .chlaa_summarise_mats(
      mats = mats_s, time = time, scenario = s$name, type = "absolute", quantiles = quantiles
    )

    mats_diff <- lapply(vars_use, function(v) mats_s[[v]] - mats_base[[v]])
    names(mats_diff) <- vars_use
    out[[length(out) + 1L]] <- .chlaa_summarise_mats(
      mats = mats_diff, time = time, scenario = s$name, type = "difference", quantiles = quantiles
    )
  }

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  if (requireNamespace("tibble", quietly = TRUE)) res <- tibble::as_tibble(res)

  attr(res, "baseline_name") <- baseline_name
  attr(res, "n_draws") <- n_draws
  attr(res, "n_particles") <- n_particles
  attr(res, "dt") <- dt
  attr(res, "draw_indices") <- idx
  attr(res, "quantiles") <- quantiles
  attr(res, "obs_interval") <- obs_interval
  res
}

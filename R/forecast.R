# Forecast and posterior predictive check utilities

#' Forecast summary from a fitted posterior
#'
#' Runs the model for multiple posterior draws and returns time-by-time summaries
#' (mean and quantiles) for selected variables.
#'
#' If `include_cases = TRUE`, this also generates a predictive distribution for observed
#' cases using the observation model:
#'   mu = reporting_rate * observed incidence
#' where observed incidence is daily `inc_symptoms` for `obs_interval = 1`, or
#' weekly `inc_symptoms_weekly` for `obs_interval = 7`.
#' and either:
#' - `obs_model = "nbinom"`: sample Negative Binomial noise
#' - `obs_model = "mean"`: use mu directly
#'
#' @param fit A `chlaa_fit` object (or compatible).
#' @param pars Baseline parameter list. If NULL, uses `attr(fit, "start_pars")`, otherwise `chlaa_parameters()`.
#' @param time Vector of times to simulate.
#' @param vars Character vector of model variables to summarise.
#' @param include_cases Logical; include predicted observed cases as variable "cases".
#' @param obs_model One of "nbinom" or "mean".
#' @param quantiles Numeric vector of quantiles to return.
#' @param n_draws Number of posterior draws to use.
#' @param burnin Burn-in proportion or integer.
#' @param thin Thinning interval.
#' @param seed Seed.
#' @param dt Model time step.
#' @param n_particles Particles per posterior draw.
#' @param n_threads Threads for dust2.
#' @param deterministic Run the process model deterministically (if supported).
#' @param obs_interval Observation interval in days for generated observed
#'   cases. If NULL, uses `attr(fit, "obs_interval")`, falling back to 1.
#' @param modify Optional named list of parameter modifications applied after each draw update.
#'
#' @return A tidy data.frame with columns: time, variable, mean, quantiles, n_samples.
#' @export
chlaa_forecast_from_fit <- function(fit,
                                      pars = NULL,
                                      time = NULL,
                                      vars = c("inc_symptoms", "cum_symptoms", "cum_deaths"),
                                      include_cases = TRUE,
                                      obs_model = c("nbinom", "mean"),
                                      quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                      n_draws = 100,
                                      burnin = 0.5,
                                      thin = 1,
                                      seed = 1,
                                      dt = 0.25,
                                      n_particles = 1,
                                      n_threads = 1,
                                      deterministic = FALSE,
                                      obs_interval = NULL,
                                      modify = NULL) {
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
  if (!is.numeric(time) || length(time) < 2) stop("time must be a numeric vector of length >= 2", call. = FALSE)

  if (is.null(pars)) {
    pars <- attr(fit, "start_pars", exact = TRUE)
    if (is.null(pars)) pars <- chlaa_parameters()
  }
  .check_named_list(pars, "pars")
  chlaa_parameters_validate(pars)

  if (!is.null(modify)) .check_named_list(modify, "modify")
  obs_interval <- .chlaa_forecast_obs_interval(obs_interval, fit)
  obs_incidence_var <- .chlaa_obs_incidence_var(obs_interval)

  draws <- .chlaa_fit_selected_draws_matrix(fit, burnin = burnin, thin = thin)
  if (nrow(draws) < 1) stop("No posterior iterations remain after burn-in/thinning", call. = FALSE)

  set.seed(seed)
  idx <- sample.int(nrow(draws), size = n_draws, replace = n_draws > nrow(draws))

  T <- length(time)
  ns <- n_draws * n_particles

  vars_use <- vars
  if (isTRUE(include_cases)) vars_use <- unique(c(vars_use, "cases"))

  mats <- lapply(vars_use, function(v) matrix(NA_real_, nrow = ns, ncol = T))
  names(mats) <- vars_use

  row0 <- 0L
  for (i in seq_len(n_draws)) {
    theta <- draws[idx[i], , drop = TRUE]
    p <- .chlaa_update_pars_from_theta(theta, pars, fit)

    if (!is.null(modify)) p <- utils::modifyList(p, modify)
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
      if (!(obs_incidence_var %in% names(sim))) {
        stop(obs_incidence_var, " required to generate observed cases", call. = FALSE)
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

  qnames <- paste0("q", gsub("\\.", "p", as.character(quantiles)))
  out_list <- vector("list", length(vars_use))

  for (j in seq_along(vars_use)) {
    v <- vars_use[j]
    m <- mats[[v]]
    q <- apply(m, 2, stats::quantile, probs = quantiles, names = FALSE)
    if (!is.matrix(q)) q <- matrix(q, nrow = length(quantiles), ncol = T)

    df <- data.frame(
      time = time,
      variable = v,
      mean = colMeans(m),
      stringsAsFactors = FALSE
    )
    for (k in seq_along(quantiles)) {
      df[[qnames[k]]] <- as.numeric(q[k, ])
    }
    df$n_samples <- ns
    out_list[[j]] <- df
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  if (requireNamespace("tibble", quietly = TRUE)) out <- tibble::as_tibble(out)

  attr(out, "n_draws") <- n_draws
  attr(out, "n_particles") <- n_particles
  attr(out, "dt") <- dt
  attr(out, "obs_interval") <- obs_interval
  out
}

.chlaa_forecast_obs_interval <- function(obs_interval, fit) {
  if (is.null(obs_interval)) {
    obs_interval <- attr(fit, "obs_interval", exact = TRUE)
    if (is.null(obs_interval)) obs_interval <- 1
  }
  .chlaa_obs_interval(obs_interval, observed_step = obs_interval)
}

.chlaa_obs_incidence_var <- function(obs_interval) {
  if (.chlaa_obs_interval(obs_interval, observed_step = obs_interval) == 7) {
    "inc_symptoms_weekly"
  } else {
    "inc_symptoms"
  }
}

#' Plot a forecast summary
#'
#' @param forecast Output from `chlaa_forecast_from_fit()`.
#' @param var Variable to plot.
#' @param data Optional observed data frame (must have a time column and a matching y column).
#' @param data_time Column name for time in `data`.
#' @param data_y Column name for y in `data`. If NULL, uses `var`.
#' @param show_mean Plot mean line as well as median.
#'
#' @return A ggplot object.
#' @export
chlaa_plot_forecast <- function(forecast,
                                  var,
                                  data = NULL,
                                  data_time = "time",
                                  data_y = NULL,
                                  show_mean = FALSE) {
  .require_suggested("ggplot2")

  if (!is.data.frame(forecast)) stop("forecast must be a data.frame", call. = FALSE)
  if (!all(c("time", "variable", "mean") %in% names(forecast))) stop("forecast missing required columns", call. = FALSE)
  if (!is.character(var) || length(var) != 1) stop("var must be a single string", call. = FALSE)

  f <- forecast[forecast$variable == var, , drop = FALSE]
  if (nrow(f) == 0) stop("No rows for var = ", var, call. = FALSE)

  req_q <- c("q0p025", "q0p25", "q0p5", "q0p75", "q0p975")
  missing <- setdiff(req_q, names(f))
  if (length(missing) > 0) {
    stop("forecast is missing quantile columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  p <- ggplot2::ggplot(f, ggplot2::aes(x = .data$time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$q0p025, ymax = .data$q0p975), alpha = 0.2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$q0p25, ymax = .data$q0p75), alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = .data$q0p5)) +
    ggplot2::labs(x = "Time (days)", y = var, title = paste("Forecast:", var)) +
    ggplot2::theme_minimal()

  if (isTRUE(show_mean)) {
    p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$mean), linetype = 2)
  }

  if (!is.null(data)) {
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

  p
}

#' Posterior predictive check plot for observed cases
#'
#' Convenience wrapper: builds a forecast over the observed data times including "cases" and overlays observations.
#'
#' @param fit A `chlaa_fit` object.
#' @param pars Baseline parameter list. If NULL uses `attr(fit, "start_pars")`.
#' @param data Observed data. If NULL uses `attr(fit, "data")`.
#' @param n_draws Number of posterior draws.
#' @param burnin Burn-in for posterior draws.
#' @param thin Thinning for posterior draws.
#' @param seed Seed.
#' @param dt Model time step.
#' @param n_particles Particles per draw.
#' @param obs_model Observation model noise: "nbinom" or "mean".
#'
#' @return A ggplot object.
#' @export
chlaa_plot_ppc <- function(fit,
                             pars = NULL,
                             data = NULL,
                             n_draws = 200,
                             burnin = 0.5,
                             thin = 1,
                             seed = 1,
                             dt = 0.25,
                             n_particles = 1,
                             obs_model = c("nbinom", "mean")) {
  obs_model <- match.arg(obs_model)
  fit <- chlaa_as_fit(fit)

  if (is.null(data)) data <- attr(fit, "data", exact = TRUE)
  if (!is.data.frame(data) || !all(c("time", "cases") %in% names(data))) {
    stop("data must be a data.frame with columns time and cases", call. = FALSE)
  }

  fc <- chlaa_forecast_from_fit(
    fit = fit,
    pars = pars,
    time = sort(unique(data$time)),
    vars = c("inc_symptoms"),
    include_cases = TRUE,
    obs_model = obs_model,
    n_draws = n_draws,
    burnin = burnin,
    thin = thin,
    seed = seed,
    dt = dt,
    n_particles = n_particles
  )

  chlaa_plot_forecast(
    forecast = fc,
    var = "cases",
    data = data,
    data_time = "time",
    data_y = "cases"
  )
}

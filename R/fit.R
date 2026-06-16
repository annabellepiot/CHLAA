# Fitting utilities (pMCMC)

#' Fit the cholera model to incidence data using pMCMC (monty + dust2)
#'
#' Expects a data.frame with columns:
#' - time: numeric times (days)
#' - cases: integer case counts at those times
#'
#' The likelihood is defined inside `inst/odin/cholera_model.R`:
#' cases ~ NegativeBinomial(mu = reporting_rate * observed incidence, size = obs_size)
#' where observed incidence is daily `inc_symptoms` when `obs_interval = 1`,
#' or weekly `inc_symptoms_weekly` when `obs_interval = 7`.
#'
#' @param data Data frame with columns time and cases.
#' @param pars Starting parameter list.
#' @param n_particles Number of particles for the dust2 filter likelihood.
#' @param n_steps Number of MCMC steps.
#' @param n_chains Number of MCMC chains.
#' @param chain_pars Optional list of per-chain starting parameter lists. If
#'   supplied, its length must equal `n_chains`; each element is packed with the
#'   same `packer` as `pars`.
#' @param seed Random seed.
#' @param prior Optional monty prior model. If NULL, uses `chlaa_default_prior()`.
#' @param packer Optional monty packer. If NULL, uses `chlaa_default_packer()`.
#' @param proposal_var Proposal variance for the random walk sampler. Can be:
#'   - A scalar: used as the diagonal value for all parameters (default 0.02).
#'   - A numeric vector (length = number of fitted parameters): per-parameter
#'     diagonal variances.
#'   - A square matrix: full variance-covariance proposal matrix.
#' @param obs_interval Observation interval in days: 1 for daily cases or 7 for
#'   weekly cases. If NULL, inferred from the smallest spacing in `data$time`.
#' @param time_start Optional filter start time. If NULL, inferred as the first
#'   observation time minus the smallest observed time step. For example,
#'   weekly data at `time = 7, 14, ...` starts from day 0, while daily data at
#'   `time = 0, 1, ...` starts from day -1 because dust requires the filter
#'   start to be strictly before the first observation.
#' @param deterministic Logical; if TRUE, use dust2's deterministic unfilter
#'   likelihood instead of the particle filter.
#'
#' @return A `chlaa_fit` object (also keeps monty class) with attributes:
#'   packer, prior, start_pars, chain_pars, data, obs_interval, time_start,
#'   deterministic.
#' @export
chlaa_fit_pmcmc <- function(data,
                            pars,
                            n_particles = 200,
                            n_steps = 2000,
                            n_chains = 1,
                            chain_pars = NULL,
                            seed = 1,
                            prior = NULL,
                            packer = NULL,
                            proposal_var = 0.02,
                            obs_interval = NULL,
                            time_start = NULL,
                            deterministic = FALSE) {
  .check_named_list(pars, "pars")
  if (!requireNamespace("monty", quietly = TRUE)) stop("monty is required for fitting", call. = FALSE)
  n_chains <- .chlaa_n_chains(n_chains)

  data <- chlaa_prepare_data(data, time_col = "time", cases_col = "cases")
  observed_step <- .chlaa_observed_step(data)
  obs_interval <- .chlaa_obs_interval(obs_interval, observed_step)
  if (nrow(data) == 1) observed_step <- obs_interval
  data$obs_interval <- obs_interval

  gen <- chlaa_generator()
  if (is.null(time_start)) {
    # dust2 expects the filter start time to be strictly earlier than the first
    # observation time. Infer the preceding observation boundary from the data
    # spacing so weekly observations at 7, 14, ... accumulate from day 0.
    time_start <- min(data$time) - observed_step
  } else if (!is.numeric(time_start) || length(time_start) != 1 || !is.finite(time_start)) {
    stop("time_start must be NULL or a single finite number", call. = FALSE)
  }
  if (time_start >= min(data$time)) {
    stop("time_start must be strictly earlier than the first observation time", call. = FALSE)
  }

  if (isTRUE(deterministic)) {
    filter <- dust2::dust_unfilter_create(gen, time_start = time_start, data = data)
  } else {
    filter <- dust2::dust_filter_create(gen, time_start = time_start, data = data, n_particles = n_particles)
  }

  if (is.null(packer)) {
    packer <- chlaa_default_packer(pars)
  }
  if (is.null(prior)) {
    prior <- chlaa_default_prior()
  }

  likelihood <- dust2::dust_likelihood_monty(filter, packer)
  posterior <- prior + likelihood

  packer_names <- packer[["names"]]()
  d <- length(packer_names)

  # Build the proposal variance-covariance matrix
  if (is.matrix(proposal_var)) {
    if (!all(dim(proposal_var) == d)) {
      stop(sprintf(
        "proposal_var matrix must be %d x %d (matching %d fitted parameters)",
        d, d, d
      ), call. = FALSE)
    }
    vcv <- proposal_var
  } else if (length(proposal_var) == d) {
    vcv <- diag(proposal_var)
  } else if (length(proposal_var) == 1) {
    vcv <- diag(d) * proposal_var
  } else {
    stop(sprintf(
      "proposal_var must be a scalar, a vector of length %d, or a %d x %d matrix",
      d, d, d
    ), call. = FALSE)
  }

  sampler <- monty::monty_sampler_random_walk(vcv)

  set.seed(seed)
  initial_vec <- .chlaa_initial_from_chain_pars(pars, chain_pars, packer, n_chains)

  res <- monty::monty_sample(posterior, sampler, n_steps, initial = initial_vec, n_chains = n_chains)

  attr(res, "packer") <- packer
  attr(res, "prior") <- prior
  attr(res, "start_pars") <- pars
  attr(res, "chain_pars") <- chain_pars
  attr(res, "data") <- data
  attr(res, "obs_interval") <- obs_interval
  attr(res, "time_start") <- time_start
  attr(res, "deterministic") <- isTRUE(deterministic)

  class(res) <- unique(c("chlaa_fit", class(res)))
  res
}

.chlaa_n_chains <- function(n_chains) {
  if (!is.numeric(n_chains) || length(n_chains) != 1 || !is.finite(n_chains) || n_chains < 1) {
    stop("n_chains must be a single positive integer", call. = FALSE)
  }
  if (abs(n_chains - round(n_chains)) > sqrt(.Machine$double.eps)) {
    stop("n_chains must be a single positive integer", call. = FALSE)
  }
  as.integer(n_chains)
}

.chlaa_initial_from_chain_pars <- function(pars, chain_pars, packer, n_chains) {
  if (is.null(chain_pars)) {
    return(packer[["pack"]](pars))
  }
  if (!is.list(chain_pars) || length(chain_pars) != n_chains) {
    stop("chain_pars must be NULL or a list of length n_chains", call. = FALSE)
  }

  pnames <- packer[["names"]]()
  initial <- vapply(seq_along(chain_pars), function(i) {
    .check_named_list(chain_pars[[i]], sprintf("chain_pars[[%d]]", i))
    x <- packer[["pack"]](chain_pars[[i]])
    if (length(x) != length(pnames)) {
      stop("Each element of chain_pars must pack to the same length as pars", call. = FALSE)
    }
    as.numeric(x)
  }, numeric(length(pnames)))
  rownames(initial) <- pnames
  colnames(initial) <- .chlaa_fit_chain_names(n_chains)
  initial
}

.chlaa_observed_step <- function(data) {
  if (nrow(data) > 1) min(diff(data$time)) else 1
}

.chlaa_obs_interval <- function(obs_interval, observed_step) {
  tol <- sqrt(.Machine$double.eps)
  if (is.null(obs_interval)) {
    if (abs(observed_step - 1) <= tol) {
      obs_interval <- 1
    } else if (abs(observed_step - 7) <= tol) {
      obs_interval <- 7
    } else {
      stop(
        "Could not infer obs_interval from the data spacing. ",
        "Please set obs_interval = 1 for daily data or obs_interval = 7 for weekly data.",
        call. = FALSE
      )
    }
  }

  if (!is.numeric(obs_interval) || length(obs_interval) != 1 || !is.finite(obs_interval)) {
    stop("obs_interval must be NULL, 1, or 7", call. = FALSE)
  }
  if (!any(abs(obs_interval - c(1, 7)) <= tol)) {
    stop("obs_interval must be 1 for daily cases or 7 for weekly cases", call. = FALSE)
  }

  as.numeric(if (abs(obs_interval - 7) <= tol) 7 else 1)
}

chlaa_default_prior <- function() {
  if (!requireNamespace("monty", quietly = TRUE)) stop("monty is required for fitting", call. = FALSE)
  monty::monty_dsl({
    trans_prob ~ Uniform(0.001, 0.2)
    contact_rate ~ Uniform(0.1, 30)
    incubation_time ~ Uniform(1, 10)
    duration_sym ~ Uniform(3, 30)
    reporting_rate ~ Uniform(0.01, 1.0)
    obs_size ~ Uniform(1, 200)
    seek_severe ~ Uniform(0.1, 1.0)
    fatality_untreated ~ Uniform(0.05, 0.9)
    fatality_treated ~ Uniform(0.0001, 0.05)
  })
}

chlaa_default_packer <- function(pars) {
  if (!requireNamespace("monty", quietly = TRUE)) stop("monty is required for fitting", call. = FALSE)
  .check_named_list(pars, "pars")

  # Keep this order consistent with `chlaa_default_prior()` (monty uses
  # positional parameter matching).
  names_fit <- c(
    "trans_prob", "contact_rate",
    "incubation_time", "duration_sym",
    "reporting_rate", "obs_size",
    "seek_severe",
    "fatality_untreated", "fatality_treated"
  )

  fixed <- pars[setdiff(names(pars), names_fit)]
  monty::monty_packer(names_fit, fixed = fixed)
}

# Budget optimisation utilities

#' Optimise intervention allocation under a budget constraint
#'
#' Simple grid-search optimiser for allocating budget across vaccination and WASH, with
#' optional constraints. This is a pragmatic starting point; swap to a more sophisticated
#' optimiser once costs and decision variables are finalised.
#'
#' The objective is to minimise expected deaths (mean across particles) over the horizon.
#'
#' @param pars Baseline parameter list.
#' @param budget Total budget.
#' @param cost Cost list (see Details).
#' @param time Simulation times.
#' @param n_particles Number of particles for each evaluation.
#' @param dt Time step.
#' @param seed Seed.
#' @param grid_size Number of grid points per decision dimension.
#' @param min_fraction Named list of minimum allocation fractions by intervention (`vax`, `wash`, `care`).
#' @param max_fraction Named list of maximum allocation fractions by intervention (`vax`, `wash`, `care`).
#' @param max_vax_doses_per_day Maximum feasible vaccination doses per day.
#' @param max_total_doses Optional upper bound on total vaccine doses.
#' @param method Optimisation method: "auto", "grid", or "continuous".
#'
#' @details
#' `cost` is a named list that can include:
#' - cost_per_vaccine_dose
#' - cost_chlorination_per_person_day
#' - cost_hygiene_per_person_day
#' - cost_latrine_per_person_day
#' - cost_cati_per_person_day
#'
#' Budget is spent on:
#' - vaccination doses (vax1_total_doses, vax1_doses_per_day within the campaign window)
#' - WASH (implemented by setting intervention effects; this is a placeholder mapping)
#'
#' @return A list with best allocation and a data.frame of evaluated allocations.
#' @export
chlaa_optimise_budget <- function(pars,
                                   budget,
                                   cost = list(
                                     cost_per_vaccine_dose = 2.0,
                                     cost_chlorination_per_person_day = 0.02,
                                     cost_hygiene_per_person_day = 0.03,
                                     cost_latrine_per_person_day = 0.01,
                                     cost_cati_per_person_day = 0.05,
                                     cost_per_orc_treatment = 10,
                                     cost_per_ctc_treatment = 80
                                   ),
                                   time = 0:180,
                                   n_particles = 100,
                                   dt = 0.25,
                                   seed = 1,
                                   grid_size = 10,
                                   min_fraction = list(vax = 0, wash = 0, care = 0),
                                   max_fraction = list(vax = 1, wash = 1, care = 1),
                                   max_vax_doses_per_day = Inf,
                                   max_total_doses = Inf,
                                   method = c("auto", "grid", "continuous")) {
  .check_named_list(pars, "pars")
  .check_named_list(cost, "cost")
  if (!is.numeric(budget) || length(budget) != 1 || budget <= 0) stop("budget must be > 0", call. = FALSE)
  method <- match.arg(method)
  if (method == "auto") method <- if (grid_size <= 20) "grid" else "continuous"

  mins <- c(vax = 0, wash = 0, care = 0)
  mins[names(min_fraction)] <- as.numeric(unlist(min_fraction))
  maxs <- c(vax = 1, wash = 1, care = 1)
  maxs[names(max_fraction)] <- as.numeric(unlist(max_fraction))

  if (any(mins < 0) || any(maxs > 1) || any(mins > maxs)) {
    stop("Invalid min/max allocation fraction bounds", call. = FALSE)
  }
  if (sum(mins) > 1 + 1e-12) stop("Sum of minimum fractions cannot exceed 1", call. = FALSE)

  .eval_allocation <- function(frac_vax, frac_wash, i_seed) {
    frac_care <- max(0, 1 - frac_vax - frac_wash)
    if (frac_vax < mins["vax"] || frac_vax > maxs["vax"]) return(NULL)
    if (frac_wash < mins["wash"] || frac_wash > maxs["wash"]) return(NULL)
    if (frac_care < mins["care"] || frac_care > maxs["care"]) return(NULL)

    b_vax <- budget * frac_vax
    b_wash <- budget * frac_wash
    b_care <- budget * frac_care

    doses <- floor(b_vax / cost$cost_per_vaccine_dose)
    doses <- min(doses, max_total_doses)
    dur <- max(time) - min(time) + 1
    feasible_doses <- floor(max_vax_doses_per_day * min(14, dur))
    doses <- min(doses, feasible_doses)

    denom <- pars$N * dur * (
      cost$cost_chlorination_per_person_day +
        cost$cost_hygiene_per_person_day +
        cost$cost_latrine_per_person_day +
        cost$cost_cati_per_person_day
    )
    intensity_wash <- if (denom > 0) min(1, b_wash / denom) else 0

    p <- pars
    dur_vax <- min(14, dur)
    p$vax1_start <- min(time)
    p$vax1_end <- min(time) + dur_vax
    p$vax1_total_doses <- doses
    p$vax1_doses_per_day <- if (dur_vax > 0) doses / dur_vax else 0

    p$chlor_start <- min(time); p$chlor_end <- max(time) + 1; p$chlor_effect <- 0.3 * intensity_wash
    p$hyg_start <- min(time); p$hyg_end <- max(time) + 1; p$hyg_effect <- 0.3 * intensity_wash
    p$lat_start <- min(time); p$lat_end <- max(time) + 1; p$lat_effect <- 0.2 * intensity_wash
    p$cati_start <- min(time); p$cati_end <- max(time) + 1; p$cati_effect <- 0.2 * intensity_wash

    care_scale <- if (b_care <= 0) 0 else sqrt(b_care / budget)
    p$orc_start <- min(time); p$orc_end <- max(time) + 1; p$orc_capacity <- pars$orc_capacity * care_scale
    p$ctc_start <- min(time); p$ctc_end <- max(time) + 1; p$ctc_capacity <- pars$ctc_capacity * care_scale

    chlaa_parameters_validate(p)

    sim <- chlaa_simulate(p, time = time, n_particles = n_particles, dt = dt, seed = i_seed)
    end_time <- max(sim$time)
    end <- sim[sim$time == end_time, , drop = FALSE]
    deaths <- mean(end$cum_deaths)
    cases <- mean(end$cum_symptoms)

    data.frame(
      frac_vax = frac_vax,
      frac_wash = frac_wash,
      frac_care = frac_care,
      budget_vax = b_vax,
      budget_wash = b_wash,
      budget_care = b_care,
      doses = doses,
      wash_intensity = intensity_wash,
      deaths = deaths,
      cases = cases,
      stringsAsFactors = FALSE
    )
  }

  eval <- list()
  if (method == "grid") {
    g <- seq(0, 1, length.out = grid_size)
    k <- 0L
    for (i in seq_along(g)) {
      for (j in seq_along(g)) {
        if (g[i] + g[j] > 1 + 1e-12) next
        tmp <- .eval_allocation(g[i], g[j], i_seed = seed + i + j)
        if (is.null(tmp)) next
        k <- k + 1L
        eval[[k]] <- tmp
      }
    }
  } else {
    # continuous optimisation over vax/wash fractions; care is residual
    starts <- list(c(0.5, 0.3), c(0.2, 0.2), c(0.7, 0.1))
    best_val <- Inf
    best_par <- NULL
    for (st in starts) {
      obj <- function(x) {
        fv <- min(max(x[1], 0), 1)
        fw <- min(max(x[2], 0), 1)
        if (fv + fw > 1) return(1e12 + 1e9 * (fv + fw - 1))
        row <- .eval_allocation(fv, fw, i_seed = seed + 1000)
        if (is.null(row)) return(1e12)
        row$deaths
      }
      opt <- stats::optim(par = st, fn = obj, method = "Nelder-Mead")
      fv <- min(max(opt$par[1], 0), 1)
      fw <- min(max(opt$par[2], 0), 1)
      if (fv + fw <= 1) {
        row <- .eval_allocation(fv, fw, i_seed = seed + 2000)
        if (!is.null(row)) {
          eval[[length(eval) + 1L]] <- row
          if (row$deaths < best_val) {
            best_val <- row$deaths
            best_par <- c(fv, fw)
          }
        }
      }
    }

    # add coarse grid support points for robustness
    g <- seq(0, 1, length.out = min(grid_size, 8))
    for (i in seq_along(g)) {
      for (j in seq_along(g)) {
        if (g[i] + g[j] > 1 + 1e-12) next
        row <- .eval_allocation(g[i], g[j], i_seed = seed + i + j + 3000)
        if (!is.null(row)) eval[[length(eval) + 1L]] <- row
      }
    }
  }

  if (length(eval) == 0) stop("No feasible allocations under given constraints", call. = FALSE)
  res <- unique(do.call(rbind, eval))
  best <- res[which.min(res$deaths), , drop = FALSE]
  list(best = best, evaluations = res, method = method)
}

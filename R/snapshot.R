# Snapshot helpers for branching counterfactual simulations

.chlaa_sim_to_df <- function(st, time, n_particles) {
  n_t <- length(time)
  out <- data.frame(
    time = rep(time, each = n_particles),
    particle = rep(seq_len(n_particles), times = n_t)
  )

  for (nm in names(st)) {
    x <- st[[nm]]
    if (is.matrix(x) && nrow(x) == n_particles && ncol(x) == n_t) {
      out[[nm]] <- as.vector(x)
    } else if (is.numeric(x) && length(x) == n_t) {
      out[[nm]] <- rep(x, each = n_particles)
    } else {
      out[[nm]] <- as.vector(x)
    }
  }

  out
}

#' Create A Simulation Snapshot At A Given Time
#'
#' Runs the model from initial state to `snapshot_time` and stores the full
#' model state plus RNG state. The snapshot can then be reused to branch
#' multiple counterfactual futures with identical history.
#'
#' @param pars Parameter list.
#' @param snapshot_time Time to snapshot.
#' @param n_particles Number of particles.
#' @param dt Model time step.
#' @param seed RNG seed.
#' @param n_threads Threads for dust2.
#' @param deterministic Deterministic mode (if supported).
#'
#' @return A list of class `chlaa_snapshot`.
#' @export
chlaa_snapshot_create <- function(pars,
                                  snapshot_time,
                                  n_particles = 1,
                                  dt = 0.25,
                                  seed = 1,
                                  n_threads = 1,
                                  deterministic = FALSE) {
  .check_named_list(pars, "pars")
  if (!is.numeric(snapshot_time) || length(snapshot_time) != 1) {
    stop("snapshot_time must be a single numeric value", call. = FALSE)
  }

  gen <- chlaa_generator()
  pars_use <- pars[names(pars) %in% attr(gen, "parameters")$name]
  sys <- dust2::dust_system_create(
    generator = gen,
    pars = pars_use,
    n_particles = n_particles,
    dt = dt,
    seed = seed,
    n_threads = n_threads,
    deterministic = deterministic,
    preserve_particle_dimension = TRUE
  )

  dust2::dust_system_set_state_initial(sys)
  if (snapshot_time > 0) {
    dust2::dust_system_run_to_time(sys, snapshot_time)
  }

  snap <- list(
    state = dust2::dust_system_state(sys),
    rng_state = dust2::dust_system_rng_state(sys),
    time = dust2::dust_system_time(sys),
    n_particles = n_particles,
    dt = dt,
    pars = pars
  )
  class(snap) <- c("chlaa_snapshot", class(snap))
  snap
}

#' Simulate Forward From A Snapshot
#'
#' @param snapshot Snapshot object from `chlaa_snapshot_create()`.
#' @param pars Optional parameter list to use from snapshot onward.
#' @param time Times to simulate to (must be >= `snapshot$time`).
#' @param n_threads Threads for dust2.
#' @param deterministic Deterministic mode (if supported).
#'
#' @return Simulation output data.frame as in `chlaa_simulate()`.
#' @export
chlaa_simulate_from_snapshot <- function(snapshot,
                                         pars = NULL,
                                         time,
                                         n_threads = 1,
                                         deterministic = FALSE) {
  if (!inherits(snapshot, "chlaa_snapshot")) {
    stop("snapshot must be a chlaa_snapshot object", call. = FALSE)
  }
  if (!is.numeric(time) || length(time) < 1) stop("time must be a non-empty numeric vector", call. = FALSE)
  time <- sort(as.numeric(time))
  if (min(time) < snapshot$time) {
    stop("all requested time points must be >= snapshot$time", call. = FALSE)
  }

  if (is.null(pars)) pars <- snapshot$pars
  .check_named_list(pars, "pars")

  gen <- chlaa_generator()
  pars_use <- pars[names(pars) %in% attr(gen, "parameters")$name]
  sys <- dust2::dust_system_create(
    generator = gen,
    pars = pars_use,
    n_particles = snapshot$n_particles,
    dt = snapshot$dt,
    seed = 1,
    n_threads = n_threads,
    deterministic = deterministic,
    preserve_particle_dimension = TRUE
  )

  dust2::dust_system_set_state(sys, snapshot$state)
  dust2::dust_system_set_rng_state(sys, snapshot$rng_state)
  dust2::dust_system_set_time(sys, snapshot$time)
  dust2::dust_system_update_pars(sys, pars_use)

  y <- dust2::dust_system_simulate(sys, time)
  st <- dust2::dust_unpack_state(sys, y)
  .chlaa_sim_to_df(st, time = time, n_particles = snapshot$n_particles)
}

#' Run Multiple Scenarios From A Shared Snapshot
#'
#' @param snapshot Snapshot from `chlaa_snapshot_create()`.
#' @param pars Baseline parameter list.
#' @param scenarios List of `chlaa_scenario` objects.
#' @param time Times to simulate to (must be >= snapshot time).
#' @param n_threads Threads for dust2.
#' @param deterministic Deterministic mode (if supported).
#'
#' @return Combined scenario simulation data.frame with `scenario` column.
#' @export
chlaa_run_scenarios_from_snapshot <- function(snapshot,
                                              pars,
                                              scenarios,
                                              time,
                                              n_threads = 1,
                                              deterministic = FALSE) {
  .check_named_list(pars, "pars")
  if (!inherits(snapshot, "chlaa_snapshot")) {
    stop("snapshot must be a chlaa_snapshot object", call. = FALSE)
  }
  if (!is.list(scenarios) || length(scenarios) == 0) {
    stop("scenarios must be a non-empty list of chlaa_scenario objects", call. = FALSE)
  }

  out <- lapply(scenarios, function(sc) {
    if (!inherits(sc, "chlaa_scenario")) stop("all scenarios must be chlaa_scenario objects", call. = FALSE)
    p <- utils::modifyList(pars, sc$modify)
    chlaa_parameters_validate(p)
    sim <- chlaa_simulate_from_snapshot(
      snapshot = snapshot,
      pars = p,
      time = time,
      n_threads = n_threads,
      deterministic = deterministic
    )
    sim$scenario <- sc$name
    sim
  })

  if (requireNamespace("dplyr", quietly = TRUE)) {
    dplyr::bind_rows(out) |>
      dplyr::relocate("scenario", "time", "particle")
  } else {
    res <- do.call(rbind, out)
    res[, c("scenario", "time", "particle", setdiff(names(res), c("scenario", "time", "particle")))]
  }
}

# Model interface: generator lookup + simulation wrapper

.as_dust_generator <- function(gen) {
  if (is.function(gen) && !inherits(gen, "dust_system_generator")) {
    gen <- gen()
  }
  if (!inherits(gen, "dust_system_generator")) {
    stop(
      "Expected a dust_system_generator. Run odin2::odin_package('.') to generate a packaged generator, ",
      "or install odin2 to compile on the fly.",
      call. = FALSE
    )
  }
  gen
}

.find_odin_file <- function(filename) {
  path <- system.file("odin", filename, package = "chlaa")
  if (is.character(path) && length(path) == 1 && nzchar(path)) {
    return(path)
  }
  cand <- file.path("inst", "odin", filename)
  if (file.exists(cand)) {
    return(cand)
  }
  ""
}

.runtime_odin2_allowed <- function() {
  opt <- getOption("chlaa.allow_runtime_odin2", FALSE)
  env <- Sys.getenv("CHLAA_ALLOW_RUNTIME_ODIN2", unset = "")
  isTRUE(opt) || tolower(env) %in% c("1", "true", "t", "yes", "y")
}

# Internal: get a dust2 generator for either simulation or fitting
chlaa_generator <- function() {

  # Keep these object names aligned with generated odin/dust symbols.
  obj <- "cholera_model"
  ns <- asNamespace("chlaa")

  # Packaged generator (present after odin2::odin_package())
  if (exists(obj, envir = ns, inherits = FALSE)) {
    gen <- get(obj, envir = ns, inherits = FALSE)
    return(.as_dust_generator(gen))
  }

  # Dev/runtime compilation fallback (explicitly opt-in).
  if (!.runtime_odin2_allowed()) {
    stop(
      "Compiled generator object not found (", obj, "). Runtime odin2 compilation is disabled.\n",
      "Enable it only for development via options(chlaa.allow_runtime_odin2 = TRUE) ",
      "or CHLAA_ALLOW_RUNTIME_ODIN2=true.",
      call. = FALSE
    )
  }

  if (!requireNamespace("odin2", quietly = TRUE)) {
    stop(
      "Compiled generator object not found (", obj, "). ",
      "Install odin2 to compile at runtime, or run odin2::odin_package('.') and rebuild.",
      call. = FALSE
    )
  }

  filename <- paste0(obj, ".R")
  path <- .find_odin_file(filename)
  if (!nzchar(path)) {
    stop("Could not find odin file '", filename, "'.", call. = FALSE)
  }

  # odin2::odin returns a dust_system_generator directly
  gen <- odin2::odin(path, input_type = "file", quiet = TRUE)
  .as_dust_generator(gen)
}

#' Simulate the cholera model
#'
#' @param pars Parameter list, typically from `chlaa_parameters()`.
#' @param time Vector of times to simulate at (days).
#' @param n_particles Number of particles.
#' @param dt Discrete time step (days).
#' @param seed RNG seed.
#' @param n_threads Threads for dust2.
#' @param deterministic Run in deterministic mode (replacing RNG draws with expectations) if supported.
#'
#' @return A data.frame with columns `time`, `particle`, and model variables.
#' @export
chlaa_simulate <- function(pars,
                            time,
                            n_particles = 1,
                            dt = 0.25,
                            seed = 1,
                            n_threads = 1,
                            deterministic = FALSE) {

  .check_named_list(pars, "pars")

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
  y <- dust2::dust_system_simulate(sys, time)
  st <- dust2::dust_unpack_state(sys, y)

  n_t <- length(time)
  out <- data.frame(
    time = rep(time, each = n_particles),
    particle = rep(seq_len(n_particles), times = n_t)
  )

  # dust_unpack_state returns matrices shaped (particle x time)
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

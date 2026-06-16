test_that("simulation runs with the bundled generator (fast smoke test)", {
  skip_if_not_installed("dust2")

  pars <- chlaa_parameters()
  time <- 0:2
  n_particles <- 2

  sim <- chlaa_simulate(pars, time = time, n_particles = n_particles, dt = 1, seed = 1)

  expect_true(is.data.frame(sim))
  expect_equal(nrow(sim), length(time) * n_particles)
  expect_true(all(c("time", "particle", "inc_symptoms", "cum_deaths") %in% names(sim)))

  # Basic sanity: no missing values in core outputs
  expect_false(anyNA(sim$inc_symptoms))
  expect_false(anyNA(sim$cum_deaths))
})

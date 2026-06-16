test_that("pmcmc runs with the bundled fit generator (tiny smoke test)", {
  skip_if_not_installed("dust2")
  skip_if_not_installed("monty")
  skip_on_cran()

  pars <- chlaa_parameters()

  # Keep cases at zero to ensure finite initial likelihood at default parameters.
  data <- data.frame(
    time = 1:3,
    cases = c(0L, 0L, 0L)
  )

  fit <- chlaa_fit_pmcmc(
    data = data,
    pars = pars,
    n_particles = 16,
    n_steps = 8,
    seed = 1,
    proposal_var = 0.01
  )

  expect_s3_class(fit, "chlaa_fit")

  draws <- chlaa_fit_draws(fit)
  expect_true(is.matrix(draws))
  expect_equal(nrow(draws), 8)
  expect_true(ncol(draws) > 0)
  expect_true(all(is.finite(draws)))
})

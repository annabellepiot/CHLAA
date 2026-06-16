test_that("chlaa_prepare_data validates and normalises input", {
  dat <- data.frame(t = c(1, 3, 2), y = c(1, 2, 3))
  out <- chlaa_prepare_data(dat, time_col = "t", cases_col = "y")
  expect_equal(names(out), c("time", "cases"))
  expect_equal(out$time, c(1, 2, 3))
  expect_equal(out$cases, c(1L, 3L, 2L))

  out2 <- chlaa_prepare_data(
    data.frame(time = c(1, 3), cases = c(2, 4)),
    expected_step = 1,
    fill_missing = TRUE
  )
  expect_equal(out2$time, c(1, 2, 3))
  expect_equal(out2$cases, c(2L, 0L, 4L))

  expect_error(chlaa_prepare_data(data.frame(time = c(1, 1), cases = c(0, 1))))
  expect_error(chlaa_prepare_data(data.frame(time = 1:3, cases = c(1, -1, 2))))
})

test_that("fit report and trace helpers return expected structures", {
  skip_if_not_installed("dust2")
  skip_if_not_installed("monty")

  pars <- chlaa_parameters()
  data <- data.frame(time = 1:3, cases = c(0L, 0L, 0L))
  fit <- chlaa_fit_pmcmc(
    data = data,
    pars = pars,
    n_particles = 16,
    n_steps = 10,
    seed = 1,
    proposal_var = 0.01
  )

  tr <- chlaa_fit_trace(fit)
  expect_true(is.data.frame(tr))
  expect_true(all(c("iteration", "parameter", "value") %in% names(tr)))

  rpt <- chlaa_fit_report(fit, burnin = 0, thin = 1)
  expect_true(is.list(rpt))
  expect_true(is.numeric(rpt$acceptance_rate))
  expect_true(rpt$acceptance_rate >= 0 && rpt$acceptance_rate <= 1)
  expect_true(is.data.frame(rpt$posterior_summary))
})

test_that("pmcmc infers daily observation interval and start time", {
  skip_if_not_installed("dust2")
  skip_if_not_installed("monty")

  pars <- chlaa_parameters()
  data <- data.frame(time = 0:3, cases = c(0L, 0L, 0L, 0L))
  fit <- chlaa_fit_pmcmc(
    data = data,
    pars = pars,
    n_particles = 16,
    n_steps = 8,
    seed = 1,
    proposal_var = 0.01
  )

  expect_equal(attr(fit, "obs_interval"), 1)
  expect_equal(attr(fit, "time_start"), -1)
  expect_equal(unique(attr(fit, "data")$obs_interval), 1)
})

test_that("pmcmc infers filter start from weekly observation spacing", {
  skip_if_not_installed("dust2")
  skip_if_not_installed("monty")

  pars <- chlaa_parameters()
  data <- data.frame(time = c(7, 14, 21), cases = c(0L, 0L, 0L))
  fit <- chlaa_fit_pmcmc(
    data = data,
    pars = pars,
    n_particles = 16,
    n_steps = 8,
    seed = 1,
    proposal_var = 0.01
  )

  expect_equal(attr(fit, "obs_interval"), 7)
  expect_equal(attr(fit, "time_start"), 0)
  expect_equal(unique(attr(fit, "data")$obs_interval), 7)
})

test_that("pmcmc allows explicit observation interval for one data point", {
  skip_if_not_installed("dust2")
  skip_if_not_installed("monty")

  pars <- chlaa_parameters()
  fit <- chlaa_fit_pmcmc(
    data = data.frame(time = 7, cases = 0L),
    pars = pars,
    n_particles = 16,
    n_steps = 8,
    seed = 1,
    proposal_var = 0.01,
    obs_interval = 7
  )

  expect_equal(attr(fit, "obs_interval"), 7)
  expect_equal(attr(fit, "time_start"), 0)
})

test_that("pmcmc supports multiple chains with per-chain starting parameters", {
  skip_if_not_installed("dust2")
  skip_if_not_installed("monty")

  pars <- chlaa_parameters()
  chain_2 <- pars
  chain_2$reporting_rate <- 0.25

  fit <- chlaa_fit_pmcmc(
    data = data.frame(time = 0:2, cases = c(0L, 0L, 0L)),
    pars = pars,
    chain_pars = list(pars, chain_2),
    n_particles = 1,
    n_steps = 6,
    n_chains = 2,
    seed = 1,
    proposal_var = 0.001,
    deterministic = TRUE
  )

  expect_equal(dim(fit$pars), c(9L, 6L, 2L))
  expect_equal(attr(fit, "chain_pars")[[2]]$reporting_rate, 0.25)
  expect_equal(unname(fit$initial["reporting_rate", 2]), 0.25)

  draws <- chlaa_fit_draws(fit)
  expect_equal(nrow(draws), 12)

  trace <- chlaa_fit_trace(fit, burnin = 0)
  expect_true(all(c("chain", "iteration", "parameter", "value") %in% names(trace)))
  expect_equal(sort(unique(trace$chain)), c("chain_1", "chain_2"))

  rpt <- chlaa_fit_report(fit, burnin = 0)
  expect_true(is.data.frame(rpt$acceptance_by_chain))
  expect_equal(nrow(rpt$acceptance_by_chain), 2)
  expect_true(is.numeric(rpt$acceptance_rate))
})

test_that("pmcmc validates chain_pars and n_chains", {
  skip_if_not_installed("dust2")
  skip_if_not_installed("monty")

  pars <- chlaa_parameters()
  dat <- data.frame(time = 0:2, cases = c(0L, 0L, 0L))

  expect_error(
    chlaa_fit_pmcmc(
      data = dat,
      pars = pars,
      chain_pars = list(pars),
      n_particles = 1,
      n_steps = 2,
      n_chains = 2,
      deterministic = TRUE
    ),
    "chain_pars"
  )

  expect_error(
    chlaa_fit_pmcmc(
      data = dat,
      pars = pars,
      n_particles = 1,
      n_steps = 2,
      n_chains = 1.5,
      deterministic = TRUE
    ),
    "n_chains"
  )
})

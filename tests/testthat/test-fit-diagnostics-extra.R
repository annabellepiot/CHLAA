test_that("trace helpers can unpack transformed posterior draws to natural scale", {
  skip_if_not_installed("monty")

  pars <- chlaa_parameters()
  pars$log_trans_prob <- log(pars$trans_prob)
  pars$logit_reporting_rate <- qlogis(pars$reporting_rate)
  pars$log_obs_size <- log(pars$obs_size)

  fit_names <- c("log_trans_prob", "logit_reporting_rate", "log_obs_size")
  fixed <- pars[setdiff(names(pars), c(fit_names, "trans_prob", "reporting_rate", "obs_size"))]
  packer <- monty::monty_packer(
    scalar = fit_names,
    fixed = fixed,
    process = function(p) {
      list(
        trans_prob = exp(p$log_trans_prob),
        reporting_rate = plogis(p$logit_reporting_rate),
        obs_size = exp(p$log_obs_size)
      )
    }
  )

  draws <- array(
    c(
      log(8e-4), qlogis(0.30), log(40),
      log(1e-3), qlogis(0.40), log(80),
      log(9e-4), qlogis(0.35), log(50),
      log(1.1e-3), qlogis(0.45), log(90)
    ),
    dim = c(3, 2, 2),
    dimnames = list(fit_names, NULL, c("chain_1", "chain_2"))
  )
  fit <- list(pars = draws)
  class(fit) <- "chlaa_fit"
  attr(fit, "packer") <- packer

  tr <- chlaa_fit_trace(fit, burnin = 0, scale = "natural")
  expect_true(all(c("trans_prob", "reporting_rate", "obs_size") %in% unique(tr$parameter)))
  expect_false(any(fit_names %in% unique(tr$parameter)))
  expect_equal(sort(unique(tr$chain)), c("chain_1", "chain_2"))
})

test_that("parameter distribution and pair plots accept truth overlays", {
  skip_if_not_installed("monty")
  skip_if_not_installed("ggplot2")

  pars <- chlaa_parameters()
  fit_names <- c("trans_prob", "reporting_rate", "obs_size")
  fixed <- pars[setdiff(names(pars), fit_names)]
  packer <- monty::monty_packer(scalar = fit_names, fixed = fixed)

  draws <- array(
    c(
      8e-4, 0.30, 40,
      9e-4, 0.35, 60,
      1e-3, 0.40, 80,
      1.1e-3, 0.45, 100
    ),
    dim = c(3, 2, 2),
    dimnames = list(fit_names, NULL, c("chain_1", "chain_2"))
  )
  fit <- list(pars = draws)
  class(fit) <- "chlaa_fit"
  attr(fit, "packer") <- packer

  truth <- c(trans_prob = 9e-4, reporting_rate = 0.35, obs_size = 60)
  p1 <- chlaa_plot_parameter_distributions(fit, parameters = fit_names, burnin = 0, truth = truth)
  p2 <- chlaa_plot_parameter_pairs(fit, parameters = fit_names, burnin = 0, truth = truth)

  expect_true(inherits(p1, "ggplot"))
  expect_true(inherits(p2, "ggplot") || inherits(p2, "ggmatrix"))
})

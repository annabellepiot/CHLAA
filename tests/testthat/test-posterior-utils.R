test_that("chlaa_update_from_fit updates known parameters", {
  pars <- chlaa_parameters()

  fake_draws <- matrix(
    c(
      0.06, 12.0, 4.0,
      0.05, 10.0, 5.0
    ),
    nrow = 2, byrow = TRUE
  )
  colnames(fake_draws) <- c("trans_prob", "contact_rate", "incubation_time")

  fit <- list(pars = fake_draws)
  class(fit) <- "chlaa_fit"

  upd <- chlaa_update_from_fit(fit, pars, draw = "mean", burnin = 0, thin = 1, validate = TRUE)

  expect_true(is.list(upd))
  expect_equal(upd$trans_prob, mean(c(0.06, 0.05)))
  expect_equal(upd$contact_rate, mean(c(12.0, 10.0)))
  expect_equal(upd$incubation_time, mean(c(4.0, 5.0)))
})

test_that("chlaa_update_from_fit respects transformed packer process", {
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

  fake_draws <- matrix(
    c(
      log(8e-4), qlogis(0.30), log(40),
      log(1e-3), qlogis(0.40), log(80)
    ),
    nrow = 2,
    byrow = TRUE
  )
  colnames(fake_draws) <- fit_names

  fit <- list(pars = fake_draws)
  class(fit) <- "chlaa_fit"
  attr(fit, "packer") <- packer

  upd <- chlaa_update_from_fit(fit, pars, draw = "mean", burnin = 0, thin = 1, validate = TRUE)

  theta <- colMeans(fake_draws)
  expect_equal(upd$trans_prob, exp(theta[["log_trans_prob"]]))
  expect_equal(upd$reporting_rate, plogis(theta[["logit_reporting_rate"]]))
  expect_equal(upd$obs_size, exp(theta[["log_obs_size"]]))
})

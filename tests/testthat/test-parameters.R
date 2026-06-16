test_that("chlaa_parameters provides required names", {
  pars <- chlaa_parameters(validate = FALSE)
  expect_true(is.list(pars))
  expect_true(!is.null(names(pars)))
  req <- chlaa_parameter_info()$name
  expect_true(all(req %in% names(pars)))
})

test_that("chlaa_parameters_validate catches missing parameters", {
  pars <- chlaa_parameters(validate = FALSE)
  pars$N <- NULL
  expect_error(chlaa_parameters_validate(pars))
})

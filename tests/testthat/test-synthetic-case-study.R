test_that("chlaa_generate_example_outbreak_data returns coherent data", {
  dat <- chlaa_generate_example_outbreak_data(time = 0:120, seed = 1, n_particles = 4)

  expect_true(is.data.frame(dat))
  expect_true(all(c("date", "time", "cases", "mu_cases", "inc_symptoms_truth", "inc_infections_truth") %in% names(dat)))
  expect_equal(dat$time, 0:120)
  expect_true(all(dat$cases >= 0))

  truth <- attr(dat, "truth_parameters", exact = TRUE)
  expect_true(is.list(truth))
  expect_true(all(c("N", "trans_prob", "reporting_rate", "obs_size") %in% names(truth)))
})

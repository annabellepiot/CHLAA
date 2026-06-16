test_that("standard scenarios include required templates", {
  pars <- chlaa_parameters()
  sc <- chlaa_standard_scenarios(pars, trigger_time = 30, horizon = 180)
  nms <- vapply(sc, `[[`, character(1), "name")
  expect_true(all(c(
    "baseline", "no_vaccination", "early_vaccination",
    "late_vaccination", "wash_only", "cm_only", "combined_package"
  ) %in% nms))
})

test_that("scenario summary and report return expected outputs", {
  skip_if_not_installed("dust2")
  pars <- chlaa_parameters()
  sc <- chlaa_standard_scenarios(pars, trigger_time = 20, horizon = 100, vax_total_doses = 20000)
  runs <- chlaa_run_scenarios(pars, sc, time = 0:30, n_particles = 3, dt = 1, seed = 1)

  summ <- chlaa_scenario_summary(runs, baseline = "baseline")
  expect_true(is.data.frame(summ))
  expect_true(all(c("scenario", "total_cases", "total_deaths", "cases_averted", "deaths_averted") %in% names(summ)))

  rep <- chlaa_scenario_report(runs, baseline = "baseline", include_econ = TRUE, wtp = 1000)
  expect_true(is.list(rep))
  expect_true(is.data.frame(rep$summary))
  expect_true(is.data.frame(rep$comparison))
  expect_true(is.list(rep$plots))
})

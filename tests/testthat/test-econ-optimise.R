test_that("economics defaults load and can be overridden", {
  d <- chlaa_econ_defaults()
  expect_true(is.list(d))
  expect_true(all(c("cost_per_vaccine_dose", "dw_symptomatic", "yll_per_death") %in% names(d)))

  d2 <- chlaa_econ_defaults(overrides = list(cost_per_vaccine_dose = 3.5))
  expect_equal(d2$cost_per_vaccine_dose, 3.5)
})

test_that("ce outputs include NMB and CEAC", {
  skip_if_not_installed("dust2")
  pars <- chlaa_parameters()
  sc <- list(
    chlaa_scenario("baseline", list()),
    chlaa_scenario("intervention", list(chlor_start = 0, chlor_end = 30, chlor_effect = 0.2))
  )
  runs <- chlaa_run_scenarios(pars, sc, time = 0:25, n_particles = 4, dt = 1, seed = 1)

  cmp <- chlaa_compare_scenarios(runs, baseline = "baseline", include_econ = TRUE, wtp = 1000)
  expect_true(all(c("nmb", "inmb", "cost_diff", "dalys_averted") %in% names(cmp)))

  ce <- chlaa_ceac(runs, baseline = "baseline", wtp = c(0, 500, 1000))
  expect_true(is.data.frame(ce))
  expect_true(all(c("scenario", "wtp", "prob_best") %in% names(ce)))
})

test_that("optimiser respects constraints and returns evaluations", {
  skip_if_not_installed("dust2")
  pars <- chlaa_parameters()

  opt <- chlaa_optimise_budget(
    pars = pars,
    budget = 100000,
    time = 0:20,
    n_particles = 5,
    dt = 1,
    grid_size = 6,
    min_fraction = list(vax = 0.1),
    max_fraction = list(wash = 0.6),
    max_vax_doses_per_day = 1000,
    method = "grid"
  )

  expect_true(is.list(opt))
  expect_true(is.data.frame(opt$best))
  expect_true(is.data.frame(opt$evaluations))
  expect_true(all(opt$evaluations$frac_vax >= 0.1))
  expect_true(all(opt$evaluations$frac_wash <= 0.6 + 1e-8))
})

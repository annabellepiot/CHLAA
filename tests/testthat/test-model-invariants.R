test_that("simulation preserves core state invariants", {
  skip_if_not_installed("dust2")

  pars <- chlaa_parameters()
  sim <- chlaa_simulate(pars, time = 0:20, n_particles = 3, dt = 1, seed = 1)

  compartments <- c("S", "E", "A", "M", "Sev", "Mu", "Mt", "Sevu", "Sevt", "Ra", "Rs", "V1", "V2", "Du", "Dt")
  expect_true(all(compartments %in% names(sim)))

  for (nm in compartments) {
    expect_true(all(sim[[nm]] >= 0), info = nm)
  }

  total <- rowSums(sim[compartments])
  expect_equal(total, rep(pars$N, length(total)))
})

test_that("daily incidence and cumulative counters are consistent", {
  skip_if_not_installed("dust2")

  pars <- chlaa_parameters()
  sim <- chlaa_simulate(pars, time = 0:20, n_particles = 2, dt = 1, seed = 2)

  map <- c(
    inc_infections = "cum_infections",
    inc_symptoms = "cum_symptoms",
    inc_deaths = "cum_deaths",
    inc_vax1 = "cum_vax1",
    inc_vax2 = "cum_vax2"
  )

  by_particle <- split(sim, sim$particle)
  for (df in by_particle) {
    df <- df[order(df$time), , drop = FALSE]
    for (inc in names(map)) {
      cum <- map[[inc]]
      expect_true(all(df[[inc]] >= 0), info = inc)
      expect_true(all(diff(df[[cum]]) >= 0), info = cum)
      expect_equal(diff(df[[cum]]), df[[inc]][-1], info = paste(inc, cum))
    }
  }
})

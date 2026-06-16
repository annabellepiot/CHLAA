# Deterministic_workflow.R — Step 1: find good PMCMC starting values
#
# Uses chlaa_simulate(deterministic = TRUE) to run the model without
# stochasticity, then optimises trans_prob, reporting_rate, seek_severe
# against weekly IDSR data for Kirotshe via optim().

library(chlaa)
library(tidyverse)

data_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/CHLAA/chlaa/analysis/data"
data_dir <- "analysis/data"
# ---- 1. Load & prepare data ----

hz_params <- read.csv(file.path(data_dir, "hz_parameters.csv"),
  stringsAsFactors = FALSE
) %>%
  pivot_wider(names_from = parameter, values_from = value)

idsr <- read.csv(file.path(data_dir, "IDSR_dataset.csv"))

i <- which(hz_params$hz == "kirotshe")
outbreak_start <- as.Date(hz_params$outbreak_start[i])
outbreak_end <- as.Date(hz_params$outbreak_end[i])

hz_weekly <- idsr %>%
  filter(hz == "kirotshe") %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= outbreak_start, date <= outbreak_end) %>%
  arrange(date) %>%
  mutate(time = seq_len(n()) * 7L)

pop_hz <- hz_weekly$population[1]

# Initial seeding (same logic as Case_Study)
seed_row <- idsr %>%
  filter(hz == "kirotshe") %>%
  mutate(date = as.Date(date)) %>%
  filter(date <= outbreak_start - 14) %>%
  arrange(desc(date)) %>%
  slice(1)

# Scale observed cases by reporting rate to estimate true unobserved burden
expected_reporting_rate <- 0.0282
E0_val <- if (nrow(seed_row) > 0 && seed_row$cases > 0) {
  ceiling(seed_row$cases / expected_reporting_rate)
} else {
  ceiling(max(1, hz_weekly$cases[1]) / expected_reporting_rate)
}

# ---- 2. Convert intervention dates to day offsets ----

to_day <- function(x) {
  d <- as.Date(x)
  if (is.na(d)) 0L else as.integer(d - outbreak_start)
}

# ---- 3. Base (fixed) parameters ----

pars_base <- chlaa_parameters(
  N = pop_hz, E0 = E0_val, Sev0 = 0, M0 = 0,
  immunity_asym = 280,
  contact_rate = 0,
  incubation_time = 4.845,
  duration_sym = 14.48,
  seek_mild = 0.1,
  fatality_treated = 0.001,
  fatality_untreated = 0.0043,
  vax2_doses_per_day = 0, vax2_total_doses = 0,
  orc_start = to_day(hz_params$orc_start[i]),
  orc_end = to_day(hz_params$orc_end[i]),
  ctc_start = to_day(hz_params$ctc_start[i]),
  ctc_end = to_day(hz_params$ctc_end[i]),
  hyg_start = to_day(hz_params$hyg_start[i]),
  hyg_end = to_day(hz_params$hyg_end[i]),
  hyg_effect = as.numeric(hz_params$hyg_effect[i]),
  cati_start = to_day(hz_params$cati_start[i]),
  cati_end = to_day(hz_params$cati_end[i]),
  cati_effect = as.numeric(hz_params$cati_effect[i]),
  chlor_start = to_day(hz_params$chlor_start[i]),
  chlor_end = to_day(hz_params$chlor_end[i]),
  chlor_effect = as.numeric(hz_params$chlor_effect[i])
)

# ---- 4. Objective: Poisson negative log-likelihood ----
# With zero_every = 7, inc_symptoms accumulates over each 7-day window.
# Simulate at weekly time points and read the accumulator directly.

sim_times <- hz_weekly$time

nll_fn <- function(theta) {
  pars_try <- pars_base
  pars_try$trans_prob <- theta[1]
  pars_try$reporting_rate <- theta[2]
  pars_try$seek_severe <- theta[3]

  sim <- tryCatch(
    chlaa_simulate(pars_try, time = sim_times, n_particles = 1,
                   deterministic = TRUE, seed = 1),
    error = function(e) NULL
  )
  if (is.null(sim)) return(1e10)

  predicted <- pmax(sim$inc_symptoms * pars_try$reporting_rate, 1e-6)
  -sum(dpois(hz_weekly$cases, lambda = predicted, log = TRUE))
}

# ---- 5. Optimise (L-BFGS-B with Hessian) ----

fit <- optim(
  par     = c(0.001, 0.20, 0.50),
  fn      = nll_fn,
  method  = "L-BFGS-B",
  lower   = c(1e-6, 0.01, 0.01),
  upper   = c(0.1,  1.00, 1.00),
  hessian = TRUE
)

best <- list(
  trans_prob     = fit$par[1],
  reporting_rate = fit$par[2],
  seek_severe    = fit$par[3]
)

# ---- 6. Extract MLE starting values & proposal covariance matrix ----

vcv <- solve(fit$hessian)
d <- length(fit$par)
proposal_matrix <- vcv * (2.38^2 / d)

# Fix any negative diagonal entries (variance must be positive).
# Replace with (10% of MLE)^2 — a safe step size for exploration.
for (k in seq_len(d)) {
  if (proposal_matrix[k, k] <= 0) {
    proposal_matrix[k, k] <- (0.1 * fit$par[k])^2
  }
}

cat("\n--- MLE starting values ---\n")
cat(sprintf("  trans_prob      = %.6f\n  reporting_rate  = %.4f\n  seek_severe     = %.4f\n",
            best$trans_prob, best$reporting_rate, best$seek_severe))
cat(sprintf("  NLL             = %.1f\n", fit$value))
cat(sprintf("  convergence     = %d (0 = success)\n", fit$convergence))
cat("\n--- Proposal covariance matrix (Gelman-scaled) ---\n")
print(proposal_matrix)

# ---- 7. Plot fit vs data ----

pars_best <- pars_base
pars_best$trans_prob <- best$trans_prob
pars_best$reporting_rate <- best$reporting_rate
pars_best$seek_severe <- best$seek_severe

sim_best <- chlaa_simulate(pars_best, time = sim_times,
                           n_particles = 1, deterministic = TRUE, seed = 1)

sim_plot <- data.frame(
  date = hz_weekly$date,
  predicted = sim_best$inc_symptoms * best$reporting_rate
)

p <- ggplot() +
  geom_col(data = hz_weekly, aes(date, cases), fill = "grey70", width = 6, alpha = 0.6) +
  geom_line(data = sim_plot, aes(date, predicted), colour = "#2c7f62", linewidth = 0.8) +
  geom_point(data = hz_weekly, aes(date, cases), size = 1.5) +
  labs(title = "Kirotshe — deterministic calibration",
       subtitle = sprintf("trans_prob=%.5f, reporting_rate=%.3f, seek_severe=%.3f",
                          best$trans_prob, best$reporting_rate, best$seek_severe),
       x = NULL, y = "Cases/week") +
  theme_minimal(base_size = 12)

fig_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/CHLAA/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(fig_dir, "kirotshe_deterministic_calibration.png"),
       p, width = 10, height = 5, dpi = 150)
cat("\nPlot saved to:", file.path(fig_dir, "kirotshe_deterministic_calibration.png"), "\n")

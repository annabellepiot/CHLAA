## odin2 cholera outbreak model (stochastic, discrete-time)
##
## Notes:
## - Designed for compilation with odin2 -> dust2 generator.


# Parameters


N <- parameter(540000)

# Initial conditions
E0 <- parameter(0)
A0 <- parameter(0)
M0 <- parameter(0)
Sev0 <- parameter(2)
Mu0 <- parameter(0)
Mt0 <- parameter(0)
Sevu0 <- parameter(0)
Sevt0 <- parameter(0)
Ra0 <- parameter(0)
Rs0 <- parameter(0)
V10 <- parameter(0)
V20 <- parameter(0)
Du0 <- parameter(0)
Dt0 <- parameter(0)
C0 <- parameter(0.0)

# Natural history
prop_asym <- parameter(0.75)
incubation_time <- parameter(4.845)
duration_asym <- parameter(5.0)
duration_sym <- parameter(14.48)
time_to_next_stage <- parameter(1.0)
p_progress_severe <- parameter(0.30)
immunity_asym <- parameter(180.0)
immunity_sym <- parameter(1095.0)

# Transmission + environment
beta_p2p <- parameter(0.05)
trans_prob <- parameter(0.055)

# Stochastic drift parameters
drift_volatility <- parameter(0.05) # sigma: How high the spikes can jump
drift_reversion <- parameter(0.05)  # theta: How fast the spikes are pulled back to baseline (half-life ~14 days)
time_to_contaminate <- parameter(19.075)
water_clearance_time <- parameter(30.0)
contam_half_sat <- parameter(1.0)
shed_asym <- parameter(90.69e3)
shed_mild <- parameter(9.5005e6)
shed_severe <- parameter(32.945e6)
contam_scale <- parameter(1.0e10)

# Care and case management
seek_mild <- parameter(0.30)
seek_severe <- parameter(0.68)
orc_capacity <- parameter(500.0)
ctc_capacity <- parameter(100.0)
treated_shed_mult_orc <- parameter(0.5)
treated_shed_mult_ctc <- parameter(0.0)
fatality_treated <- parameter(0.0021)
fatality_untreated <- parameter(0.5)

# Intervention windows and effects (step functions)
chlor_start <- parameter(0.0)
chlor_end <- parameter(0.0)
chlor_effect <- parameter(0.0)

hyg_start <- parameter(0.0)
hyg_end <- parameter(0.0)
hyg_effect <- parameter(0.0)

lat_start <- parameter(0.0)
lat_end <- parameter(0.0)
lat_effect <- parameter(0.0)

cati_start <- parameter(0.0)
cati_end <- parameter(0.0)
cati_effect <- parameter(0.0)

orc_start <- parameter(0.0)
orc_end <- parameter(0.0)
ctc_start <- parameter(0.0)
ctc_end <- parameter(0.0)

# Vaccination campaigns
vax1_start <- parameter(0.0)
vax1_end <- parameter(0.0)
vax1_total_doses <- parameter(0.0)

vax2_start <- parameter(0.0)
vax2_end <- parameter(0.0)
vax2_total_doses <- parameter(0.0)

#Time-varying vaccination rates (interpolated from data)
vax1_doses_daily <- interpolate(vax1_schedule_time, vax1_schedule_doses, "constant")
vax2_doses_daily <- interpolate(vax2_schedule_time, vax2_schedule_doses, "constant")

#Data arrays for vaccination schedule
vax1_schedule_time <- parameter()
vax1_schedule_doses <- parameter()
dim(vax1_schedule_time) <- n_vax1_schedule
dim(vax1_schedule_doses) <- n_vax1_schedule
n_vax1_schedule <- parameter()

vax2_schedule_time <- parameter()
vax2_schedule_doses <- parameter()
dim(vax2_schedule_time) <- n_vax2_schedule
dim(vax2_schedule_doses) <- n_vax2_schedule
n_vax2_schedule <- parameter()


ve_1 <- parameter(0.4)
ve_2 <- parameter(0.7)
vax_immunity_1 <- parameter(180.0)
vax_immunity_2 <- parameter(1095.0)

# Derived "active" flags
chlor_active <- if (time >= chlor_start && time < chlor_end) 1.0 else 0.0
hyg_active <- if (time >= hyg_start && time < hyg_end) 1.0 else 0.0
lat_active <- if (time >= lat_start && time < lat_end) 1.0 else 0.0
cati_active <- if (time >= cati_start && time < cati_end) 1.0 else 0.0
orc_active <- if (time >= orc_start && time < orc_end) 1.0 else 0.0
ctc_active <- if (time >= ctc_start && time < ctc_end) 1.0 else 0.0

vax1_active <- if (time >= vax1_start && time < vax1_end) 1.0 else 0.0
vax2_active <- if (time >= vax2_start && time < vax2_end) 1.0 else 0.0

# Transmission modifier (bounded)
trans_mult <- max(0.0, 1.0 - (chlor_active * chlor_effect + hyg_active * hyg_effect + cati_active * cati_effect))
shed_mult <- max(0.0, 1.0 - (lat_active * lat_effect))

# State variables (integer-valued compartments; C is continuous)
initial(S) <- N - E0 - A0 - M0 - Sev0 - Mu0 - Mt0 - Sevu0 - Sevt0 - Ra0 - Rs0 - V10 - V20
initial(E) <- E0
initial(A) <- A0
initial(M) <- M0
initial(Sev) <- Sev0
initial(Mu) <- Mu0
initial(Mt) <- Mt0
initial(Sevu) <- Sevu0
initial(Sevt) <- Sevt0
initial(Ra) <- Ra0
initial(Rs) <- Rs0
initial(V1) <- V10
initial(V2) <- V20

initial(Du) <- Du0
initial(Dt) <- Dt0
initial(C) <- C0

# Stochastic transmission drift (OU mean-reverting process in log-space)
initial(log_trans_drift) <- 0.0
update(log_trans_drift) <- log_trans_drift - (drift_reversion * log_trans_drift * dt) + Normal(0.0, drift_volatility * sqrt(dt))
trans_drift <- exp(log_trans_drift)

# Incidence accumulators. The daily names are used by the scenario and
# historical daily-data workflows. Weekly names support IDSR fitting without
# dividing weekly counts into pseudo-daily observations.
initial(inc_infections, zero_every = 1) <- 0
initial(inc_symptoms, zero_every = 1) <- 0
initial(inc_deaths, zero_every = 1) <- 0
initial(inc_vax1, zero_every = 1) <- 0
initial(inc_vax2, zero_every = 1) <- 0

initial(inc_infections_weekly, zero_every = 7) <- 0
initial(inc_symptoms_weekly, zero_every = 7) <- 0
initial(inc_deaths_weekly, zero_every = 7) <- 0
initial(inc_vax1_weekly, zero_every = 7) <- 0
initial(inc_vax2_weekly, zero_every = 7) <- 0

# Cumulative outputs
initial(cum_infections) <- 0
initial(cum_symptoms) <- 0
initial(cum_deaths) <- Du0 + Dt0
initial(cum_vax1) <- 0
initial(cum_vax2) <- 0
initial(cum_orc_treated) <- 0
initial(cum_ctc_treated) <- 0

# Helper transition probabilities
p_EI <- 1.0 - exp(-dt / incubation_time)
p_rec_asym <- 1.0 - exp(-dt / duration_asym)
p_stage <- 1.0 - exp(-dt / time_to_next_stage)
p_leave_sym <- 1.0 - exp(-dt / duration_sym)

p_wane_Ra <- 1.0 - exp(-dt / immunity_asym)
p_wane_Rs <- 1.0 - exp(-dt / immunity_sym)
p_wane_V1 <- 1.0 - exp(-dt / vax_immunity_1)
p_wane_V2 <- 1.0 - exp(-dt / vax_immunity_2)

# Environmental force (saturating)
env_force <- trans_prob * (C / (C + contam_half_sat))

# Effective infectious pool for person-to-person (simple)
I_eff <- A + M + Sev + Mu + Mt + Sevu + Sevt
p2p_force <- beta_p2p * (I_eff / N)

lambda <- trans_drift * trans_mult * (p2p_force + env_force)
p_inf <- 1.0 - exp(-lambda * dt)

# New infections by group (susceptibility reduced in vaccinated)
new_E_S <- Binomial(S, p_inf)
new_E_V1 <- Binomial(V1, min(1.0, p_inf * (1.0 - ve_1)))
new_E_V2 <- Binomial(V2, min(1.0, p_inf * (1.0 - ve_2)))
new_E <- new_E_S + new_E_V1 + new_E_V2

# Progression from E to infectious
new_I <- Binomial(E, p_EI)
new_A <- Binomial(new_I, prop_asym)
new_symp <- new_I - new_A

# Split symptomatic into mild vs severe at onset (simple)
new_Sev <- Binomial(new_symp, p_progress_severe)
new_M <- new_symp - new_Sev

# Asymptomatic recovery
rec_A <- Binomial(A, p_rec_asym)

# Symptomatic triage/progression from pre-triage stages
prog_M <- Binomial(M, p_stage)
prog_Sev <- Binomial(Sev, p_stage)

# Care seeking and capacity constraints
seek_M <- Binomial(prog_M, seek_mild)
seek_Sev <- Binomial(prog_Sev, seek_severe)

orc_cap_step <- floor(orc_active * orc_capacity * dt)
ctc_cap_step <- floor(ctc_active * ctc_capacity * dt)

treat_orc <- min(seek_M, orc_cap_step)
treat_ctc <- min(seek_Sev, ctc_cap_step)

# Remaining progression goes untreated
to_Mu <- prog_M - treat_orc
to_Sevu <- prog_Sev - treat_ctc

# Leave treated/untreated symptomatic compartments
leave_Mu <- Binomial(Mu, p_leave_sym)
leave_Mt <- Binomial(Mt, p_leave_sym)

leave_Sevu <- Binomial(Sevu, p_leave_sym)
leave_Sevt <- Binomial(Sevt, p_leave_sym)

death_Sevu <- Binomial(leave_Sevu, fatality_untreated)
death_Sevt <- Binomial(leave_Sevt, fatality_treated)

rec_Sevu <- leave_Sevu - death_Sevu
rec_Sevt <- leave_Sevt - death_Sevt

# Waning immunity
wane_Ra <- Binomial(Ra, p_wane_Ra)
wane_Rs <- Binomial(Rs, p_wane_Rs)
wane_V1 <- Binomial(V1, p_wane_V1)
wane_V2 <- Binomial(V2, p_wane_V2)

# Vaccination administration (bounded by supply and eligible people)
# Use interpolated daily doses instead of constant rate
vax1_daily_doses <- if (vax1_active > 0) vax1_doses_daily else 0.0
vax1_cap_step <- floor(vax1_daily_doses * dt)
vax1_remain <- max(0, vax1_total_doses - cum_vax1)
vax1_admin <- min(S, min(vax1_remain, vax1_cap_step))

vax2_daily_doses <- if (vax2_active > 0) vax2_doses_daily else 0.0
vax2_cap_step <- floor(vax2_daily_doses * dt)
vax2_remain <- max(0, vax2_total_doses - cum_vax2)
vax2_admin <- min(V1, min(vax2_remain, vax2_cap_step))

# Environment update (simple input/output dynamics)
shed_cfu <- shed_mult * (
  shed_asym * A +
    shed_mild * (M + Mu) +
    shed_mild * treated_shed_mult_orc * Mt +
    shed_severe * (Sev + Sevu) +
    shed_severe * treated_shed_mult_ctc * Sevt
)
shed_index <- shed_cfu / contam_scale
dC <- (shed_index / max(1e-9, time_to_contaminate)) - (C / max(1e-9, water_clearance_time))

update(C) <- max(0.0, C + dt * dC)

# Updates for compartments
# max(0, ...) guards on S and V1 prevent negative compartments when
# independent binomial outflows (infection + vaccination + waning)
# occasionally exceed the compartment size in the same sub-step.
update(S) <- max(0, S - new_E_S - vax1_admin + wane_Ra + wane_Rs + wane_V1 + wane_V2)
update(E) <- E + new_E - new_I
update(A) <- A + new_A - rec_A

update(M) <- M + new_M - prog_M
update(Sev) <- Sev + new_Sev - prog_Sev

update(Mu) <- Mu + to_Mu - leave_Mu
update(Mt) <- Mt + treat_orc - leave_Mt

update(Sevu) <- Sevu + to_Sevu - leave_Sevu
update(Sevt) <- Sevt + treat_ctc - leave_Sevt

update(Ra) <- Ra + rec_A - wane_Ra
update(Rs) <- Rs + leave_Mu + leave_Mt + rec_Sevu + rec_Sevt - wane_Rs

update(V1) <- max(0, V1 + vax1_admin - vax2_admin - wane_V1 - new_E_V1)
update(V2) <- max(0, V2 + vax2_admin - wane_V2 - new_E_V2)

update(Du) <- Du + death_Sevu
update(Dt) <- Dt + death_Sevt

# Incidence and cumulative outputs
update(inc_infections) <- inc_infections + new_E
update(inc_symptoms) <- inc_symptoms + new_symp
update(inc_deaths) <- inc_deaths + death_Sevu + death_Sevt
update(inc_vax1) <- inc_vax1 + vax1_admin
update(inc_vax2) <- inc_vax2 + vax2_admin

update(inc_infections_weekly) <- inc_infections_weekly + new_E
update(inc_symptoms_weekly) <- inc_symptoms_weekly + new_symp
update(inc_deaths_weekly) <- inc_deaths_weekly + death_Sevu + death_Sevt
update(inc_vax1_weekly) <- inc_vax1_weekly + vax1_admin
update(inc_vax2_weekly) <- inc_vax2_weekly + vax2_admin

update(cum_infections) <- cum_infections + new_E
update(cum_symptoms) <- cum_symptoms + new_symp
update(cum_deaths) <- cum_deaths + death_Sevu + death_Sevt
update(cum_vax1) <- cum_vax1 + vax1_admin
update(cum_vax2) <- cum_vax2 + vax2_admin
update(cum_orc_treated) <- cum_orc_treated + treat_orc
update(cum_ctc_treated) <- cum_ctc_treated + treat_ctc


# Data and observation model (for filtering / likelihood)

# Observation model parameters
reporting_rate <- parameter(0.2)
obs_size <- parameter(25.0)
death_reporting_rate <- parameter(0.5)
obs_size_deaths <- parameter(5.0)

cases <- data()
deaths <- data()
obs_interval <- data()
obs_inc_symptoms <- if (obs_interval <= 1.5) inc_symptoms else inc_symptoms_weekly
obs_inc_deaths <- if (obs_interval <= 1.5) inc_deaths else inc_deaths_weekly
cases ~ NegativeBinomial(mu = reporting_rate * obs_inc_symptoms, size = obs_size)
deaths ~ NegativeBinomial(mu = death_reporting_rate * obs_inc_deaths, size = obs_size_deaths)

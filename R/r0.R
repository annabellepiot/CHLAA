# Basic reproduction number via the next-generation matrix (NGM)

#' Compute R0 from the model's next-generation matrix
#'
#' Derives the basic reproduction number analytically from the odin cholera
#' model (see `inst/odin/cholera_model.R`) using the van den Driessche &
#' Watmough (2002) next-generation-matrix framework, linearised at the
#' disease-free equilibrium (S = N, all infected compartments and the
#' environmental reservoir `C` at 0). This replaces the previous ad hoc
#' `R0 = frac_neff * K_R0 * trans_prob` formula, whose constant `K_R0` had no
#' in-repo derivation and which omitted `contam_half_sat` (so it could not
#' reflect the environmental transmission route that in practice carries all
#' transmission once `contact_rate = 0`, as used throughout the fitting
#' pipeline).
#'
#' The infected subsystem is `E, A, M, Sev, Mu, Mt, Sevu, Sevt, C` (the
#' environmental reservoir is included as in standard treatments of
#' water-borne/environmentally-transmitted pathogens, e.g. Tien & Earn 2010).
#' R0 is the spectral radius of `F %*% solve(V)`, where `F` collects new
#' infections (only the `E` row is non-zero: person-to-person contacts and
#' the linearised environmental force `trans_prob * C / contam_half_sat`) and
#' `V` collects all other transitions (progression, care-seeking, treatment,
#' recovery/death exits, and the environmental shedding/decay dynamics).
#' Interventions (chlorination, hygiene, CATI, latrines) and vaccination are
#' intentionally not included in `pars` used here: R0 is always evaluated at
#' baseline, i.e. as though `trans_mult = shed_mult = 1`, matching the
#' standard definition of the basic (intervention-free) reproduction number.
#'
#' @param pars A named parameter list (or list with vector-valued
#'   `trans_prob`/`N` for evaluating R0 at many posterior draws at once). Must
#'   contain `trans_prob`, `N`, `contact_rate`, `contam_half_sat`,
#'   `incubation_time`, `prop_asym`, `p_progress_severe`, `duration_asym`,
#'   `time_to_next_stage`, `seek_mild`, `seek_severe`, `duration_sym`,
#'   `contam_scale`, `time_to_contaminate`, `water_clearance_time`,
#'   `shed_asym`, `shed_mild`, `shed_severe`, `treated_shed_mult_orc`,
#'   `treated_shed_mult_ctc`.
#'
#' @return A numeric vector of R0 values, one per element of
#'   `trans_prob`/`N` (recycled to a common length); a single value if both
#'   are scalar.
#' @export
chlaa_r0 <- function(pars) {
  .check_named_list(pars, "pars")

  fixed_scalar <- c(
    "contact_rate", "contam_half_sat",
    "incubation_time", "prop_asym", "p_progress_severe",
    "duration_asym", "time_to_next_stage", "seek_mild", "seek_severe",
    "duration_sym", "contam_scale", "time_to_contaminate",
    "water_clearance_time", "shed_asym", "shed_mild", "shed_severe",
    "treated_shed_mult_orc", "treated_shed_mult_ctc"
  )
  required <- c("trans_prob", "N", fixed_scalar)
  missing <- setdiff(required, names(pars))
  if (length(missing) > 0) {
    stop("chlaa_r0: pars is missing required fields: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  for (nm in fixed_scalar) {
    v <- pars[[nm]]
    if (length(unique(v)) > 1) {
      stop("chlaa_r0: pars$", nm, " must be a single (scalar) value, not a vector of varying values", call. = FALSE)
    }
  }

  g <- function(nm) as.numeric(pars[[nm]][1])
  contact_rate <- g("contact_rate")
  contam_half_sat <- g("contam_half_sat")
  incubation_time <- g("incubation_time")
  prop_asym <- g("prop_asym")
  p_progress_severe <- g("p_progress_severe")
  duration_asym <- g("duration_asym")
  time_to_next_stage <- g("time_to_next_stage")
  seek_mild <- g("seek_mild")
  seek_severe <- g("seek_severe")
  duration_sym <- g("duration_sym")
  contam_scale <- g("contam_scale")
  time_to_contaminate <- g("time_to_contaminate")
  water_clearance_time <- g("water_clearance_time")
  shed_asym <- g("shed_asym")
  shed_mild <- g("shed_mild")
  shed_severe <- g("shed_severe")
  treated_shed_mult_orc <- g("treated_shed_mult_orc")
  treated_shed_mult_ctc <- g("treated_shed_mult_ctc")

  # Underlying continuous-time rates recovered from the odin model's
  # 1 - exp(-dt / timescale) discretisations (inst/odin/cholera_model.R).
  sigma <- 1 / incubation_time
  sigma_A <- sigma * prop_asym
  sigma_Sev <- sigma * (1 - prop_asym) * p_progress_severe
  sigma_M <- sigma * (1 - prop_asym) * (1 - p_progress_severe)
  gamma_A <- 1 / duration_asym
  k <- 1 / time_to_next_stage
  gamma_sym <- 1 / duration_sym
  c_in <- 1 / (contam_scale * time_to_contaminate)
  c_out <- 1 / water_clearance_time

  states <- c("E", "A", "M", "Sev", "Mu", "Mt", "Sevu", "Sevt", "C")
  n <- length(states)
  idx <- setNames(seq_len(n), states)

  V <- matrix(0, n, n, dimnames = list(states, states))
  V["E", "E"] <- sigma

  V["A", "E"] <- -sigma_A
  V["A", "A"] <- gamma_A

  V["M", "E"] <- -sigma_M
  V["M", "M"] <- k

  V["Sev", "E"] <- -sigma_Sev
  V["Sev", "Sev"] <- k

  V["Mu", "M"] <- -k * (1 - seek_mild)
  V["Mu", "Mu"] <- gamma_sym

  V["Mt", "M"] <- -k * seek_mild
  V["Mt", "Mt"] <- gamma_sym

  V["Sevu", "Sev"] <- -k * (1 - seek_severe)
  V["Sevu", "Sevu"] <- gamma_sym

  V["Sevt", "Sev"] <- -k * seek_severe
  V["Sevt", "Sevt"] <- gamma_sym

  V["C", "A"] <- -c_in * shed_asym
  V["C", "M"] <- -c_in * shed_mild
  V["C", "Mu"] <- -c_in * shed_mild
  V["C", "Mt"] <- -c_in * shed_mild * treated_shed_mult_orc
  V["C", "Sev"] <- -c_in * shed_severe
  V["C", "Sevu"] <- -c_in * shed_severe
  V["C", "Sevt"] <- -c_in * shed_severe * treated_shed_mult_ctc
  V["C", "C"] <- c_out

  Vinv <- solve(V)

  trans_prob <- as.numeric(pars$trans_prob)
  N <- as.numeric(pars$N)
  m <- max(length(trans_prob), length(N))
  trans_prob <- rep_len(trans_prob, m)
  N <- rep_len(N, m)

  new_infection_cols <- c("A", "M", "Sev", "Mu", "Mt", "Sevu", "Sevt")

  r0 <- vapply(seq_len(m), function(i) {
    F_i <- matrix(0, n, n, dimnames = list(states, states))
    F_i["E", new_infection_cols] <- contact_rate
    F_i["E", "C"] <- trans_prob[i] * N[i] / contam_half_sat
    K_i <- F_i %*% Vinv
    max(Mod(eigen(K_i, only.values = TRUE)$values))
  }, numeric(1))

  if (m == 1) r0[[1]] else r0
}

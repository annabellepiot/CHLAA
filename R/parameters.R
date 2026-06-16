# Parameters

.param <- function(name, default, units, description, source) {
  data.frame(
    name = as.character(name),
    default = default,
    units = as.character(units),
    description = as.character(description),
    source = as.character(source),
    stringsAsFactors = FALSE
  )
}

#' Parameter metadata table for the cholera model
#'
#' Returns metadata for all parameters required by the odin model and observation model.
#' Defaults here are aligned to the defaults in `inst/odin/cholera_model.R`
#'
#' @return A data.frame with columns: name, default, units, description, source.
#' @export
chlaa_parameter_info <- function() {
  src <- "Aligned to odin model defaults (review against paper and Table S1)"

  rows <- list(
    # Population and initial conditions
    .param("N", 540000, "persons", "Total population size.", src),
    .param("E0", 0L, "persons", "Initial exposed (latent).", src),
    .param("A0", 0L, "persons", "Initial asymptomatic infectious.", src),
    .param("M0", 0L, "persons", "Initial mild symptomatic (pre-triage stage).", src),
    .param("Sev0", 1L, "persons", "Initial severe symptomatic (pre-triage stage).", src),
    .param("Mu0", 0L, "persons", "Initial mild untreated.", src),
    .param("Mt0", 0L, "persons", "Initial mild treated.", src),
    .param("Sevu0", 0L, "persons", "Initial severe untreated.", src),
    .param("Sevt0", 0L, "persons", "Initial severe treated.", src),
    .param("Ra0", 0L, "persons", "Initial recovered after asymptomatic infection.", src),
    .param("Rs0", 0L, "persons", "Initial recovered after symptomatic infection.", src),
    .param("V10", 0L, "persons", "Initial vaccinated (1 dose protected).", src),
    .param("V20", 0L, "persons", "Initial vaccinated (2 dose protected).", src),
    .param("Du0", 0L, "persons", "Initial cumulative deaths (untreated).", src),
    .param("Dt0", 0L, "persons", "Initial cumulative deaths (treated).", src),
    .param("C0", 0.0, "index", "Initial environmental contamination state.", src),

    # Natural history
    .param("prop_asym", 0.75, "probability", "Proportion of infections asymptomatic.", src),
    .param("incubation_time", 4.845, "days", "Mean incubation time.", src),
    .param("duration_asym", 5.0, "days", "Duration of asymptomatic infectiousness.", src),
    .param("duration_sym", 14.48, "days", "Duration of symptomatic infection.", src),
    .param("time_to_next_stage", 1.0, "days", "Time to next symptomatic stage (triage/progression).", src),
    .param("p_progress_severe", 0.30, "probability", "Probability mild progresses to severe.", src),
    .param("immunity_asym", 280, "days", "Immunity duration after asymptomatic infection.", src),
    .param("immunity_sym", 1095.0, "days", "Immunity duration after symptomatic infection.", src),

    # Transmission and environment
    .param("contact_rate", 10.01, "contacts/person/day", "Effective contact rate.", src),
    .param("trans_prob", 0.127, "probability/contact", "Per-contact transmission probability.", src),
    .param("time_to_contaminate", 19.075, "days", "Time-scale for contamination dynamics.", src),
    .param("water_clearance_time", 30.0, "days", "Environmental clearance time.", src),
    .param("contam_half_sat", 1.0, "index", "Half-saturation constant for contamination effect.", src),
    .param("shed_asym", 90.69e3, "CFU/person/day", "Shedding rate asymptomatic.", src),
    .param("shed_mild", 9.5005e6, "CFU/person/day", "Shedding rate mild.", src),
    .param("shed_severe", 32.945e6, "CFU/person/day", "Shedding rate severe.", src),
    .param("contam_scale", 1.0e10, "CFU/index", "Scaling from CFU to contamination index.", src),

    # Care and case management
    .param("seek_mild", 0.1, "probability", "Care seeking probability (mild).", src),
    .param("seek_severe", 0.2, "probability", "Care seeking probability (severe).", src),
    .param("orc_capacity", 500.0, "persons/day", "ORC capacity while active.", src),
    .param("ctc_capacity", 100.0, "persons/day", "CTC capacity while active.", src),
    .param("treated_shed_mult_orc", 0.5, "multiplier", "Shedding multiplier in ORC.", src),
    .param("treated_shed_mult_ctc", 0.0, "multiplier", "Shedding multiplier in CTC.", src),
    .param("fatality_treated", 0.0021, "probability", "Fatality probability treated (severe).", src),
    .param("fatality_untreated", 0.5, "probability", "Fatality probability untreated (severe).", src),

    # Intervention windows and effects
    .param("chlor_start", 0.0, "day", "Start time for chlorination.", src),
    .param("chlor_end", 0.0, "day", "End time for chlorination.", src),
    .param("chlor_effect", 0.0, "fraction", "Effect size for chlorination (transmission reduction).", src),

    .param("hyg_start", 0.0, "day", "Start time for hygiene promotion.", src),
    .param("hyg_end", 0.0, "day", "End time for hygiene promotion.", src),
    .param("hyg_effect", 0.0, "fraction", "Effect size for hygiene (transmission reduction).", src),

    .param("lat_start", 0.0, "day", "Start time for latrine intervention.", src),
    .param("lat_end", 0.0, "day", "End time for latrine intervention.", src),
    .param("lat_effect", 0.0, "fraction", "Effect size for latrines (shedding reduction).", src),

    .param("cati_start", 0.0, "day", "Start time for CATI.", src),
    .param("cati_end", 0.0, "day", "End time for CATI.", src),
    .param("cati_effect", 0.0, "fraction", "Effect size for CATI (transmission reduction).", src),

    .param("orc_start", 0.0, "day", "Start time for ORC availability.", src),
    .param("orc_end", 0.0, "day", "End time for ORC availability.", src),
    .param("ctc_start", 0.0, "day", "Start time for CTC availability.", src),
    .param("ctc_end", 0.0, "day", "End time for CTC availability.", src),

    # Vaccination
    .param("vax1_start", 0.0, "day", "Start time dose 1 campaign.", src),
    .param("vax1_end", 0.0, "day", "End time dose 1 campaign.", src),
    .param("vax1_doses_per_day", 0.0, "doses/day", "Dose 1 delivery rate.", src),
    .param("vax1_total_doses", 0.0, "doses", "Total dose 1 supply.", src),

    .param("vax2_start", 0.0, "day", "Start time dose 2 campaign.", src),
    .param("vax2_end", 0.0, "day", "End time dose 2 campaign.", src),
    .param("vax2_doses_per_day", 0.0, "doses/day", "Dose 2 delivery rate.", src),
    .param("vax2_total_doses", 0.0, "doses", "Total dose 2 supply.", src),

    .param("ve_1", 0.4, "fraction", "Efficacy after dose 1 (susceptibility reduction).", src),
    .param("ve_2", 0.7, "fraction", "Efficacy after dose 2 (susceptibility reduction).", src),
    .param("vax_immunity_1", 180.0, "days", "Protection duration after dose 1.", src),
    .param("vax_immunity_2", 1095.0, "days", "Protection duration after dose 2.", src),

    # Observation model
    .param("reporting_rate", 0.2, "fraction", "Reporting fraction for observed cases.", src),
    .param("obs_size", 25.0, "size", "Negative binomial size parameter.", src)
  )

  do.call(rbind, rows)
}

#' Create a parameter list for the cholera model
#'
#' Uses defaults from `chlaa_parameter_info()` and allows overriding via `...`.
#'
#' @param ... Named parameter overrides, e.g. `Sev0 = 3`, `N = 250000`.
#' @param validate Logical; validate the resulting parameter list.
#'
#' @return A named list of parameters.
#' @export
chlaa_parameters <- function(..., validate = TRUE) {
  info <- chlaa_parameter_info()
  out <- as.list(info$default)
  names(out) <- info$name

  overrides <- list(...)
  if (length(overrides) > 0) {
    if (is.null(names(overrides)) || any(names(overrides) == "")) {
      stop("All overrides must be named", call. = FALSE)
    }
    unknown <- setdiff(names(overrides), names(out))
    if (length(unknown) > 0) {
      stop("Unknown parameters in overrides: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    out[names(overrides)] <- overrides
  }

  if (isTRUE(validate)) chlaa_parameters_validate(out)
  out
}

#' Validate a cholera parameter list
#'
#' @param pars Parameter list.
#'
#' @return Invisibly returns `pars` if valid, otherwise errors.
#' @export
chlaa_parameters_validate <- function(pars) {
  .check_named_list(pars, "pars")
  req <- chlaa_parameter_info()$name
  missing <- setdiff(req, names(pars))
  if (length(missing) > 0) {
    stop("Missing parameters: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  prob_pars <- c(
    "prop_asym", "p_progress_severe",
    "seek_mild", "seek_severe",
    "chlor_effect", "hyg_effect", "lat_effect", "cati_effect",
    "ve_1", "ve_2", "reporting_rate",
    "fatality_treated", "fatality_untreated"
  )
  for (nm in intersect(prob_pars, names(pars))) {
    v <- pars[[nm]]
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < 0 || v > 1) {
      stop("Parameter ", nm, " must be in [0, 1]", call. = FALSE)
    }
  }

  nonneg <- c(
    "N","E0","A0","M0","Sev0","Mu0","Mt0","Sevu0","Sevt0","Ra0","Rs0","V10","V20","Du0","Dt0",
    "C0","incubation_time","duration_asym","duration_sym","time_to_next_stage",
    "immunity_asym","immunity_sym","contact_rate","trans_prob","time_to_contaminate",
    "water_clearance_time","contam_half_sat","shed_asym","shed_mild","shed_severe",
    "contam_scale","orc_capacity","ctc_capacity",
    "vax1_doses_per_day","vax1_total_doses","vax2_doses_per_day","vax2_total_doses",
    "vax_immunity_1","vax_immunity_2","obs_size"
  )
  for (nm in intersect(nonneg, names(pars))) {
    v <- pars[[nm]]
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < 0) {
      stop("Parameter ", nm, " must be non-negative", call. = FALSE)
    }
  }

  # Time window sanity: end >= start
  time_pairs <- list(
    c("chlor_start", "chlor_end"),
    c("hyg_start", "hyg_end"),
    c("lat_start", "lat_end"),
    c("cati_start", "cati_end"),
    c("orc_start", "orc_end"),
    c("ctc_start", "ctc_end"),
    c("vax1_start", "vax1_end"),
    c("vax2_start", "vax2_end")
  )
  for (p in time_pairs) {
    a <- as.numeric(pars[[p[1]]])
    b <- as.numeric(pars[[p[2]]])
    if (is.finite(a) && is.finite(b) && b < a) {
      stop("Invalid time window: ", p[1], " > ", p[2], call. = FALSE)
    }
  }

  invisible(pars)
}

#' Print a parameter table
#'
#' @param pars Parameter list.
#'
#' @return A data.frame joining values to metadata.
#' @export
chlaa_parameters_print <- function(pars) {
  .check_named_list(pars, "pars")
  info <- chlaa_parameter_info()
  info$value <- vapply(info$name, function(nm) as.numeric(pars[[nm]]), numeric(1))
  info
}

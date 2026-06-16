# Economics defaults and decision-analysis helpers

.chlaa_extdata_file <- function(path) {
  pkg_path <- system.file("extdata", path, package = "chlaa")
  if (nzchar(pkg_path) && file.exists(pkg_path)) return(pkg_path)

  local_path <- file.path("inst", "extdata", path)
  if (file.exists(local_path)) return(local_path)

  ""
}

#' Load default economics parameters
#'
#' Reads default unit costs and DALY parameters from `inst/extdata/econ/`.
#' These defaults are intended for demonstration and rapid scenario analysis;
#' they should be replaced with context-specific values for real decisions.
#'
#' @param overrides Optional named list of values overriding defaults.
#'
#' @return A named list of economics parameters. Includes a `sources`
#'   attribute containing source metadata (see `chlaa_econ_sources()`).
#' @export
chlaa_econ_defaults <- function(overrides = NULL) {
  costs_file <- .chlaa_extdata_file(file.path("econ", "unit_costs.csv"))
  daly_file <- .chlaa_extdata_file(file.path("econ", "daly_params.csv"))
  if (!nzchar(costs_file) || !nzchar(daly_file)) {
    stop("Could not find default economics files under inst/extdata/econ", call. = FALSE)
  }

  costs <- utils::read.csv(costs_file, stringsAsFactors = FALSE)
  daly <- utils::read.csv(daly_file, stringsAsFactors = FALSE)

  if (!all(c("name", "value") %in% names(costs))) stop("unit_costs.csv must contain name,value", call. = FALSE)
  if (!all(c("name", "value") %in% names(daly))) stop("daly_params.csv must contain name,value", call. = FALSE)

  vals <- c(stats::setNames(as.list(as.numeric(costs$value)), costs$name),
            stats::setNames(as.list(as.numeric(daly$value)), daly$name))

  if (!is.null(overrides)) {
    .check_named_list(overrides, "overrides")
    vals[names(overrides)] <- overrides
  }

  src <- try(chlaa_econ_sources(), silent = TRUE)
  if (!inherits(src, "try-error")) {
    attr(vals, "sources") <- src
  }

  vals
}

#' Load Economics Assumption Source Metadata
#'
#' Reads source notes and reference links for each economics default.
#'
#' @return A data.frame with one row per economics parameter and columns
#'   including `name`, `source_type`, `citation`, and `source_url`.
#' @export
chlaa_econ_sources <- function() {
  ref_file <- .chlaa_extdata_file(file.path("econ", "references.csv"))
  if (!nzchar(ref_file)) {
    stop("Could not find economics references file under inst/extdata/econ", call. = FALSE)
  }
  utils::read.csv(ref_file, stringsAsFactors = FALSE)
}

#' Cost-effectiveness acceptability table from scenario simulations
#'
#' @param scenario_runs Output from `chlaa_run_scenarios()`.
#' @param baseline Baseline scenario name.
#' @param wtp Numeric vector of willingness-to-pay per DALY averted.
#' @param econ Optional economics override list.
#'
#' @return A data.frame with CEAC probabilities by scenario and WTP.
#' @export
chlaa_ceac <- function(scenario_runs,
                         baseline = "baseline",
                         wtp = seq(0, 5000, by = 250),
                         econ = NULL) {
  .require_suggested("dplyr")

  pair <- .chlaa_particle_econ_delta(
    scenario_runs = scenario_runs,
    baseline = baseline,
    econ = econ
  )

  out <- vector("list", length(wtp))
  for (i in seq_along(wtp)) {
    w <- wtp[[i]]
    tmp <- pair |>
      dplyr::mutate(nmb = w * .data$dalys_averted - .data$cost_diff) |>
      dplyr::group_by(.data$particle) |>
      dplyr::filter(.data$nmb == max(.data$nmb, na.rm = TRUE)) |>
      dplyr::ungroup() |>
      dplyr::group_by(.data$scenario) |>
      dplyr::summarise(prob_best = dplyr::n() / dplyr::n_distinct(pair$particle), .groups = "drop")
    tmp$wtp <- w
    out[[i]] <- tmp
  }

  dplyr::bind_rows(out) |>
    dplyr::select("scenario", "wtp", "prob_best") |>
    dplyr::arrange(.data$wtp, .data$scenario)
}

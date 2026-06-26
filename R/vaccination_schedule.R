#' Generate vaccination schedule with empirical delivery profile
#'
#' Creates a daily vaccination schedule that distributes total doses across
#' a campaign window using an empirical delivery profile. The profile is based
#' on observed OCV campaign patterns showing non-uniform daily delivery rates.
#'
#' @param total_doses Numeric. Total number of vaccine doses to be delivered
#' @param start_date Date. Start date of vaccination campaign
#' @param end_date Date. End date of vaccination campaign (exclusive - last day is end_date - 1)
#' @param outbreak_start Date. Outbreak start date for converting to model time
#' @param profile Numeric vector. Empirical weights for 6-day profile.
#'   Default: c(0.305, 0.377, 0.227, 0.074, 0.014, 0.003) representing
#'   30.5%, 37.7%, 22.7%, 7.4%, 1.4%, 0.3% across days 1-6.
#'
#' @return A data.frame with columns:
#'   \item{date}{Date of vaccination}
#'   \item{time}{Days since outbreak start (for model input)}
#'   \item{doses}{Number of doses to administer on this day}
#'   \item{day_of_campaign}{Day number within campaign (1, 2, 3, ...)}
#'
#' @details
#' The function handles campaigns of different lengths:
#' \itemize{
#'   \item \strong{6-day campaigns}: Uses empirical profile exactly
#'   \item \strong{Shorter campaigns}: Uses first n days of profile, renormalized
#'   \item \strong{Longer campaigns}: Interpolates profile over campaign duration
#' }
#'
#' The empirical profile (30.5%, 37.7%, 22.7%, 7.4%, 1.4%, 0.3%) is based on
#' observed OCV campaign delivery patterns, showing front-loading of doses with
#' peak delivery on day 2.
#'
#' @examples
#' \dontrun{
#' # Generate schedule for a 6-day campaign
#' schedule <- generate_vax_schedule(
#'   total_doses = 343619,
#'   start_date = as.Date("2025-07-14"),
#'   end_date = as.Date("2025-07-20"),
#'   outbreak_start = as.Date("2025-02-24")
#' )
#'
#' # Verify total doses
#' sum(schedule$doses)  # Should equal 343619
#'
#' # For a longer campaign (23 days)
#' schedule_long <- generate_vax_schedule(
#'   total_doses = 423814,
#'   start_date = as.Date("2025-04-07"),
#'   end_date = as.Date("2025-04-30"),
#'   outbreak_start = as.Date("2025-01-06")
#' )
#' }
#'
#' @export
generate_vax_schedule <- function(
    total_doses,
    start_date,
    end_date,
    outbreak_start,
    profile = c(0.305, 0.377, 0.227, 0.074, 0.014, 0.003)
) {

  # Validate inputs
  if (!inherits(start_date, "Date") || !inherits(end_date, "Date")) {
    stop("start_date and end_date must be Date objects")
  }
  if (!inherits(outbreak_start, "Date")) {
    stop("outbreak_start must be a Date object")
  }
  if (start_date >= end_date) {
    stop("start_date must be before end_date")
  }
  if (total_doses <= 0) {
    stop("total_doses must be positive")
  }
  if (length(profile) != 6 || abs(sum(profile) - 1.0) > 1e-6) {
    stop("profile must be a vector of 6 weights summing to 1")
  }

  # Generate campaign dates (end_date is exclusive)
  campaign_dates <- seq(start_date, end_date - 1, by = "day")
  n_days <- length(campaign_dates)

  if (n_days == 0) {
    stop("Campaign duration is zero days")
  }

  # Calculate daily doses based on campaign length
  if (n_days == 6) {
    # Exact match to empirical profile
    daily_doses <- total_doses * profile

  } else if (n_days < 6) {
    # Shorter campaign: use first n days of profile, renormalize
    daily_weights <- profile[1:n_days] / sum(profile[1:n_days])
    daily_doses <- total_doses * daily_weights

  } else {
    # Longer campaign: smooth profile over duration
    # Map each campaign day to position in 6-day template
    template_positions <- seq(0, 1, length.out = n_days)
    template_indices <- template_positions * 5 + 1  # 1 to 6

    # Interpolate weights from template using linear interpolation
    daily_weights <- stats::approx(
      x = 1:6,
      y = profile,
      xout = template_indices,
      method = "linear",
      rule = 2
    )$y

    # Renormalize to ensure they sum to 1
    daily_weights <- daily_weights / sum(daily_weights)
    daily_doses <- total_doses * daily_weights
  }

  # Calculate model time (days since outbreak start)
  model_time <- as.numeric(campaign_dates - outbreak_start)

  # Create output data frame
  schedule <- data.frame(
    date = campaign_dates,
    time = model_time,
    doses = daily_doses,
    day_of_campaign = 1:n_days,
    stringsAsFactors = FALSE
  )

  return(schedule)
}


#' Prepare vaccination schedule arrays for odin model
#'
#' Converts a vaccination schedule data frame into the format needed for
#' passing to the odin cholera model as interpolation arrays.
#'
#' @param schedule Data frame from \code{generate_vax_schedule()}
#'
#' @return A list with elements:
#'   \item{time}{Numeric vector of time points}
#'   \item{doses}{Numeric vector of doses at each time point}
#'
#' @details
#' This function extracts the time and doses columns from a schedule and
#' returns them as a named list suitable for passing to the odin model's
#' interpolation parameters (vax1_schedule_time, vax1_schedule_doses, etc.)
#'
#' @examples
#' \dontrun{
#' schedule <- generate_vax_schedule(
#'   total_doses = 343619,
#'   start_date = as.Date("2025-07-14"),
#'   end_date = as.Date("2025-07-20"),
#'   outbreak_start = as.Date("2025-02-24")
#' )
#'
#' arrays <- prepare_vax_arrays(schedule)
#' # Use in model:
#' # model$new(vax1_schedule_time = arrays$time,
#' #           vax1_schedule_doses = arrays$doses, ...)
#' }
#'
#' @export
prepare_vax_arrays <- function(schedule) {
  if (!is.data.frame(schedule)) {
    stop("schedule must be a data frame")
  }
  if (!all(c("time", "doses") %in% names(schedule))) {
    stop("schedule must contain 'time' and 'doses' columns")
  }

  list(
    time = schedule$time,
    doses = schedule$doses
  )
}


#' Generate empty vaccination schedule for zones without vaccination
#'
#' Creates a dummy vaccination schedule with zero doses for health zones
#' that did not receive vaccination. This ensures the model receives valid
#' array inputs even when no vaccination occurred.
#'
#' @param outbreak_start Date. Outbreak start date
#'
#' @return A data.frame with a single row of zeros
#'
#' @details
#' The odin model requires vaccination schedule arrays even for zones without
#' vaccination. This function creates a minimal valid schedule with zero doses
#' at time 0.
#'
#' @examples
#' \dontrun{
#' # For a zone without vaccination
#' empty_schedule <- generate_empty_vax_schedule(
#'   outbreak_start = as.Date("2025-01-06")
#' )
#' }
#'
#' @export
generate_empty_vax_schedule <- function(outbreak_start) {
  if (!inherits(outbreak_start, "Date")) {
    stop("outbreak_start must be a Date object")
  }

  data.frame(
    date = outbreak_start,
    time = 0,
    doses = 0,
    day_of_campaign = 0,
    stringsAsFactors = FALSE
  )
}

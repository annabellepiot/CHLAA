# =====================================================================
# 00_setup.R  —  Shared configuration and helpers for the CHLAA audit
# =====================================================================
# Source this at the top of every V*.R script:  source("00_setup.R")
# Edit ONLY the PATHS block below to match your machine / HPC.
# The scripts target your REAL fits and private data as requested.
# ---------------------------------------------------------------------

suppressMessages({
  library(chlaa)
  library(dplyr)
})

## ---- PATHS (edit these) --------------------------------------------
# Point these at the same locations used by analysis/01_02_fitting_all_HZs.R
PROJ     <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA"
DATA_DIR <- file.path(PROJ, "analysis", "data")
FIG_DIR  <- file.path(PROJ, "figures")
RDS_DIR  <- file.path(FIG_DIR, ".rds files")
TAB_DIR  <- file.path(FIG_DIR, "tables")
OUT_DIR  <- file.path("output")            # where these scripts write
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- Global constants (must match the fitting script) --------------
H_REF   <- 1.0
POP_REF <- 516000
RR_FIXED <- 0.30
K_R0    <- POP_REF * 5.1446e-3 / H_REF     # the constant under audit (C4)

## ---- Discover fitted health zones ----------------------------------
list_hz_fits <- function() {
  f <- list.files(RDS_DIR, pattern = "_fit\\.rds$", full.names = FALSE)
  f <- setdiff(f, list.files(RDS_DIR, pattern = "_FAILED\\.rds$"))
  sort(gsub("_fit\\.rds$", "", f))
}

load_fit_obj <- function(hz) {
  p <- file.path(RDS_DIR, sprintf("%s_fit.rds", hz))
  if (!file.exists(p)) stop("Fit not found: ", p)
  readRDS(p)
}

## ---- Posterior helpers ---------------------------------------------
# Median of the SAMPLED (transformed) coordinates, in packer order.
theta_median_sampled <- function(fit, burnin = 0.25) {
  tr <- chlaa_fit_trace(fit, burnin = burnin, scale = "sampled")
  pk <- attr(fit, "packer", exact = TRUE)
  nm <- pk$names()
  vapply(nm, function(p) median(tr$value[tr$parameter == p], na.rm = TRUE), numeric(1))
}

# Turn a sampled-coordinate vector into a full model parameter list
# (runs the packer's process function, so transforms + fixed params apply).
theta_to_pars <- function(fit, theta) {
  pk <- attr(fit, "packer", exact = TRUE)
  pk$unpack(as.numeric(theta))
}

## ---- Rebuild the exact particle filter used at fit time ------------
# Uses the stored data / time_start / packer so the likelihood matches.
rebuild_filter <- function(fit, n_particles, seed = 1L, dt = NULL) {
  gen  <- chlaa:::chlaa_generator()
  dat  <- attr(fit, "data",       exact = TRUE)
  ts   <- attr(fit, "time_start", exact = TRUE)
  args <- list(generator = gen, time_start = ts, data = dat,
               n_particles = n_particles, seed = seed)
  if (!is.null(dt)) args$dt <- dt          # only if your dust2 supports it
  do.call(dust2::dust_filter_create, args)
}

# Evaluate the marginal log-likelihood at a sampled-coord theta.
loglik_at <- function(fit, theta, n_particles, seed = 1L, dt = NULL) {
  filt <- rebuild_filter(fit, n_particles, seed = seed, dt = dt)
  pars <- theta_to_pars(fit, theta)
  # function name differs across dust2 versions; try the common ones
  if (exists("dust_likelihood_run", where = asNamespace("dust2")))
    dust2::dust_likelihood_run(filt, pars)
  else
    dust2::dust_filter_run(filt, pars)
}

verdict <- function(pass, msg) cat(if (pass) "PASS   " else "CONCERN", "| ", msg, "\n")
cat("Setup loaded. HZs available:", paste(list_hz_fits(), collapse = ", "), "\n")

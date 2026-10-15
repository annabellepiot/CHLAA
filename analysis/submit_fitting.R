#=========================================================================
#Submit per-HZ PMCMC fits to the DIDE cluster via hipercow
#=========================================================================
#
#Replaces run_All_HZ_fitting.pbs: one task per health zone, each running
#fit_hz() from 01_02_fitting_all_HZs.R. Run from the repo root:
#
#  source("analysis/submit_fitting.R")
#
#Set `hzs` first to submit a subset, e.g. hzs <- "kirotshe".
#
#The cluster installs chlaa from GitHub (pkgdepends.txt), so after pushing
#changes to R/ or src/ re-run hipercow::hipercow_provision() before
#submitting.
#
#=========================================================================

library(hipercow)

if (!exists("hzs")) {
    hzs <- c(
        "alunguli", "bumbu", "goma", "kindu", "kingabwa", "kirotshe",
        "kokolo", "limete", "maluku_i", "ngiri_ngiri", "nyiragongo", "opala"
    )
}

hipercow_environment_create("fitting",
    sources = c("analysis/hipercow_setup.R", "analysis/01_02_fitting_all_HZs.R")
)

resources <- hipercow_resources(cores = 4, memory_per_node = "16G", max_runtime = "24h")

fit_bundle <- task_create_bulk_expr(fit_hz(hz, verbose = TRUE),
    data.frame(hz = hzs),
    environment = "fitting", resources = resources, bundle_name = "chlaa_fit"
)

cat("\nSubmitted", length(hzs), "fit task(s): bundle", fit_bundle$name, "\n")
cat("Monitor with hipercow_bundle_status(fit_bundle) and task_log_show(<task id>)\n")

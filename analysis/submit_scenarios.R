#=========================================================================
#Submit per-HZ scenario analysis to the DIDE cluster via hipercow
#=========================================================================
#
#Replaces run_scenario_analysis.pbs: one task per health zone with a
#completed fit, each running run_scenarios_hz() from
#02_02_scenario_analysis_all_HZs.R. Run from the repo root:
#
#  source("analysis/submit_scenarios.R")
#
#Set `hzs` first to submit a subset. Once every task has finished, build
#the cross-HZ figures locally from the cached .rds files with
#
#  options(chlaa.skip_main = TRUE)
#  source("analysis/02_02_scenario_analysis_all_HZs.R")
#
#=========================================================================

library(hipercow)

rds_dir <- file.path("figures", ".rds files")

if (!exists("hzs")) {
    all_fits <- list.files(rds_dir, pattern = "_fit\\.rds$")
    hzs <- sort(gsub("_fit\\.rds$", "", all_fits[!grepl("comparative", all_fits)]))
}

missing_fits <- hzs[!file.exists(file.path(rds_dir, sprintf("%s_fit.rds", hzs)))]
if (length(missing_fits) > 0) {
    stop("Fit artifact not found for: ", paste(missing_fits, collapse = ", "),
         " - run analysis/submit_fitting.R for these HZs first.")
}

hipercow_environment_create("scenarios",
    sources = c("analysis/hipercow_setup.R", "analysis/02_02_scenario_analysis_all_HZs.R")
)

resources <- hipercow_resources(cores = 4, memory_per_node = "16G", max_runtime = "4h")

scenario_bundle <- task_create_bulk_expr(run_scenarios_hz(hz, verbose = TRUE),
    data.frame(hz = hzs),
    environment = "scenarios", resources = resources, bundle_name = "chlaa_scenarios"
)

cat("\nSubmitted", length(hzs), "scenario task(s): bundle", scenario_bundle$name, "\n")
cat("Monitor with hipercow_bundle_status(scenario_bundle) and task_log_show(<task id>)\n")

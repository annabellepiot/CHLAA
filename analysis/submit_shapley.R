#=========================================================================
#Submit per-HZ Shapley decomposition to the DIDE cluster via hipercow
#=========================================================================
#
#Replaces run_shapley_analysis.pbs: one task per health zone with a
#completed fit, each running compute_shapley_hz() from
#05_Shapley_analysis.R (N_THREADS = 4, matching cores below). Run from
#the repo root:
#
#  source("analysis/submit_shapley.R")
#
#Set `hzs` first to submit a subset. Once every task has finished, build
#the cross-HZ figures locally from the cached .rds files with
#
#  options(chlaa.skip_main = TRUE)
#  source("analysis/05_Shapley_analysis.R")
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

hipercow_environment_create("shapley",
    sources = c("analysis/hipercow_setup.R", "analysis/05_Shapley_analysis.R")
)

resources <- hipercow_resources(cores = 4, memory_per_node = "16G", max_runtime = "30")

shapley_bundle <- task_create_bulk_expr(compute_shapley_hz(hz, verbose = TRUE),
    data.frame(hz = hzs),
    environment = "shapley", resources = resources, bundle_name = "chlaa_shapley"
)

cat("\nSubmitted", length(hzs), "Shapley task(s): bundle", shapley_bundle$name, "\n")
cat("Monitor with hipercow_bundle_status(shapley_bundle) and task_log_show(<task id>)\n")

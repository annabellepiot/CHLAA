# Multi-Site PMCMC Fitting

This folder contains the all-health-zone fitting workflow for the CHLAA model.
It follows the same pMCMC conventions as `01_fitting.R`:

- weekly data are fit with `obs_interval = 7` and `time_start = 0`
- the fitted model uses transformed parameters
- `seek_severe` stays fixed
- `E0` is seeded from outbreak size and pre-outbreak information where possible
- the posterior predictive plot uses the same `plot_case_fit` styling as the master script

## Main files

- `01_v2_fitting_all_HZs.R`: core fitter for one health zone at a time
- `run_All_HZ_fitting.pbs`: PBS array wrapper for batch fitting
- `check_fitting_status.R`: checks which zones succeeded or failed
- `collect_fitted_parameters.R`: gathers posterior summaries across completed fits

## Run

Submit the array job:

```bash
qsub run_All_HZ_fitting.pbs
```

Check progress:

```bash
Rscript check_fitting_status.R
```

Collect summaries and figures:

```bash
Rscript collect_fitted_parameters.R
```

## Statistical notes

- Each health zone is fit independently; there is no hierarchical shrinkage across zones.
- Small outbreaks may be weakly identified, especially for `reporting_rate` and `obs_size`.
- If a zone is sparse or noisy, tighter priors or fewer free parameters may be needed.
- Cross-zone comparison is only valid when each fit is interpreted as a separate outbreak.

## Outputs

- `<hz>_fit.rds`: fit object and posterior summary
- `<hz>_FAILED.rds`: error record if a fit fails
- `diagnosis_<hz>_production_fit.png`: posterior predictive fit plot

## Data source

Intervention dates are read from `analysis/data/hz_parameters.csv`.
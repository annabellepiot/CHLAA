# Cholera Outbreak Simulation and Anticipatory Action Modelling

<!-- Logo banner - add your collaborator logos here -->
<p align="center">
  <img src="man/figures/logo1.png" height="80" />
  &nbsp;&nbsp;&nbsp;
  <img src="man/figures/logo2.png" height="80" />
</p>

---

## Overview

Stochastic cholera outbreak simulation and anticipatory action (AA) modelling tools.

The compartmental model behind this chlaa package is inspired by the published results from “A simulation-based policy analysis of anticipatory action for cholera outbreaks, Democratic Republic of the Congo” by Loo, P. S., Rajah, J. K., de Leon, H. J. H., Kopainsky, B., & Milano, L. (2025). https://doi.org/10.2471/BLT.25.293226. 

### Fitting the model

The model is implemented as a stochastic compartmental model using the odin2/dust2/monty stack. This provides a framework that can be applied to any healthzone in the DRC (varying population sizes, endemic and epidemic), and captures inherent randomness - an important methodological advancement over previous deterministic work. 

The main fitting script is `analysis/01_02_fitting_all_HZs.R`, which can be submitted to an HPC scheduler as an array job via `analysis/run_All_HZ_fitting.pbs`. Other scripts are available for single HZs (Kirotshe, Goma) to understand the workflow, and can be run directly in R. 


### Analysing the effect of anticipatory action
Eight intervention types are modelled: chlorination, hygiene kits, latrines, cholera treatment centres (CTC), oral rehydration corners (ORC), oral cholera vaccination (OCV), surveillance, and case-area targeted interventions (CATI).

The main scenario analysis script is 'analysis/02_02_scenario_analysis_all_HZs.R'. This includes comparison of : baseline (no AA activations), real world AA activation, AA plus vaccination (1 and 2 dose). 

*Note*: The Modelling Assumptions script is not up to date, disregard. 

## Repository Layout

```r
chlaa/
├── R/                 # Reusable R functions (model, fitting, scenarios)
├── analysis/          # Analysis scripts and reports (.R/.qmd)
├── analysis/data/     # Input data (cleaned outbreak data and parameters) ###This will be made available in the future upon partner approvals. 
├── inst/extdata/      # Minimal package example data ##Sufficient to run vignette examples
├── figures/           # Generated figures
├── output/            # Model output and results
├── inst/              # odin/dust model definitions
├── man/               # Function documentation
├── tests/             # Unit tests
├── vignettes_src/     # Source for extended pre-rendered vignettes
└── DESCRIPTION        # Package metadata and dependencies
```

## Installation

```r
# Install dependencies
remotes::install_github("annabellepiot/chlaa")

```

### Rendering

Some analysis reports are written as Quarto (`.qmd`) documents in the `analysis/` directory. To render them locally:

1. Install [Quarto](https://quarto.org/docs/get-started/)
2. Install the R package and its dependencies (see Installation above)
3. Render the reports:

```bash
quarto render "analysis/Cases and Deaths.qmd"
quarto render "analysis/Model Assumptions.qmd"
```

## Data Sources
The data is proprietaty of United Nations Office for the Coordination of Humanitarian Affairs (OCHA), Centre for Humanitarian Data.

## Dependencies
Key dependencies: odin2, dust2, monty (from mrc-ide).

This package vendors the generated dust2 C++ model code, so you can run simulation and fitting without having **odin2** installed at runtime. (odin2 is only needed to regenerate the bundled code; see `docs/packaging.md`.)

## License


## Citation
How to cite this work.

## Minimal simulation

```r
library(chlaa)

pars <- chlaa_parameters(Sev0 = 2)
sim <- chlaa_simulate(pars, time = 0:180, n_particles = 50, dt = 0.25, seed = 1)
sim 
```

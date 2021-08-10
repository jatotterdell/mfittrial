# mfittrial

R functions used for MFIT trial simulations.

Currently four versions:

- `simulate_trial.R`: original functions used for grant submission simulations March 2020
- `simulate_trial_model2.R`: updated functions
- `simulate_trial_longitudinal.R`: functions for simulating longitudinal design
- `simulate_trial_with_control.R`: (current version), simulate trial treating one arm as fixed control

The script `simulate_trial_with_control.R` is the current version 
for the latest protocol updates following discussions at the stats meeting on 2021-08-02.

## Installation

```
remotes::install_github("jatotterdell/mfittrial")
```
# Run from the extracted package directory after installing version 0.5.3.
# The four scripts can also be sourced separately, pausing between stages.

source("analysis/01_calibrate_prior.R")

source("analysis/02_historical_meta.R")

source("analysis/03_hope_update.R")

source("analysis/04_programs.R")

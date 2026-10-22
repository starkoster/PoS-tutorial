# Optional manual review. This script is deliberately outside run_all.R.
# Stan's own sampling warnings remain visible in the main workflow.

output_dir <- "results_R_readable"

meta <- readRDS(file.path(output_dir, "stan_meta", "fit.rds"))
hope <- readRDS(file.path(output_dir, "stan_hope", "fit.rds"))

# Inspect divergences, treedepth and energy diagnostics by chain.
meta$diagnostic_summary()
hope$diagnostic_summary()

# Inspect mixing and effective sample sizes for the main model parameters.
print(meta$summary(variables = c("mu_eta", "tau_eta")))
print(hope$summary(variables = c("eta", "log_ed50")))

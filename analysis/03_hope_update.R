library(MigraineEmaxSim)

# Check the version actually loaded in this R session, before fitting Stan.
if (utils::compareVersion(as.character(getNamespaceVersion("MigraineEmaxSim")),
                          "0.5.3") < 0) {
  stop("Install MigraineEmaxSim >= 0.5.3 and restart R before running this script.",
       call. = FALSE)
}

output_dir <- "results_R_readable"
predictive <- readRDS(file.path(output_dir, "02_predictive_prior.rds"))
arms <- read_analysis_arms()

# STEP 1: enter a separate, coherent lognormal ED50 elicitation.
# These values give median 20 mg and central 90% interval 6.32-63.25 mg.
ed50_prior <- make_ed50_prior(
  q05 = sqrt(40),
  q50 = 20,
  q95 = sqrt(4000)
)

# Retain the chosen elicitation so output-only recovery uses the same ED50 prior.
saveRDS(ed50_prior, file.path(output_dir, "03_ed50_prior.rds"))

# STEP 2: HOPE updates the joint historical prior and the independent ED50 prior.
hope <- fit_emax(
  arms = arms,
  predictive_prior = predictive,
  ed50_prior = ed50_prior,
  output_dir = file.path(output_dir, "stan_hope")
)

# STEP 3: save WHOLE posterior rows for the trial simulations.
# Use the installed extractor explicitly, avoiding older global patch functions.
draws <- MigraineEmaxSim::emax_draws(hope)

# Each clinical quantity must be an ordinary numeric vector before file export.
draws <- MigraineEmaxSim::as_numeric_draws(draws)
saveRDS(draws, file.path(output_dir, "03_emax_draws.rds"))
readr::write_csv(draws, file.path(output_dir, "emax_posterior.csv"))

print(clinical_quantiles(draws))
readr::write_csv(clinical_quantiles(draws), file.path(output_dir, "emax_summary.csv"))

# Draw the actual Gaussian transfer prior used by the HOPE model for comparison.
set.seed(4400)
prior_eta <- matrix(stats::rnorm(50000 * 6), ncol = 6) %*% chol(predictive$cov)
prior_eta <- sweep(prior_eta, 2, predictive$mean, "+")

save_analysis_plot(
  MigraineEmaxSim::plot_distributions(list(
    Prior = MigraineEmaxSim::clinical_from_eta(prior_eta),
    HOPE = draws
  )),
  "03_hope_update",
  file.path(output_dir, "figures")
)

save_analysis_plot(
  plot_emax_curves(draws),
  "03_placebo_adjusted_emax",
  file.path(output_dir, "figures")
)

save_analysis_plot(
  plot_emax_curves(draws, adjusted = FALSE),
  "03_total_change_emax",
  file.path(output_dir, "figures")
)

save_analysis_plot(
  plot_ed50_update(draws, ed50_prior),
  "03_ed50_update",
  file.path(output_dir, "figures"),
  width = 6,
  height = 4
)

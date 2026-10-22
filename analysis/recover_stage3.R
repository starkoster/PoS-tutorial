# Recover the Stage 3 CSVs and figures without rerunning either Stan model.
# Install version 0.5.3 and restart R first; do not source the old patch scripts.
library(MigraineEmaxSim)

if (utils::compareVersion(as.character(getNamespaceVersion("MigraineEmaxSim")),
                          "0.5.3") < 0) {
  stop("Install MigraineEmaxSim >= 0.5.3 and restart R before running this script.",
       call. = FALSE)
}

# Point this at your EXISTING results directory if it is in a different folder.
output_dir <- "results_R_readable"

# Prefer the saved draw table: Stage 3 saved it before the CSV error occurred.
# If that table is absent, extract from the saved fit instead. Neither path fits.
draw_file <- file.path(output_dir, "03_emax_draws.rds")
fit_file <- file.path(output_dir, "stan_hope", "fit.rds")

if (file.exists(draw_file)) {

  draws <- readRDS(draw_file)

} else if (file.exists(fit_file)) {

  hope <- readRDS(fit_file)
  draws <- MigraineEmaxSim::emax_draws(hope)

} else {

  stop("No saved Emax draws or fit found. Set output_dir to your existing results folder.",
       call. = FALSE)
}

# Remove legacy matrix columns while preserving every value and posterior row.
draws <- MigraineEmaxSim::as_numeric_draws(draws)

# Historical prediction supplies the same transferred prior used in Stage 3.
predictive <- readRDS(file.path(output_dir, "02_predictive_prior.rds"))

# New runs save the ED50 elicitation. For earlier runs, the fallback below is
# the original default: edit its quantiles if you used a different elicitation.
ed50_file <- file.path(output_dir, "03_ed50_prior.rds")

if (file.exists(ed50_file)) {

  ed50_prior <- readRDS(ed50_file)

} else {

  ed50_prior <- MigraineEmaxSim::make_ed50_prior(
    q05 = sqrt(40),
    q50 = 20,
    q95 = sqrt(4000)
  )
}

# Repeat only the export and figure section of Stage 3, using the existing draws.
saveRDS(draws, file.path(output_dir, "03_emax_draws.rds"))
readr::write_csv(draws, file.path(output_dir, "emax_posterior.csv"))

print(MigraineEmaxSim::clinical_quantiles(draws))
readr::write_csv(MigraineEmaxSim::clinical_quantiles(draws), file.path(output_dir, "emax_summary.csv"))

# Draw the actual Gaussian transfer prior used by the HOPE model for comparison.
set.seed(4400)
prior_eta <- matrix(stats::rnorm(50000 * 6), ncol = 6) %*% chol(predictive$cov)
prior_eta <- sweep(prior_eta, 2, predictive$mean, "+")

MigraineEmaxSim::save_analysis_plot(
  MigraineEmaxSim::plot_distributions(list(
    Prior = MigraineEmaxSim::clinical_from_eta(prior_eta),
    HOPE = draws
  )),
  "03_hope_update",
  file.path(output_dir, "figures")
)

MigraineEmaxSim::save_analysis_plot(
  MigraineEmaxSim::plot_emax_curves(draws),
  "03_placebo_adjusted_emax",
  file.path(output_dir, "figures")
)

MigraineEmaxSim::save_analysis_plot(
  MigraineEmaxSim::plot_emax_curves(draws, adjusted = FALSE),
  "03_total_change_emax",
  file.path(output_dir, "figures")
)

MigraineEmaxSim::save_analysis_plot(
  MigraineEmaxSim::plot_ed50_update(draws, ed50_prior),
  "03_ed50_update",
  file.path(output_dir, "figures"),
  width = 6,
  height = 4
)

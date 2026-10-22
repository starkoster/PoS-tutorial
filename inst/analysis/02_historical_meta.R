library(MigraineEmaxSim)

# Check the version actually loaded in this R session, before doing any work.
if (utils::compareVersion(as.character(getNamespaceVersion("MigraineEmaxSim")),
                          "0.5.3") < 0) {
  stop("Install MigraineEmaxSim >= 0.5.3 and restart R before running this script.",
       call. = FALSE)
}

output_dir <- "results_R_readable"
calibration <- readRDS(file.path(output_dir, "01_calibration.rds"))

# STEP 1: read arm summaries and apply the documented HOPE input corrections.
# HOPE is then excluded from this historical fit and reserved for Stage 3.
arms <- read_analysis_arms()

readr::write_csv(arms, file.path(output_dir, "analysis_arms.csv"))

# STEP 2: fit the hierarchical historical model in Stan.
# Only active change likelihoods are discounted; baseline/placebo retain full weight.
meta <- fit_meta(
  arms = arms,
  prior = calibration$prior,
  output_dir = file.path(output_dir, "stan_meta"),
  active_weight = 0.5
)

# STEP 3: use the posterior prediction for a NEW study as the next model's prior.
# Its joint covariance retains uncertainty and between-study heterogeneity.
predictive <- predictive_prior_from_meta(meta)

saveRDS(predictive, file.path(output_dir, "02_predictive_prior.rds"))

meta_summary <- clinical_quantiles(predictive$clinical_draws)
print(meta_summary)
readr::write_csv(meta_summary, file.path(output_dir, "historical_prediction.csv"))

save_analysis_plot(
  plot_distributions(list(
    Prior = calibration$draws,
    Historical = predictive$clinical_draws
  )),
  "02_meta_update",
  file.path(output_dir, "figures")
)

save_analysis_plot(
  plot_meta_intervals(calibration$draws, meta, dplyr::filter(arms, title != "HOPE")),
  "02_historical_trials",
  file.path(output_dir, "figures"),
  height = 7
)

save_analysis_plot(
  plot_em_cm(list(Prior = calibration$draws, Historical = predictive$clinical_draws)),
  "02_meta_em_cm",
  file.path(output_dir, "figures")
)

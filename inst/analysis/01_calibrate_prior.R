library(MigraineEmaxSim)

# Check the version actually loaded in this R session, before doing any work.
if (utils::compareVersion(as.character(getNamespaceVersion("MigraineEmaxSim")),
                          "0.5.3") < 0) {
  stop("Install MigraineEmaxSim >= 0.5.3 and restart R before running this script.",
       call. = FALSE)
}

output_dir <- "results_R_readable"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# STEP 1: enter the elicited clinical quantiles, using negative effects for benefit.
# Edit this table directly, or replace it with readr::read_csv("my_targets.csv").
targets <- default_targets()

print(targets)

# STEP 2: calibrate the latent hyperpriors to these new-study clinical targets.
# The heterogeneity fractions are separate assumptions, not additional targets.
calibration <- calibrate_priors(
  targets = targets,
  heterogeneity_fraction = c(4 / 13, 4 / 13, rep(0.36, 4)),
  lkj_shape = 2
)

# STEP 3: inspect achieved quantiles and save the prior for the historical model.
print(calibration$comparison)

saveRDS(calibration, file.path(output_dir, "01_calibration.rds"))
readr::write_csv(calibration$comparison, file.path(output_dir, "prior_calibration.csv"))

save_analysis_plot(
  plot_calibration(calibration),
  "01_prior_calibration",
  file.path(output_dir, "figures")
)

save_analysis_plot(
  plot_em_cm(list(Prior = calibration$draws)),
  "01_prior_em_cm",
  file.path(output_dir, "figures")
)

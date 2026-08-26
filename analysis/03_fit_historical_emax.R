source(file.path(if (dir.exists("analysis")) "analysis" else ".", "_setup.R"))
dat <- load_trial_csv(file.path(pkg_root, "data", "hist_trials.csv"))
split <- prepare_trial_data(dat, hope_title = "HOPE")
fit_hist <- readRDS(file.path(out_dir, "fit_historical_meta.rds"))
eta_prior <- extract_eta_prior(fit_hist, variance_inflation = 1)
ed50_prior <- make_ed50_prior(median = 100, lower = 20, upper = 300, mass = 0.90)
fit_emax <- fit_historical_emax(
  hope = split$hope,
  eta_prior = eta_prior,
  ed50_prior = ed50_prior,
  stan_file = file.path(pkg_root, "stan", "historical_emax.stan"),
  iter_warmup = 1000,
  iter_sampling = 3000,
  seed = 456
)
saveRDS(fit_emax, file.path(out_dir, "fit_historical_emax.rds"))
print(fit_emax$summary(c("e0_em", "e0_cm", "emax_em", "emax_cm", "ed50")))

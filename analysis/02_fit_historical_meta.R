source(file.path(if (dir.exists("analysis")) "analysis" else ".", "_setup.R"))
dat <- load_trial_csv(file.path(pkg_root, "data", "hist_trials.csv"))
split <- prepare_trial_data(dat, hope_title = "HOPE")
priors <- default_eta_priors()
fit_hist <- fit_historical_meta(
  historical = split$historical,
  priors = priors,
  use_historical_active = TRUE,
  active_borrow_weight = 0.25,
  stan_file = file.path(pkg_root, "stan", "historical_meta.stan"),
  iter_warmup = 2000,
  iter_sampling = 4000,
  seed = 123
)
saveRDS(fit_hist, file.path(out_dir, "fit_historical_meta.rds"))
print(fit_hist$summary(c("e0_em_new", "e0_cm_new", "emax_em_new", "emax_cm_new")))

#' Run the full historical-meta, Emax, and program-simulation workflow
#'
#' This convenience wrapper reads the uploaded CSV, fits the constrained
#' historical meta-analysis, fits the HOPE ED50 Emax model, and simulates both
#' development programs.
#'
#' @param root Package source root.
#' @param csv_path Path to the uploaded CSV.
#' @param output_dir Directory for output CSV and RDS files.
#' @param n_program_sims Number of program simulations. Defaults to 5000.
#' @param use_historical_active Include historical active-reference rows?
#' @param active_borrow_weight Likelihood weight for historical active-reference rows.
#' @param ed50_median Prior median ED50 for the new drug.
#' @param ed50_lower Lower plausible ED50.
#' @param ed50_upper Upper plausible ED50.
#' @param hist_iter_warmup Historical model warmup iterations.
#' @param hist_iter_sampling Historical model sampling iterations.
#' @param emax_iter_warmup Emax model warmup iterations.
#' @param emax_iter_sampling Emax model sampling iterations.
#' @return Invisibly, a list of fits, draws, and simulations.
#' @export
run_full_analysis <- function(root = ".",
                              csv_path = file.path(root, "data", "hist_trials.csv"),
                              output_dir = file.path(root, "analysis", "output"),
                              n_program_sims = 5000,
                              use_historical_active = TRUE,
                              active_borrow_weight = 0.25,
                              ed50_median = 100,
                              ed50_lower = 20,
                              ed50_upper = 300,
                              hist_iter_warmup = 2000,
                              hist_iter_sampling = 4000,
                              emax_iter_warmup = 1000,
                              emax_iter_sampling = 3000) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  dat <- load_trial_csv(csv_path)
  split <- prepare_trial_data(dat, hope_title = "HOPE")

  priors <- default_eta_priors()
  prior_spec <- eta_prior_table(priors)
  utils::write.csv(prior_spec, file.path(output_dir, "meta_analysis_prior_specification.csv"), row.names = FALSE)

  prior_pred <- simulate_eta_prior(priors, n = 100000, seed = 11)
  prior_checks <- prior_predictive_summary(prior_pred, threshold = 15)
  utils::write.csv(prior_checks, file.path(output_dir, "prior_predictive_checks.csv"), row.names = FALSE)

  fit_hist <- fit_historical_meta(
    historical = split$historical,
    priors = priors,
    use_historical_active = use_historical_active,
    active_borrow_weight = active_borrow_weight,
    stan_file = file.path(root, "stan", "historical_meta.stan"),
    iter_warmup = hist_iter_warmup,
    iter_sampling = hist_iter_sampling,
    seed = 123
  )
  saveRDS(fit_hist, file.path(output_dir, "fit_historical_meta.rds"))

  hist_summary <- fit_hist$summary(c(
    "baseline_em_new", "baseline_cm_new",
    "e0_em_new", "e0_cm_new",
    "emax_em_new", "emax_cm_new",
    "pbo_red_em_new", "pbo_red_cm_new",
    "drug_red_em_new", "drug_red_cm_new",
    "active_red_em_new", "active_red_cm_new"
  ))
  utils::write.csv(hist_summary, file.path(output_dir, "historical_meta_summary.csv"), row.names = FALSE)

  eta_prior <- extract_eta_prior(fit_hist, variance_inflation = 1)
  ed50_prior <- make_ed50_prior(ed50_median, ed50_lower, ed50_upper, mass = 0.90)
  dose_grid <- seq(0, max(split$hope$dose, 20, 200, 300), length.out = 150)

  fit_emax <- fit_historical_emax(
    hope = split$hope,
    eta_prior = eta_prior,
    ed50_prior = ed50_prior,
    dose_grid = dose_grid,
    stan_file = file.path(root, "stan", "historical_emax.stan"),
    iter_warmup = emax_iter_warmup,
    iter_sampling = emax_iter_sampling,
    seed = 456
  )
  saveRDS(fit_emax, file.path(output_dir, "fit_historical_emax.rds"))

  emax_summary <- fit_emax$summary(c(
    "baseline_em", "baseline_cm", "e0_em", "e0_cm",
    "emax_em", "emax_cm", "ed50", "pbo_red_em", "pbo_red_cm",
    "drug_red_em", "drug_red_cm", "active_red_em", "active_red_cm"
  ))
  utils::write.csv(emax_summary, file.path(output_dir, "historical_emax_summary.csv"), row.names = FALSE)

  emax_draws <- extract_emax_draws(fit_emax)
  utils::write.csv(emax_draws, file.path(output_dir, "emax_posterior_draws.csv"), row.names = FALSE)

  prog1 <- simulate_program1_direct_phase3(
    draws = emax_draws,
    n_sim = n_program_sims,
    seed = 101,
    sigma = 5,
    doses = c(0, 20, 200),
    n_per_arm = 233,
    thresholds = c(2, 3, 4),
    alpha_p3 = 0.05
  )

  prog2 <- simulate_program2_phase2_then_phase3(
    draws = emax_draws,
    n_sim = n_program_sims,
    seed = 202,
    sigma = 5,
    p2_doses = c(0, 20, 100, 300),
    p2_n_per_arm = 140,
    p2_prop_em = 0.50,
    p2_fixed_split = TRUE,
    p2_effect_go = 2,
    p2_alpha_dr = 0.05,
    p3_n_per_arm = 174,
    thresholds = c(2, 3, 4),
    alpha_p3 = 0.05
  )

  saveRDS(prog1, file.path(output_dir, "program1_simulation.rds"))
  saveRDS(prog2, file.path(output_dir, "program2_simulation.rds"))

  utils::write.csv(prog1$by_pop_dose, file.path(output_dir, "program1_by_pop_dose.csv"), row.names = FALSE)
  utils::write.csv(prog1$significance_only, file.path(output_dir, "program1_significance_only.csv"), row.names = FALSE)
  utils::write.csv(prog1$threshold_program_success, file.path(output_dir, "program1_threshold_success.csv"), row.names = FALSE)
  utils::write.csv(prog2$phase2_summary, file.path(output_dir, "program2_phase2_summary.csv"), row.names = FALSE)
  utils::write.csv(prog2$selection_summary, file.path(output_dir, "program2_selection_summary.csv"), row.names = FALSE)
  utils::write.csv(prog2$phase3_by_pop_conditional_on_go, file.path(output_dir, "program2_phase3_conditional_on_go.csv"), row.names = FALSE)
  utils::write.csv(prog2$program_success, file.path(output_dir, "program2_program_success.csv"), row.names = FALSE)

  invisible(list(
    data = dat,
    split = split,
    prior_checks = prior_checks,
    fit_hist = fit_hist,
    eta_prior = eta_prior,
    fit_emax = fit_emax,
    emax_draws = emax_draws,
    prog1 = prog1,
    prog2 = prog2
  ))
}

historical_stan_data <- function(arms, prior, active_weight = 0.5) {

  stopifnot(active_weight >= 0, active_weight <= 1)

  # HOPE is reserved for the next stage because it studies the same drug.
  historical <- arms |>
    dplyr::filter(title != "HOPE")

  studies <- unique(historical$nct)

  baseline <- historical |>
    dplyr::filter(use_baseline)

  # Historical doses do not enter this likelihood. Selected active arms are
  # discounted proxies for maximum response, rather than a pooled dose curve.
  list(
    N = nrow(historical),
    J = length(studies),
    study = match(historical$nct, studies),
    y = historical$chg,
    se = historical$chg_se,
    w_em = historical$w_em,
    active = as.integer(historical$arm == "active"),
    weight = ifelse(historical$arm == "active", active_weight, 1),
    NB = nrow(baseline),
    study_b = match(baseline$nct, studies),
    y_b = baseline$base_mmd,
    se_b = baseline$baseline_se,
    w_em_b = baseline$w_em,
    m = prior$m,
    s = prior$s,
    q = prior$q,
    lkj_shape = prior$lkj_shape
  )
}


fit_meta <- function(
    arms,
    prior,
    output_dir = "results_R/stan_meta",
    active_weight = 0.5,
    seed = 4301,
    ...) {

  # Prepare the arm-summary likelihood and fit the hierarchical Stan model.
  data <- historical_stan_data(arms, prior, active_weight)

  fit_stan("historical_meta.stan", data, output_dir, seed, ...)
}


predictive_prior_from_meta <- function(fit, covariance_multiplier = 1) {

  stopifnot(covariance_multiplier > 0)

  # eta_new is a NEW study draw. It includes between-study heterogeneity,
  # rather than only posterior uncertainty about the historical mean.
  eta <- draw_matrix(fit, paste0("eta_new[", 1:6, "]"))

  # Transfer the joint prediction through one Gaussian moment approximation.
  # Do not independently resample parameter marginals: dependence matters.
  list(
    mean = colMeans(eta),
    cov = stats::cov(eta) * covariance_multiplier + diag(1e-8, 6),
    clinical_draws = clinical_from_eta(eta)
  )
}

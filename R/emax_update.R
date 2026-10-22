make_ed50_prior <- function(q05 = sqrt(40), q50 = 20, q95 = sqrt(4000)) {

  stopifnot(0 < q05, q05 < q50, q50 < q95)

  # A lognormal median must be the geometric mean of symmetric tail quantiles.
  if (abs(log(q50) - mean(log(c(q05, q95)))) > 1e-6) {
    stop("For a lognormal prior, q50 must equal sqrt(q05 * q95).")
  }

  list(
    logmean = log(q50),
    logsd = log(q95 / q05) / (2 * stats::qnorm(0.95)),
    q05 = q05,
    q50 = q50,
    q95 = q95
  )
}


hope_stan_data <- function(arms, predictive_prior, ed50_prior) {

  hope <- arms |>
    dplyr::filter(title == "HOPE")

  baseline <- hope |>
    dplyr::filter(use_baseline)

  # Unlike historical cross-compound doses, these doses identify the new
  # drug's potency. Use each HOPE arm's own EM/CM composition.
  list(
    N = nrow(hope),
    y = hope$chg,
    se = hope$chg_se,
    dose = hope$dose,
    w_em = hope$w_em,
    NB = nrow(baseline),
    y_b = baseline$base_mmd,
    se_b = baseline$baseline_se,
    w_em_b = baseline$w_em,
    prior_mean = predictive_prior$mean,
    prior_cov = predictive_prior$cov,
    ed50_logmean = ed50_prior$logmean,
    ed50_logsd = ed50_prior$logsd
  )
}


fit_emax <- function(
    arms,
    predictive_prior,
    ed50_prior = make_ed50_prior(),
    output_dir = "results_R/stan_hope",
    seed = 4401,
    ...) {

  # HOPE jointly updates the historical clinical prior and separate ED50 prior.
  data <- hope_stan_data(arms, predictive_prior, ed50_prior)

  fit_stan("hope_emax.stan", data, output_dir, seed, ...)
}


emax_draws <- function(fit) {

  # Retain whole posterior rows so potency--effect dependence survives simulation.
  eta <- draw_matrix(fit, paste0("eta[", 1:6, "]"))

  clinical_from_eta(eta) |>
    dplyr::mutate(ed50 = as.numeric(draw_matrix(fit, "ed50"))) |>
    # Guarantee a CSV-ready result, including all eight clinical columns.
    as_numeric_draws()
}

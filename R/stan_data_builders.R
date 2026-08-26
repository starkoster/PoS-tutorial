#' Build Stan data for the historical meta-analysis
#'
#' This is a list-builder wrapper for users who want direct access to Stan data.
#' It accepts either the package's historical change data (`chg`, `chg_se`) or a
#' data frame with aliases `y`, `se`. Baseline data may use `base_y`, `base_se`
#' or the original CSV columns `base_mmd`, `base_mmd_sd`.
#'
#' @param hist_long Historical data from [make_historical_change_data()].
#' @param base_long Optional baseline data from [make_baseline_data()].
#' @param mu_eta_prior Eta prior mean vector of length 6.
#' @param sd_mu_eta_prior Eta prior SD vector of length 6.
#' @param tau_eta_scale Half-normal scales for between-study SDs.
#' @param lkj_shape LKJ correlation prior shape.
#' @return A list suitable for `historical_meta.stan`.
#' @export
build_historical_stan_data <- function(hist_long,
                                       base_long = NULL,
                                       mu_eta_prior,
                                       sd_mu_eta_prior,
                                       tau_eta_scale,
                                       lkj_shape = 2) {
  if (is.null(base_long)) {
    base_long <- data.frame(
      study_id = integer(0),
      base_y = numeric(0),
      base_se = numeric(0),
      w_em = numeric(0)
    )
  }

  y <- if ("y" %in% names(hist_long)) hist_long$y else hist_long$chg
  se <- if ("se" %in% names(hist_long)) hist_long$se else hist_long$chg_se

  list(
    N_change = nrow(hist_long),
    J = max(hist_long$study_id),
    study_change = as.integer(hist_long$study_id),
    y_change = as.numeric(y),
    se_change = as.numeric(se),
    w_em = as.numeric(hist_long$w_em),
    active_ref = as.numeric(hist_long$active_ref),
    weight_change = as.numeric(hist_long$weight),
    N_base = nrow(base_long),
    study_base = as.integer(base_long$study_id),
    y_base = as.numeric(base_long$base_y),
    se_base = as.numeric(base_long$base_se),
    w_em_base = as.numeric(base_long$w_em),
    mu_eta_prior = as.numeric(mu_eta_prior),
    sd_mu_eta_prior = as.numeric(sd_mu_eta_prior),
    tau_eta_scale = as.numeric(tau_eta_scale),
    lkj_shape = lkj_shape
  )
}

#' Build Stan data for the new-drug ED50 Emax model
#'
#' Baseline rows from the new trial are included when available. They inform the
#' new-trial baseline parameters and therefore the maximum possible placebo and
#' active reductions under the 0--30 MMD constraints.
#'
#' @param new_long New-study data from [make_hope_emax_data()].
#' @param eta_prior_mu Eta prior mean vector from [extract_eta_prior()].
#' @param eta_prior_Sigma Eta prior covariance matrix from [extract_eta_prior()].
#' @param ed50_prior List from [make_ed50_prior()].
#' @param dose_grid Prediction grid.
#' @param base_long Optional baseline data from [make_baseline_data()].
#' @return A list suitable for `historical_emax.stan`.
#' @export
build_emax_stan_data <- function(new_long,
                                 eta_prior_mu,
                                 eta_prior_Sigma,
                                 ed50_prior,
                                 dose_grid = NULL,
                                 base_long = NULL) {
  if (is.null(dose_grid)) {
    max_dose <- max(new_long$dose, na.rm = TRUE)
    dose_grid <- seq(0, max_dose * 1.5, length.out = 100)
  }
  if (is.null(base_long)) {
    base_long <- make_baseline_data(new_long)
  }

  y <- if ("y" %in% names(new_long)) new_long$y else new_long$chg
  se <- if ("se" %in% names(new_long)) new_long$se else new_long$chg_se

  list(
    N = nrow(new_long),
    y = as.numeric(y),
    se = as.numeric(se),
    dose = as.numeric(new_long$dose),
    w_em = as.numeric(new_long$w_em),
    N_base = nrow(base_long),
    y_base = as.numeric(base_long$base_y),
    se_base = as.numeric(base_long$base_se),
    w_em_base = as.numeric(base_long$w_em),
    eta_prior_mu = as.numeric(eta_prior_mu),
    eta_prior_Sigma = eta_prior_Sigma,
    ed50_prior_logmean = ed50_prior$logmean,
    ed50_prior_logsd = ed50_prior$logsd,
    N_pred = length(dose_grid),
    dose_pred = as.numeric(dose_grid)
  )
}

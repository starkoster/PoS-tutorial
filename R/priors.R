#' Normalize and validate clinical prior medians
#'
#' This helper is used by [default_eta_priors()] and [eta_mu_from_medians()].
#' Priors may be specified either with total active/reference reductions
#' (`active_red_em`, `active_red_cm`) or directly with placebo-adjusted drug
#' reductions (`drug_red_em`, `drug_red_cm`). If drug reductions are supplied,
#' total active reductions are computed as `pbo_red + drug_red`.
#'
#' Required clinical quantities are baseline MMDs and placebo reductions.
#' The returned list always contains both `drug_red_*` and `active_red_*`.
#'
#' @param med Named list with clinical prior medians.
#' @return A normalized named list of prior medians.
#' @export
normalize_prior_medians <- function(med) {
  required_basic <- c("baseline_em", "baseline_cm", "pbo_red_em", "pbo_red_cm")
  missing_basic <- setdiff(required_basic, names(med))
  if (length(missing_basic) > 0) {
    stop(
      "Missing prior median entries: ", paste(missing_basic, collapse = ", "),
      call. = FALSE
    )
  }

  has_active <- all(c("active_red_em", "active_red_cm") %in% names(med))
  has_drug <- all(c("drug_red_em", "drug_red_cm") %in% names(med))

  if (!has_active && !has_drug) {
    stop(
      "Specify either active_red_em/active_red_cm or drug_red_em/drug_red_cm.",
      call. = FALSE
    )
  }

  out <- as.list(med)

  numeric_names <- intersect(
    c(
      "baseline_em", "baseline_cm",
      "pbo_red_em", "pbo_red_cm",
      "active_red_em", "active_red_cm",
      "drug_red_em", "drug_red_cm"
    ),
    names(out)
  )
  for (nm in numeric_names) {
    out[[nm]] <- as.numeric(out[[nm]])
  }

  if (has_drug && !has_active) {
    out$active_red_em <- out$pbo_red_em + out$drug_red_em
    out$active_red_cm <- out$pbo_red_cm + out$drug_red_cm
  } else if (has_active && !has_drug) {
    out$drug_red_em <- out$active_red_em - out$pbo_red_em
    out$drug_red_cm <- out$active_red_cm - out$pbo_red_cm
  } else if (has_active && has_drug) {
    inferred_drug_em <- out$active_red_em - out$pbo_red_em
    inferred_drug_cm <- out$active_red_cm - out$pbo_red_cm
    tol <- 1e-8
    if (abs(inferred_drug_em - out$drug_red_em) > tol ||
        abs(inferred_drug_cm - out$drug_red_cm) > tol) {
      stop(
        "Inconsistent prior medians: active_red must equal pbo_red + drug_red. ",
        "Either remove active_red_* and let the function compute it, or make ",
        "active_red_em = pbo_red_em + drug_red_em and active_red_cm = pbo_red_cm + drug_red_cm.",
        call. = FALSE
      )
    }
  }

  vals <- unlist(out[c(
    "baseline_em", "baseline_cm",
    "pbo_red_em", "pbo_red_cm",
    "drug_red_em", "drug_red_cm",
    "active_red_em", "active_red_cm"
  )])
  if (any(!is.finite(vals))) {
    stop("All clinical prior medians must be finite numeric values.", call. = FALSE)
  }

  b_em <- out$baseline_em
  b_cm <- out$baseline_cm
  p_em <- out$pbo_red_em
  p_cm <- out$pbo_red_cm
  a_em <- out$drug_red_em
  a_cm <- out$drug_red_cm
  r_em <- out$active_red_em
  r_cm <- out$active_red_cm

  if (!(0 < b_em && b_em < b_cm && b_cm < 30)) {
    stop(
      "Need 0 < baseline_em < baseline_cm < 30. Current values: ",
      "baseline_em = ", b_em, ", baseline_cm = ", b_cm, ".",
      call. = FALSE
    )
  }
  if (!(0 < p_em && p_em < p_cm)) {
    stop(
      "Need 0 < pbo_red_em < pbo_red_cm. Current values: ",
      "pbo_red_em = ", p_em, ", pbo_red_cm = ", p_cm, ".",
      call. = FALSE
    )
  }
  if (!(0 < a_em && a_em < a_cm)) {
    stop(
      "Need 0 < drug_red_em < drug_red_cm, i.e. CM has the larger placebo-adjusted drug reduction. ",
      "Current values: drug_red_em = ", a_em, ", drug_red_cm = ", a_cm, ". ",
      "Remember drug_red = active_red - pbo_red if active_red_* is used.",
      call. = FALSE
    )
  }
  if (!(0 < r_em && r_em < b_em)) {
    stop(
      "Need 0 < active_red_em < baseline_em. Current values: ",
      "active_red_em = ", r_em, ", baseline_em = ", b_em, ".",
      call. = FALSE
    )
  }
  if (!(r_em < r_cm && r_cm < b_cm)) {
    stop(
      "Need active_red_em < active_red_cm < baseline_cm. Current values: ",
      "active_red_em = ", r_em, ", active_red_cm = ", r_cm,
      ", baseline_cm = ", b_cm, ".",
      call. = FALSE
    )
  }
  if (!(p_em < r_em && p_cm < r_cm)) {
    stop(
      "Need placebo reductions smaller than total active reductions. Current values: ",
      "pbo_red_em = ", p_em, ", active_red_em = ", r_em,
      ", pbo_red_cm = ", p_cm, ", active_red_cm = ", r_cm, ".",
      call. = FALSE
    )
  }

  out
}

#' Convert interpretable prior medians to unconstrained eta means
#'
#' The historical meta-analysis is parameterized through six unconstrained
#' quantities, `eta[1:6]`, but all subjective elicitation is done on a clinical
#' monthly-migraine-day scale.
#'
#' Medians are defined on a positive-reduction scale. For population `k`,
#' `pbo_red_k` is the placebo reduction, `drug_red_k` is the placebo-adjusted
#' drug reduction, and `active_red_k = pbo_red_k + drug_red_k` is the total
#' active-arm reduction from baseline. The function accepts either active
#' reductions or drug reductions; see [normalize_prior_medians()].
#'
#' @param med Named list of clinical prior medians.
#' @return Numeric eta mean vector of length 6.
#' @export
eta_mu_from_medians <- function(med) {
  med <- normalize_prior_medians(med)

  b_em <- med$baseline_em
  b_cm <- med$baseline_cm
  r_em <- med$active_red_em
  r_cm <- med$active_red_cm
  p_em <- med$pbo_red_em
  p_cm <- med$pbo_red_cm

  c(
    stats::qlogis(b_em / 30),
    stats::qlogis((b_cm - b_em) / (30 - b_em)),
    stats::qlogis(r_em / b_em),
    stats::qlogis((r_cm - r_em) / (b_cm - r_em)),
    stats::qlogis(p_em / r_em),
    stats::qlogis((p_cm - p_em) / (r_cm - r_em))
  )
}

#' Create eta priors from clinical medians and weakly informative scales
#'
#' This is the only subjective prior layer in the historical meta-analysis.
#' The medians provide a clinically interpretable center, while
#' `sd_mu_eta_prior` controls how strongly the meta-analytic average is pulled
#' toward that center and `tau_eta_scale` controls between-study heterogeneity.
#'
#' The default medians encode:
#' - EM baseline MMD around 8 and CM baseline MMD around 18;
#' - EM placebo reduction around 1 MMD and CM placebo reduction around 2 MMD;
#' - EM placebo-adjusted drug reduction around 2 MMD and CM placebo-adjusted
#'   drug reduction around 3 MMD.
#'
#' These defaults imply total active/reference reductions of 3 MMD in EM and
#' 5 MMD in CM. The hard transformation constrains baseline MMD to 0--30,
#' CM baseline > EM baseline, CM placebo response > EM placebo response, CM
#' drug effect > EM drug effect, and total reductions below baseline.
#'
#' @param prior_medians Optional list of interpretable medians. You may specify
#' `drug_red_em` and `drug_red_cm` directly, or specify `active_red_em` and
#' `active_red_cm`; if both are supplied they must be consistent.
#' @param sd_mu_eta_prior Prior SDs for the meta-analytic eta means. Values
#' around 0.6 are weak-to-moderately informative on the logit scale; 0.75--1.0
#' is weaker and should be checked with prior predictive simulation.
#' @param tau_eta_scale Half-normal scales for between-study eta SDs. Values
#' around 0.4--0.6 usually allow meaningful heterogeneity without putting too
#' much prior predictive mass near impossible extremes.
#' @return A list of Stan prior inputs and clinical medians.
#' @export
default_eta_priors <- function(prior_medians = NULL,
                               sd_mu_eta_prior = rep(0.60, 6),
                               tau_eta_scale = c(0.40, 0.40, 0.45, 0.45, 0.45, 0.45)) {
  if (is.null(prior_medians)) {
    prior_medians <- list(
      baseline_em = 10,
      baseline_cm = 15,
      pbo_red_em = 0.5,
      pbo_red_cm = 1.0,
      drug_red_em = 0.5,
      drug_red_cm = 1.0
    )
  }

  prior_medians <- normalize_prior_medians(prior_medians)

  list(
    prior_medians = prior_medians,
    mu_eta_prior = eta_mu_from_medians(prior_medians),
    sd_mu_eta_prior = as.numeric(sd_mu_eta_prior),
    tau_eta_scale = as.numeric(tau_eta_scale)
  )
}

#' Create a very weak eta prior template
#'
#' This helper is useful for sensitivity analyses. It preserves the same
#' clinical medians and hard constraints as [default_eta_priors()], but uses
#' wider eta-scale priors. Always run [simulate_eta_prior()] and
#' [prior_predictive_summary()] after using this helper.
#'
#' @param prior_medians Optional list of interpretable medians.
#' @return A prior list for Stan.
#' @export
weak_eta_priors <- function(prior_medians = NULL) {
  default_eta_priors(
    prior_medians = prior_medians,
    sd_mu_eta_prior = rep(0.75, 6),
    tau_eta_scale = rep(0.60, 6)
  )
}

#' Human-readable eta prior specification
#'
#' @param priors Output from [default_eta_priors()] or [weak_eta_priors()].
#' @return A data frame describing the six eta prior components and the core
#' clinical prior medians.
#' @export
eta_prior_table <- function(priors) {
  med <- normalize_prior_medians(priors$prior_medians)
  clinical_median <- c(
    med$baseline_em,
    med$baseline_cm - med$baseline_em,
    med$active_red_em / med$baseline_em,
    (med$active_red_cm - med$active_red_em) / (med$baseline_cm - med$active_red_em),
    med$pbo_red_em / med$active_red_em,
    (med$pbo_red_cm - med$pbo_red_em) / (med$active_red_cm - med$active_red_em)
  )

  data.frame(
    eta = paste0("eta", 1:6),
    clinical_quantity = c(
      "EM baseline MMD as fraction of 30",
      "CM baseline excess as fraction of remaining MMD range",
      "EM total active/reference reduction as fraction of EM baseline",
      "CM total active/reference excess as fraction of available CM range",
      "EM placebo reduction as fraction of EM total active/reference reduction",
      "CM placebo excess as fraction of EM-to-CM total active/reference excess"
    ),
    clinical_median_input = clinical_median,
    eta_mean = as.numeric(priors$mu_eta_prior),
    eta_mean_sd = as.numeric(priors$sd_mu_eta_prior),
    between_study_eta_sd_halfnormal_scale = as.numeric(priors$tau_eta_scale),
    baseline_em_median = med$baseline_em,
    baseline_cm_median = med$baseline_cm,
    pbo_red_em_median = med$pbo_red_em,
    pbo_red_cm_median = med$pbo_red_cm,
    drug_red_em_median = med$drug_red_em,
    drug_red_cm_median = med$drug_red_cm,
    active_red_em_median = med$active_red_em,
    active_red_cm_median = med$active_red_cm
  )
}

#' Transform eta draws to interpretable parameters
#'
#' @param eta Vector, matrix, or data frame with six eta columns.
#' @return A data frame with baseline, placebo reduction, drug reduction, Emax,
#' and post-treatment MMD quantities.
#' @export
eta_to_pars <- function(eta) {
  eta <- as.matrix(eta)
  if (ncol(eta) != 6) {
    stop("eta must have six columns.", call. = FALSE)
  }

  baseline_em <- 30 * stats::plogis(eta[, 1])
  baseline_cm <- baseline_em + (30 - baseline_em) * stats::plogis(eta[, 2])

  active_red_em <- baseline_em * stats::plogis(eta[, 3])
  active_red_cm <- active_red_em + (baseline_cm - active_red_em) * stats::plogis(eta[, 4])

  pbo_red_em <- active_red_em * stats::plogis(eta[, 5])
  pbo_red_cm <- pbo_red_em + (active_red_cm - active_red_em) * stats::plogis(eta[, 6])

  drug_red_em <- active_red_em - pbo_red_em
  drug_red_cm <- active_red_cm - pbo_red_cm

  data.frame(
    baseline_em = baseline_em,
    baseline_cm = baseline_cm,
    pbo_red_em = pbo_red_em,
    pbo_red_cm = pbo_red_cm,
    drug_red_em = drug_red_em,
    drug_red_cm = drug_red_cm,
    active_red_em = active_red_em,
    active_red_cm = active_red_cm,
    e0_em = -pbo_red_em,
    e0_cm = -pbo_red_cm,
    emax_em = -drug_red_em,
    emax_cm = -drug_red_cm,
    post_pbo_em = baseline_em - pbo_red_em,
    post_pbo_cm = baseline_cm - pbo_red_cm,
    post_active_em = baseline_em - active_red_em,
    post_active_cm = baseline_cm - active_red_cm
  )
}

#' Simulate from the eta prior predictive distribution
#'
#' This uses independent half-normal between-study standard deviations on the
#' eta scale for quick calibration. The Stan model adds an LKJ correlation prior.
#'
#' @param priors Output from [default_eta_priors()] or [weak_eta_priors()].
#' @param n Number of prior-predictive draws.
#' @param seed Random seed.
#' @return Interpretable prior-predictive draws.
#' @export
simulate_eta_prior <- function(priors, n = 200000, seed = 123) {
  set.seed(seed)
  mu_eta_prior <- priors$mu_eta_prior
  sd_mu_eta_prior <- priors$sd_mu_eta_prior
  tau_eta_scale <- priors$tau_eta_scale
  K <- length(mu_eta_prior)

  mu_eta <- matrix(stats::rnorm(n * K), nrow = n, ncol = K)
  mu_eta <- sweep(mu_eta, 2, sd_mu_eta_prior, "*")
  mu_eta <- sweep(mu_eta, 2, mu_eta_prior, "+")

  tau_eta <- abs(matrix(stats::rnorm(n * K), nrow = n, ncol = K))
  tau_eta <- sweep(tau_eta, 2, tau_eta_scale, "*")

  z <- matrix(stats::rnorm(n * K), nrow = n, ncol = K)
  eta_new <- mu_eta + tau_eta * z

  eta_to_pars(eta_new)
}

#' Summarize prior-predictive tail probabilities
#'
#' @param pars Output of [simulate_eta_prior()] or [eta_to_pars()].
#' @param threshold Reduction threshold, default 15 MMD.
#' @return A one-row data frame of key prior checks.
#' @export
prior_predictive_summary <- function(pars, threshold = 15) {
  data.frame(
    pr_drug_red_em_gt_threshold = mean(pars$drug_red_em > threshold),
    pr_drug_red_cm_gt_threshold = mean(pars$drug_red_cm > threshold),
    pr_active_red_em_gt_threshold = mean(pars$active_red_em > threshold),
    pr_active_red_cm_gt_threshold = mean(pars$active_red_cm > threshold),
    min_post_active_em = min(pars$post_active_em),
    min_post_active_cm = min(pars$post_active_cm),
    max_baseline_em = max(pars$baseline_em),
    max_baseline_cm = max(pars$baseline_cm)
  )
}

#' Make an informative ED50 prior
#'
#' @param median Prior median ED50.
#' @param lower Lower plausible ED50 bound.
#' @param upper Upper plausible ED50 bound.
#' @param mass Probability mass between lower and upper.
#' @return A lognormal prior specification.
#' @export
make_ed50_prior <- function(median, lower, upper, mass = 0.90) {
  if (!(median > 0 && lower > 0 && upper > 0)) {
    stop("median, lower, and upper must be positive.", call. = FALSE)
  }
  if (!(lower < median && median < upper)) {
    stop("Need lower < median < upper.", call. = FALSE)
  }
  if (!(mass > 0 && mass < 1)) {
    stop("mass must be between 0 and 1.", call. = FALSE)
  }

  z <- stats::qnorm((1 + mass) / 2)
  logsd_lower <- (log(median) - log(lower)) / z
  logsd_upper <- (log(upper) - log(median)) / z

  list(
    logmean = log(median),
    logsd = mean(c(logsd_lower, logsd_upper)),
    median = median,
    lower = lower,
    upper = upper,
    mass = mass
  )
}

#' ED50 Emax dose score
#'
#' @param dose Dose value.
#' @param ed50 ED50 value.
#' @return Dose score `dose / (dose + ed50)`.
#' @export
emax_g <- function(dose, ed50) {
  dose / (dose + ed50)
}

#' Population-specific mean change from baseline
#'
#' @param draw One posterior draw with `e0_em`, `emax_em`, `e0_cm`, `emax_cm`, and `ed50`.
#' @param pop Either `"EM"` or `"CM"`.
#' @param dose Dose value.
#' @return Mean change from baseline.
#' @export
mean_emax_pop <- function(draw, pop, dose) {
  g <- emax_g(dose, draw[["ed50"]])
  if (pop == "EM") {
    draw[["e0_em"]] + draw[["emax_em"]] * g
  } else if (pop == "CM") {
    draw[["e0_cm"]] + draw[["emax_cm"]] * g
  } else {
    stop("pop must be 'EM' or 'CM'.", call. = FALSE)
  }
}

#' Mixed-population mean change from baseline
#'
#' @param draw One posterior draw.
#' @param dose Dose value.
#' @param prop_em Proportion EM.
#' @return Mixed-population mean change from baseline.
#' @export
mean_emax_mix <- function(draw, dose, prop_em) {
  prop_em * mean_emax_pop(draw, "EM", dose) +
    (1 - prop_em) * mean_emax_pop(draw, "CM", dose)
}

#' Population-specific positive drug effect versus placebo
#'
#' @param draw One posterior draw.
#' @param pop Either `"EM"` or `"CM"`.
#' @param dose Dose value.
#' @return Positive placebo-adjusted reduction.
#' @export
true_effect_pop <- function(draw, pop, dose) {
  mean_emax_pop(draw, pop, 0) - mean_emax_pop(draw, pop, dose)
}

#' Mixed-population positive drug effect versus placebo
#'
#' @param draw One posterior draw.
#' @param dose Dose value.
#' @param prop_em Proportion EM.
#' @return Positive placebo-adjusted reduction.
#' @export
true_effect_mix <- function(draw, dose, prop_em) {
  mean_emax_mix(draw, 0, prop_em) - mean_emax_mix(draw, dose, prop_em)
}

#' Active-versus-placebo normal approximation test
#'
#' Negative active-minus-placebo difference is beneficial because the endpoint is
#' change from baseline in monthly migraine days.
#'
#' @param y_pbo Observed placebo arm mean.
#' @param y_act Observed active arm mean.
#' @param n_pbo Placebo sample size.
#' @param n_act Active sample size.
#' @param sigma Subject-level standard deviation.
#' @param alpha_two_sided Two-sided alpha.
#' @return List with statistics and success flag.
#' @export
active_vs_placebo_test <- function(y_pbo,
                                   y_act,
                                   n_pbo,
                                   n_act,
                                   sigma = 5,
                                   alpha_two_sided = 0.05) {
  diff <- y_act - y_pbo
  se_diff <- sigma * sqrt(1 / n_pbo + 1 / n_act)
  z <- diff / se_diff
  p_two <- 2 * stats::pnorm(-abs(z))
  list(
    diff = diff,
    obs_effect = y_pbo - y_act,
    z = z,
    p_two = p_two,
    significant = isTRUE(diff < 0 && p_two < alpha_two_sided)
  )
}

#' Make an Emax-shaped MCP-style dose-response test
#'
#' The test takes a vector of arm means and evaluates a maximum over contrasts
#' corresponding to candidate ED50 values. The null critical value is simulated.
#'
#' @param doses Dose vector including placebo.
#' @param n_per_arm Subjects per arm.
#' @param sigma Subject-level standard deviation.
#' @param ed50_candidates Positive ED50 candidates.
#' @param alpha One-sided dose-response alpha.
#' @param n_null Number of null simulations.
#' @param seed Random seed.
#' @return A function that tests observed arm means.
#' @export
make_emax_mcp_test <- function(doses,
                               n_per_arm,
                               sigma = 5,
                               ed50_candidates,
                               alpha = 0.05,
                               n_null = 200000,
                               seed = 991) {
  if (length(doses) < 3) {
    stop("At least three dose levels are needed.", call. = FALSE)
  }
  if (any(ed50_candidates <= 0)) {
    stop("ed50_candidates must be positive.", call. = FALSE)
  }

  cmat <- sapply(ed50_candidates, function(ed50) {
    x <- doses / (doses + ed50)
    cvec <- -(x - mean(x))
    cvec / sqrt(sum(cvec^2))
  })
  if (is.null(dim(cmat))) {
    cmat <- matrix(cmat, ncol = 1)
  }
  colnames(cmat) <- paste0("ed50_", signif(ed50_candidates, 4))

  corr <- crossprod(cmat)
  corr <- (corr + t(corr)) / 2
  eig <- eigen(corr, symmetric = TRUE)
  vals <- pmax(eig$values, 0)
  root <- eig$vectors %*% diag(sqrt(vals), nrow = length(vals))

  set.seed(seed)
  z <- matrix(stats::rnorm(n_null * ncol(cmat)), nrow = n_null, ncol = ncol(cmat)) %*% t(root)
  max_null <- apply(z, 1, max)
  crit <- unname(stats::quantile(max_null, probs = 1 - alpha, names = FALSE))
  se_arm <- sigma / sqrt(n_per_arm)

  function(ybar) {
    if (length(ybar) != length(doses)) {
      stop("ybar must have the same length as doses.", call. = FALSE)
    }
    tvals <- as.numeric(crossprod(cmat, ybar) / se_arm)
    names(tvals) <- colnames(cmat)
    tmax <- max(tvals)
    best <- which.max(tvals)
    list(
      significant = tmax > crit,
      p_adj = mean(max_null >= tmax),
      tmax = tmax,
      best_ed50_candidate = ed50_candidates[best],
      Tvals = tvals,
      critical_value = crit
    )
  }
}

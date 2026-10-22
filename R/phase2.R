make_mcp_test <- function(
    doses,
    ed50_candidates,
    n_per_arm = 140L,
    sigma = 5,
    alpha = 0.05,
    seed = 4501) {

  stopifnot(
    length(doses) >= 3,
    sum(doses == 0) == 1,
    !anyDuplicated(doses),
    all(is.finite(doses) & doses >= 0),
    all(is.finite(ed50_candidates) & ed50_candidates > 0),
    length(ed50_candidates) >= 1,
    n_per_arm > 0,
    sigma > 0,
    alpha > 0 & alpha < 1
  )

  doses <- sort(doses)

  # Negative maxEff tells DoseFinding that a DECREASING response is beneficial.
  models <- DoseFinding::Mods(
    emax = unique(as.numeric(ed50_candidates)),
    doses = doses,
    placEff = 0,
    maxEff = -1
  )

  # Simulated arm means have known variance sigma^2 / n and independent errors.
  covariance <- diag(sigma^2 / n_per_arm, length(doses))

  # The critical value depends on the design and candidate shapes, not on the
  # future observations. Calculate it ONCE using the package's MCP procedure.
  set.seed(seed)

  reference_test <- DoseFinding::MCTtest(
    dose = doses,
    resp = rep(0, length(doses)),
    models = models,
    S = covariance,
    type = "general",
    df = Inf,
    alpha = alpha,
    alternative = "one.sided",
    critV = TRUE,
    pVal = FALSE
  )

  list(
    doses = doses,
    models = models,
    contrasts = reference_test$contMat,
    covariance = covariance,
    critical_value = as.numeric(reference_test$critVal),
    n_per_arm = n_per_arm,
    sigma = sigma,
    alpha = alpha
  )
}


simulate_phase2 <- function(draw, mcp, prop_em = 0.5, go_delta = -2) {

  stopifnot(prop_em >= 0, prop_em <= 1, go_delta < 0)

  # Each arm has fixed EM/CM counts. Rounding makes the actual fraction explicit.
  n_em <- round(mcp$n_per_arm * prop_em)
  weight_em <- n_em / mcp$n_per_arm

  expected <- weight_em * mean_change(draw, mcp$doses, "EM") +
    (1 - weight_em) * mean_change(draw, mcp$doses, "CM")

  # Simulate one mean for each arm. Do not reset the RNG inside trial functions.
  observed <- stats::rnorm(
    n = length(mcp$doses),
    mean = expected,
    sd = mcp$sigma / sqrt(mcp$n_per_arm)
  )

  # DoseFinding performs the multiple-contrast test. Reuse its design-specific
  # contrasts and critical value; no custom null simulation is required.
  test <- DoseFinding::MCTtest(
    dose = mcp$doses,
    resp = observed,
    models = mcp$models,
    S = mcp$covariance,
    type = "general",
    df = Inf,
    alpha = mcp$alpha,
    alternative = "one.sided",
    contMat = mcp$contrasts,
    critV = mcp$critical_value,
    pVal = FALSE
  )

  significant <- max(test$tStat) > mcp$critical_value

  # Every active comparison uses the same observed placebo mean.
  arms <- tibble::tibble(
    dose = mcp$doses,
    expected = expected,
    observed = observed,
    delta = observed - observed[mcp$doses == 0]
  )

  qualifying <- arms |>
    dplyr::filter(dose > 0, delta <= go_delta)

  # Go requires statistical evidence AND a sufficiently large observed effect.
  go <- significant && nrow(qualifying) > 0
  selected_dose <- if (go) min(qualifying$dose) else NA_real_

  list(
    go = go,
    selected_dose = selected_dose,
    significant = significant,
    t_max = max(test$tStat),
    arms = arms
  )
}

simulate_phase3 <- function(
    draw,
    doses,
    n_per_arm,
    sigma = 5,
    alpha = 0.05,
    multiplicity = c("none", "bonferroni")) {

  multiplicity <- match.arg(multiplicity)

  stopifnot(
    sum(doses == 0) == 1,
    any(doses > 0),
    !anyDuplicated(doses),
    all(is.finite(doses) & doses >= 0),
    n_per_arm > 0,
    sigma > 0,
    alpha > 0 & alpha < 1
  )

  active <- doses > 0

  # The default preserves the previous nominal benchmark. Bonferroni is an
  # explicit alternative for the active-dose comparisons within each trial.
  test_alpha <- if (multiplicity == "bonferroni") alpha / sum(active) else alpha

  purrr::map_dfr(c("EM", "CM"), function(population) {

    # EM and CM are separate trials with independent sampling errors given truth.
    expected <- mean_change(draw, doses, population)

    observed <- stats::rnorm(
      n = length(doses),
      mean = expected,
      sd = sigma / sqrt(n_per_arm)
    )

    # Reuse one placebo mean within a trial, preserving dose-comparison dependence.
    placebo_mean <- observed[doses == 0]
    delta <- observed[active] - placebo_mean
    se_delta <- sigma * sqrt(2 / n_per_arm)

    p_two <- 2 * stats::pnorm(-abs(delta / se_delta))

    tibble::tibble(
      population = population,
      dose = doses[active],
      placebo_mean = placebo_mean,
      active_mean = observed[active],
      delta = delta,
      true_delta = expected[active] - expected[doses == 0],
      p_two = p_two,
      significant = delta < 0 & p_two < test_alpha
    )
  })
}


phase3_success <- function(phase3, threshold = -2) {

  # A population succeeds if any active dose is significant and meets the TPP.
  em_doses <- phase3 |>
    dplyr::filter(population == "EM", significant, delta <= threshold) |>
    dplyr::pull(dose)

  cm_doses <- phase3 |>
    dplyr::filter(population == "CM", significant, delta <= threshold) |>
    dplyr::pull(dose)

  tibble::tibble(
    EM = length(em_doses) > 0,
    CM = length(cm_doses) > 0,
    both = length(em_doses) > 0 && length(cm_doses) > 0,
    both_same_dose = length(intersect(em_doses, cm_doses)) > 0
  )
}

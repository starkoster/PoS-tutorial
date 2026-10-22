probability_summary <- function(success) {

  # MCSE describes simulation noise conditional on the available posterior draws.
  n <- length(success)
  probability <- if (n > 0) mean(success) else NA_real_

  tibble::tibble(
    n = n,
    successes = sum(success),
    probability = probability,
    mcse = if (n > 0) sqrt(probability * (1 - probability) / n) else NA_real_
  )
}


summarize_programs <- function(simulation) {

  # A long table lets the same summary cover EM, CM and joint success.
  events <- simulation$events |>
    tidyr::pivot_longer(
      cols = c(EM, CM, both, both_same_dose),
      names_to = "population",
      values_to = "success"
    )

  probabilities <- events |>
    dplyr::group_by(strategy, threshold, population) |>
    dplyr::group_modify(function(data, key) {

      # Keep the denominator explicit: all programs versus advancing programs.
      dplyr::bind_rows(
        dplyr::mutate(probability_summary(data$success), estimand = "from_start"),
        dplyr::mutate(
          probability_summary(data$success[data$go]),
          estimand = "given_go"
        )
      )
    }) |>
    dplyr::ungroup()

  mean_or_na <- function(x) if (length(x) > 0) mean(x) else NA_real_

  # The paired screening comparison holds Program 1's Phase 3 design fixed.
  derisking <- events |>
    dplyr::filter(strategy == "P1_direct") |>
    dplyr::group_by(threshold, population) |>
    dplyr::summarise(
      p_go = mean(go),
      direct = mean(success),
      given_go = mean_or_na(success[go]),
      given_no_go = mean_or_na(success[!go]),
      gated_from_start = mean(success & go),
      .groups = "drop"
    ) |>
    dplyr::mutate(conditional_gain = given_go - direct)

  # Expected enrollment includes Phase 2 and only the Phase 3 trials that run.
  design <- simulation$design
  p_go <- mean(simulation$phase2$go)

  phase2_n <- length(design$p2_doses) * design$p2_n
  direct_n <- 2 * length(design$p1_doses) * design$p1_n
  selected_n <- 4 * design$p2_program_p3_n

  participants <- tibble::tibble(
    strategy = c("P1_direct", "P1_screened", "P2_selected"),
    expected_n = c(
      direct_n,
      phase2_n + p_go * direct_n,
      phase2_n + p_go * selected_n
    )
  )

  selection <- simulation$phase2 |>
    dplyr::filter(go) |>
    dplyr::count(selected_dose, name = "count") |>
    tidyr::complete(
      selected_dose = design$p2_doses[design$p2_doses > 0],
      fill = list(count = 0L)
    ) |>
    dplyr::mutate(fraction_of_go = if (sum(count) > 0) count / sum(count) else NA_real_)

  list(
    probabilities = probabilities,
    derisking = derisking,
    phase2 = probability_summary(simulation$phase2$go),
    selection = selection,
    participants = participants
  )
}

default_program_design <- function() {

  # Keep all prospective choices together so a reviewer can change the design.
  list(
    p2_doses = c(0, 20, 100, 300),
    p2_n = 140L,
    p2_em = 0.5,
    p2_go_delta = -2,
    p2_alpha = 0.05,

    p1_doses = c(0, 20, 200),
    p1_n = 233L,
    p2_program_p3_n = 174L,

    sigma = 5,
    p3_alpha = 0.05,
    p3_multiplicity = "none",
    thresholds = c(0, -2, -3, -4)
  )
}


program_events <- function(trials, phase2, thresholds, strategy) {

  # Evaluate every requested TPP threshold against the simulated Phase 3 results.
  # This is bookkeeping after trial generation, rather than another simulation.
  if (nrow(trials) > 0) {

    expanded_trials <- purrr::map_dfr(thresholds, function(target) {
      dplyr::mutate(trials, threshold = target)
    })

    observed_success <- expanded_trials |>
      dplyr::mutate(success = significant & delta <= threshold) |>
      dplyr::group_by(program, threshold) |>
      dplyr::summarise(
        EM = any(success[population == "EM"]),
        CM = any(success[population == "CM"]),
        both_same_dose = length(intersect(
          dose[success & population == "EM"],
          dose[success & population == "CM"]
        )) > 0,
        .groups = "drop"
      )

  } else {

    observed_success <- tibble::tibble(
      program = integer(),
      threshold = double(),
      EM = logical(),
      CM = logical(),
      both_same_dose = logical()
    )
  }

  # Include stopped programs in the denominator, with all success flags FALSE.
  tidyr::expand_grid(program = phase2$program, threshold = thresholds) |>
    dplyr::left_join(observed_success, by = c("program", "threshold")) |>
    dplyr::left_join(dplyr::select(phase2, program, go), by = "program") |>
    dplyr::mutate(
      dplyr::across(c(EM, CM, both_same_dose), ~ tidyr::replace_na(.x, FALSE)),
      both = EM & CM,
      strategy = strategy
    )
}


simulate_programs <- function(
    draws,
    mcp,
    n_sim = 50000L,
    seed = 4601,
    design = default_program_design()) {

  stopifnot(
    n_sim >= 1,
    n_sim == as.integer(n_sim),
    nrow(draws) >= 1,
    identical(as.numeric(mcp$doses), as.numeric(sort(design$p2_doses))),
    mcp$n_per_arm == design$p2_n,
    mcp$sigma == design$sigma,
    mcp$alpha == design$p2_alpha
  )

  stopifnot(
    all(c("placebo_em", "placebo_cm", "emax_em", "emax_cm", "ed50") %in% names(draws)),
    all(is.finite(draws$ed50) & draws$ed50 > 0)
  )

  # Sample whole posterior rows. One truth is shared by all trials in a program.
  set.seed(seed)
  posterior_rows <- sample.int(nrow(draws), n_sim, replace = TRUE)

  phase2_records <- vector("list", n_sim)
  phase2_arms <- vector("list", n_sim)
  direct_records <- vector("list", n_sim)
  selected_records <- vector("list", n_sim)

  for (program in seq_len(n_sim)) {

    truth <- draws[posterior_rows[program], , drop = FALSE]

    # STEP 1: simulate Phase 2 for Program 2 and Program 1's screening analysis.
    # The primary direct-to-Phase-3 Program 1 does not use this result.
    phase2 <- simulate_phase2(truth, mcp, design$p2_em, design$p2_go_delta)

    # STEP 2A: generate the fixed Program 1 Phase 3 outcomes for every truth.
    # After No Go these are counterfactual outcomes used only for paired screening.
    direct <- simulate_phase3(
      truth,
      doses = design$p1_doses,
      n_per_arm = design$p1_n,
      sigma = design$sigma,
      alpha = design$p3_alpha,
      multiplicity = design$p3_multiplicity
    )

    # STEP 2B: Program 2 runs Phase 3 only if Phase 2 receives Go.
    if (phase2$go) {

      selected <- simulate_phase3(
        truth,
        doses = c(0, phase2$selected_dose),
        n_per_arm = design$p2_program_p3_n,
        sigma = design$sigma,
        alpha = design$p3_alpha,
        multiplicity = design$p3_multiplicity
      )

      selected_records[[program]] <- dplyr::mutate(selected, program = program)
    }

    # Save the decisions and observations so individual programs can be inspected.
    phase2_records[[program]] <- tibble::tibble(
      program = program,
      posterior_draw = posterior_rows[program],
      go = phase2$go,
      selected_dose = phase2$selected_dose,
      significant = phase2$significant,
      t_max = phase2$t_max
    )

    phase2_arms[[program]] <- dplyr::mutate(phase2$arms, program = program)
    direct_records[[program]] <- dplyr::mutate(direct, program = program)

    if (program %% 5000L == 0) message("Simulated ", program, " programs")
  }

  phase2 <- dplyr::bind_rows(phase2_records)
  phase3_direct <- dplyr::bind_rows(direct_records)
  phase3_selected <- dplyr::bind_rows(selected_records)

  # Give an empty selected-trial table the correct columns when every program stops.
  if (nrow(phase3_selected) == 0) phase3_selected <- phase3_direct[0, ]

  direct_events <- program_events(
    phase3_direct, phase2, design$thresholds, "P1_direct"
  )

  # Apply the gate to the SAME direct Phase 3 outcomes. This isolates screening
  # from changes in doses or sample sizes and avoids resimulating Program 1.
  screened_events <- direct_events |>
    dplyr::mutate(
      strategy = "P1_screened",
      dplyr::across(c(EM, CM, both, both_same_dose), ~ .x & go)
    )

  selected_events <- program_events(
    phase3_selected, phase2, design$thresholds, "P2_selected"
  )

  list(
    events = dplyr::bind_rows(direct_events, screened_events, selected_events),
    phase2 = phase2,
    phase2_arms = dplyr::bind_rows(phase2_arms),
    phase3_direct = phase3_direct,
    phase3_selected = phase3_selected,
    truth = draws[posterior_rows, , drop = FALSE] |>
      dplyr::mutate(program = seq_len(n_sim), posterior_draw = posterior_rows),
    design = design,
    mcp = mcp,
    seed = seed
  )
}

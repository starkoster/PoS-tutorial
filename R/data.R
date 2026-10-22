read_analysis_arms <- function(
    path = system.file("extdata", "supplied_arms.csv", package = "MigraineEmaxSim"),
    correct_hope_inputs = TRUE) {

  # Keep the supplied CSV unchanged; corrections are explicit derived columns.
  arms <- readr::read_csv(path, show_col_types = FALSE)

  required <- c(
    "nct", "title", "CM", "arm", "dose", "n",
    "base_mmd", "base_mmd_sd", "chg", "chg_se"
  )

  stopifnot(
    all(required %in% names(arms)),
    all(arms$n > 0),
    all(arms$chg_se > 0),
    all(arms$CM >= 0 & arms$CM <= 1)
  )

  arms <- arms |>
    dplyr::mutate(
      input_cm = CM,
      input_baseline_sd = base_mmd_sd
    )

  if (correct_hope_inputs) {

    # HOPE Table 1 gives arm-specific CM counts and placebo baseline SD 4.6.
    hope_counts <- tibble::tibble(
      dose = c(0, 100, 750),
      cm_count = c(62, 32, 68)
    )

    unknown_doses <- arms |>
      dplyr::filter(title == "HOPE", !dose %in% hope_counts$dose)

    if (nrow(unknown_doses) > 0) {
      stop("Supply the actual CM count for each additional HOPE dose.")
    }

    arms <- arms |>
      dplyr::left_join(hope_counts, by = "dose") |>
      dplyr::mutate(
        CM = dplyr::if_else(title == "HOPE", cm_count / n, CM),
        base_mmd_sd = dplyr::if_else(
          title == "HOPE" & dose == 0,
          4.6,
          base_mmd_sd
        )
      ) |>
      dplyr::select(-cm_count)
  }

  # LY baselines were reported for pooled arms, not separately for EM and CM.
  # Exclude their repeated baseline entries; keep every change observation.
  arms |>
    dplyr::mutate(
      w_em = 1 - CM,
      baseline_se = base_mmd_sd / sqrt(n),
      use_baseline = title != "LY3451838" &
        is.finite(baseline_se) & baseline_se > 0
    )
}

analysis_theme <- function() {

  # Keep the white background and restrained colors of the presentation.
  ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", colour = "#001965"),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(colour = "#001965")
    )
}


save_analysis_plot <- function(plot, name, directory, width = 10, height = 6) {

  # Save an editable vector figure and a convenient image preview.
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    file.path(directory, paste0(name, ".pdf")),
    plot, width = width, height = height
  )

  ggplot2::ggsave(
    file.path(directory, paste0(name, ".png")),
    plot, width = width, height = height, dpi = 180
  )

  invisible(plot)
}


plot_distributions <- function(
    draw_sets,
    quantities = default_targets()$quantity,
    title = "Clinical parameter distributions") {

  # Stack named sets of posterior/prior draws into one tidy plotting table.
  data <- purrr::imap_dfr(draw_sets, function(stage_draws, stage_name) {

    stage_draws |>
      dplyr::select(dplyr::all_of(quantities)) |>
      # Normalize each quantity before stacking stages or pivoting them.
      as_numeric_draws() |>
      dplyr::mutate(stage = stage_name)
  }) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(quantities),
      names_to = "quantity",
      values_to = "value"
    ) |>
    dplyr::mutate(stage = factor(stage, levels = names(draw_sets)))

  ggplot2::ggplot(data, ggplot2::aes(value, colour = stage)) +
    ggplot2::geom_density(linewidth = 0.7) +
    ggplot2::facet_wrap(~ quantity, scales = "free", ncol = 3) +
    ggplot2::scale_colour_manual(values = c("#8E99AB", "#001965", "#C66A2E")) +
    ggplot2::labs(x = "MMDs", y = "Density", colour = NULL, title = title) +
    analysis_theme()
}


plot_calibration <- function(calibration) {

  # Show the achieved distribution alongside all three elicited target quantiles.
  data <- calibration$draws |>
    dplyr::select(dplyr::all_of(calibration$targets$quantity)) |>
    as_numeric_draws() |>
    tidyr::pivot_longer(dplyr::everything(), names_to = "quantity", values_to = "value")

  targets <- calibration$targets |>
    tidyr::pivot_longer(c(q05, q50, q95), names_to = "quantile", values_to = "target")

  ggplot2::ggplot(data, ggplot2::aes(value)) +
    ggplot2::geom_density(colour = "#001965", fill = "#E8EEF9") +
    ggplot2::geom_vline(
      data = targets,
      ggplot2::aes(xintercept = target, linetype = quantile),
      colour = "#C66A2E"
    ) +
    ggplot2::facet_wrap(~ quantity, scales = "free", ncol = 3) +
    ggplot2::labs(x = "MMDs", y = "Density", title = "Clinical prior calibration") +
    analysis_theme()
}


plot_em_cm <- function(draw_sets) {

  # Use matched EM/CM draws. Thinning is for plotting speed, not for inference.
  data <- purrr::imap_dfr(draw_sets, function(stage_draws, stage_name) {

    keep <- unique(round(seq(
      1, nrow(stage_draws), length.out = min(3000, nrow(stage_draws))
    )))

    purrr::map_dfr(c("baseline", "placebo", "active"), function(quantity_name) {

      # Extract both vectors before building the tibble. Its sequential column
      # evaluation must not replace the scalar name used to select a draw column.
      em_values <- as.numeric(stage_draws[[paste0(quantity_name, "_em")]])[keep]
      cm_values <- as.numeric(stage_draws[[paste0(quantity_name, "_cm")]])[keep]

      tibble::tibble(
        stage = stage_name,
        quantity = quantity_name,
        EM = em_values,
        CM = cm_values
      )
    })
  }) |>
    dplyr::mutate(stage = factor(stage, levels = names(draw_sets)))

  ggplot2::ggplot(data, ggplot2::aes(EM, CM)) +
    ggplot2::geom_point(alpha = 0.08, size = 0.6, colour = "#001965") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#8E99AB") +
    ggplot2::facet_wrap(~ stage + quantity, scales = "free", ncol = 3) +
    ggplot2::labs(x = "EM (MMDs)", y = "CM (MMDs)", title = "EM/CM clinical relationships") +
    analysis_theme()
}


plot_meta_intervals <- function(prior_draws, fit, historical) {

  # Combine initial prediction, individual studies and new-study prediction.
  studies <- unique(historical$title)

  study_summaries <- purrr::map_dfr(seq_along(studies), function(study_index) {

    eta <- draw_matrix(fit, paste0("eta_study[", study_index, ",", 1:6, "]"))

    clinical_quantiles(clinical_from_eta(eta)) |>
      dplyr::mutate(source = studies[study_index])
  })

  prediction <- predictive_prior_from_meta(fit)

  data <- dplyr::bind_rows(
    dplyr::mutate(clinical_quantiles(prior_draws), source = "Initial new-study prior"),
    study_summaries,
    dplyr::mutate(
      clinical_quantiles(prediction$clinical_draws),
      source = "Posterior new-study prediction"
    )
  ) |>
    dplyr::filter(quantity %in% default_targets()$quantity) |>
    dplyr::mutate(source = factor(
      source,
      levels = rev(c("Initial new-study prior", studies, "Posterior new-study prediction"))
    ))

  ggplot2::ggplot(data, ggplot2::aes(q50, source)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = q05, xend = q95, yend = source), colour = "#001965"
    ) +
    ggplot2::geom_point(colour = "#C66A2E", size = 1.8) +
    ggplot2::facet_wrap(~ quantity, scales = "free_x", ncol = 3) +
    ggplot2::labs(
      x = "Median and central 90% interval (MMDs)", y = NULL,
      title = "Historical trials, priors and posterior predictions"
    ) +
    analysis_theme()
}


plot_emax_curves <- function(draws, doses = seq(0, 750, 5), adjusted = TRUE) {

  # Calculate pointwise posterior intervals, retaining joint posterior rows.
  data <- purrr::map_dfr(c("EM", "CM"), function(population_name) {

    purrr::map_dfr(doses, function(dose_value) {

      effect <- mean_change(draws, dose_value, population_name)

      if (adjusted) effect <- effect - mean_change(draws, 0, population_name)

      interval <- stats::quantile(effect, c(0.05, 0.50, 0.95), names = FALSE)

      tibble::tibble(
        population = population_name,
        dose = dose_value,
        lower = interval[1],
        median = interval[2],
        upper = interval[3]
      )
    })
  })

  ggplot2::ggplot(data, ggplot2::aes(dose, median)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "#B8C7E6") +
    ggplot2::geom_line(colour = "#001965", linewidth = 0.8) +
    ggplot2::facet_wrap(~ population) +
    ggplot2::labs(
      x = "Dose (mg)",
      y = if (adjusted) "Active minus placebo (MMDs)" else "Change from baseline (MMDs)",
      title = "Updated Emax dose response"
    ) +
    analysis_theme()
}


plot_ed50_update <- function(draws, prior) {

  # The prior density is analytic, so plotting it does not consume random draws.
  grid <- seq(0.01, 100, length.out = 500)

  prior_density <- tibble::tibble(
    ed50 = grid,
    density = stats::dlnorm(grid, prior$logmean, prior$logsd)
  )

  ggplot2::ggplot(draws, ggplot2::aes(ed50)) +
    ggplot2::geom_density(ggplot2::aes(colour = "HOPE posterior")) +
    ggplot2::geom_line(
      data = prior_density,
      ggplot2::aes(ed50, density, colour = "Separate lognormal prior")
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 100)) +
    ggplot2::scale_colour_manual(values = c(
      "HOPE posterior" = "#001965", "Separate lognormal prior" = "#8E99AB"
    )) +
    ggplot2::labs(x = "ED50 (mg)", y = "Density", colour = NULL) +
    analysis_theme()
}


plot_programs <- function(summary) {

  # Compare complete strategies using the same from-start denominator.
  data <- summary$probabilities |>
    dplyr::filter(estimand == "from_start", population %in% c("EM", "CM", "both")) |>
    dplyr::mutate(threshold = factor(
      threshold, levels = c(0, -2, -3, -4),
      labels = c("Statistical", "Delta <= -2", "Delta <= -3", "Delta <= -4")
    ))

  ggplot2::ggplot(data, ggplot2::aes(threshold, probability, fill = strategy)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~ population) +
    ggplot2::scale_fill_manual(values = c("#001965", "#6F8DCC", "#C66A2E")) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "Success from program start", fill = NULL) +
    analysis_theme()
}


plot_derisking <- function(summary) {

  # Separate conditional improvement after Go from success at program start.
  data <- summary$derisking |>
    dplyr::filter(population %in% c("EM", "CM", "both")) |>
    tidyr::pivot_longer(
      c(direct, given_go, gated_from_start),
      names_to = "estimand", values_to = "probability"
    ) |>
    dplyr::mutate(threshold = factor(threshold, levels = c(0, -2, -3, -4)))

  ggplot2::ggplot(data, ggplot2::aes(threshold, probability, colour = estimand, group = estimand)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ population) +
    ggplot2::scale_colour_manual(values = c(
      "direct" = "#8E99AB", "given_go" = "#001965", "gated_from_start" = "#C66A2E"
    )) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Observed-effect threshold", y = "Probability", colour = NULL) +
    analysis_theme()
}


plot_selection <- function(summary) {

  # Selection frequencies are conditional on receiving Go in Phase 2.
  ggplot2::ggplot(summary$selection, ggplot2::aes(factor(selected_dose), fraction_of_go)) +
    ggplot2::geom_col(fill = "#001965") +
    ggplot2::labs(x = "Selected dose (mg)", y = "Fraction among Go programs") +
    analysis_theme()
}

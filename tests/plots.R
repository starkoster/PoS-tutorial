library(MigraineEmaxSim)

# This test isolates plotting from calibration, Stan fitting and trial simulation.
# More than 3,000 rows reproduces the data size in the reported masking failure.
set.seed(901)
draws <- clinical_from_eta(matrix(rnorm(6002 * 6), ncol = 6))

for (n_rows in c(1L, 2L, 3001L, 6002L)) {

  subset <- draws[seq_len(n_rows), ]
  keep <- unique(round(seq(1, n_rows, length.out = min(3000, n_rows))))

  for (stage_names in list("Prior", c("Prior", "Historical"))) {

    draw_sets <- stats::setNames(rep(list(subset), length(stage_names)), stage_names)
    plot <- plot_em_cm(draw_sets)

    # Each stage and clinical quantity must retain the paired EM/CM observations.
    stopifnot(nrow(plot$data) == length(stage_names) * 3L * length(keep))

    for (stage_name in stage_names) {
      for (quantity_name in c("baseline", "placebo", "active")) {

        data <- plot$data[
          as.character(plot$data$stage) == stage_name &
            plot$data$quantity == quantity_name,
        ]

        stopifnot(
          isTRUE(all.equal(data$EM, subset[[paste0(quantity_name, "_em")]][keep])),
          isTRUE(all.equal(data$CM, subset[[paste0(quantity_name, "_cm")]][keep]))
        )
      }
    }

    # Build the plot so facet and geometry problems also surface before saving.
    invisible(ggplot2::ggplot_build(plot))
  }
}

# Exercise the sibling plotting callbacks without any fitted Stan object.
draw_sets <- list(Prior = draws[1:100, ], Historical = draws[101:200, ])
distribution_plot <- plot_distributions(draw_sets)

stopifnot(setequal(as.character(distribution_plot$data$stage), names(draw_sets)))
invisible(ggplot2::ggplot_build(distribution_plot))

posterior <- dplyr::mutate(draws[1:100, ], ed50 = 20)
curve_plot <- plot_emax_curves(posterior, doses = c(0, 20, 200))

stopifnot(
  nrow(curve_plot$data) == 6,
  all(curve_plot$data$median[curve_plot$data$dose == 0] == 0)
)

invisible(ggplot2::ggplot_build(curve_plot))

message("Plotting regression checks passed.")

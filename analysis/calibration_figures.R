# Redraw the two prior calibration figures from saved validation draws.
# Run this script from the directory containing results_R_readable/.
# Change output_dir below if your results folder has another name.
#
# Install these plotting packages once if needed:
# install.packages(c("tidyverse", "patchwork"))

library(tidyverse)
library(patchwork)


# STEP 1: Read the completed calibration. Plotting does not rerun the search.

output_dir <- "results_R_readable"

figure_dir <- file.path(output_dir, "figures")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

calibration <- readRDS(file.path(output_dir, "01_calibration.rds"))


# OPTIONAL: To run a new calibration with the 50/50 variance allocation,
# uncomment the following block before continuing with the plots.
# This replaces the calibration object and saves the new results.
# The original analysis script used different fractions in some versions,
# so the 50/50 choice is specified explicitly here.
#
# library(MigraineEmaxSim)
#
# targets <- default_targets()
#
# calibration <- calibrate_priors(
#   targets = targets,
#   heterogeneity_fraction = rep(0.5, 6),
#   lkj_shape = 2,
#   n_train = 20000L,
#   n_validate = 100000L,
#   seed = 4210,
#   maxit = 3500L,
#   n_starts = 3L
# )
#
# saveRDS(calibration, file.path(output_dir, "01_calibration.rds"))
#
# write_csv(
#   calibration$comparison,
#   file.path(output_dir, "prior_calibration.csv")
# )


# STEP 2: Keep ordinary numeric columns, preserving the pairing of rows.
# Older saved outputs can contain one-column draws_matrix objects.
# Converting each scalar column avoids pivoting and CSV class errors.

required_columns <- c(
  "baseline_em", "baseline_cm",
  "placebo_em", "placebo_cm",
  "emax_em", "emax_cm",
  "active_em", "active_cm"
)

raw_draws <- calibration$draws

if (!all(required_columns %in% names(raw_draws))) {
  stop("The calibration draws must contain all eight clinical columns.")
}

draws <- map_dfc(required_columns, function(column_name) {

  values <- raw_draws[[column_name]]

  # A matrix column is accepted only when it has one value per row.
  if (!is.numeric(values) || length(values) != nrow(raw_draws)) {
    stop("Expected one numeric value per draw in column: ", column_name)
  }

  tibble(!!column_name := as.numeric(values))
})

targets <- as_tibble(calibration$targets)

# Show the variance allocation recorded in this saved run.
# Reading an older RDS does not change its allocation to 50/50.
print(calibration$heterogeneity_fraction)


# STEP 3: Define labels and a generic plotting style.

quantity_order <- c(
  "baseline_em", "placebo_em", "emax_em",
  "baseline_cm", "placebo_cm", "emax_cm"
)

panel_labels <- c(
  baseline_em = "Baseline: EM",
  placebo_em = "Placebo change: EM",
  emax_em = "Active change: EM",
  baseline_cm = "Baseline: CM",
  placebo_cm = "Placebo change: CM",
  emax_cm = "Active change: CM"
)

prior_color <- "tomato2"
text_color <- "black"

calibration_theme <- function() {

  theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(
        colour = "#E8EBEE",
        linewidth = 0.3
      ),
      strip.text = element_text(colour = text_color, size = 10),
      axis.text = element_text(colour = text_color),
      axis.title = element_text(colour = text_color),
      panel.spacing = grid::unit(0.8, "cm"),
      plot.margin = margin(8, 8, 8, 8)
    )
}


# STEP 4: Calculate the density curves for the six calibrated quantities.
# Every density uses all saved draws. Only the displayed horizontal range
# is restricted, so long tails do not obscure the main part of the curve.
# The range is extended where needed to include the elicited targets.

density_data <- map_dfr(quantity_order, function(quantity_name) {

  values <- draws[[quantity_name]]

  target_values <- targets |>
    filter(.data$quantity == .env$quantity_name) |>
    select(q05, q50, q95) |>
    unlist(use.names = FALSE)

  display_limits <- range(
    quantile(values, probs = c(0.0025, 0.9975), names = FALSE),
    target_values
  )

  # Smooth the complete sample, evaluating the curve on this display grid.
  curve <- density(
    values,
    from = display_limits[1],
    to = display_limits[2],
    n = 512
  )

  tibble(
    quantity = quantity_name,
    value = curve$x,
    density = curve$y
  )
}) |>
  mutate(quantity = factor(quantity, levels = quantity_order))


target_lines <- targets |>
  select(quantity, q05, q50, q95) |>
  pivot_longer(
    cols = c(q05, q50, q95),
    names_to = "quantile",
    values_to = "target"
  ) |>
  mutate(quantity = factor(quantity, levels = quantity_order))


# STEP 5: Draw marginal densities and the three elicited target lines.

prior_calibration_plot <- ggplot(
  density_data,
  aes(x = value, y = density)
) +
  geom_area(fill = prior_color, alpha = 0.12) +
  geom_line(colour = prior_color, linewidth = 0.8) +
  geom_vline(
    data = target_lines,
    aes(xintercept = target),
    colour = text_color,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  facet_wrap(
    vars(quantity),
    ncol = 3,
    scales = "free",
    labeller = as_labeller(panel_labels)
  ) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "MMDs", y = "Density") +
  calibration_theme()

print(prior_calibration_plot)


# STEP 6: Select matched rows for EM/CM comparisons.
# The same row numbers are used in every panel. Selection is for display
# only and does not change the saved calibration or its quantile summaries.

keep_rows <- unique(round(seq(
  from = 1,
  to = nrow(draws),
  length.out = min(1800L, nrow(draws))
)))

pair_titles <- c(
  baseline = "Baseline",
  placebo = "Placebo change",
  active = "Active change"
)


# Each panel has its own limits, but equal units on its EM and CM axes.
# Separate plots avoid the restriction on fixed aspect ratios with free facets.

pair_plots <- map(names(pair_titles), function(quantity_name) {

  # Extract vectors before creating a table, avoiding name masking.
  em_values <- draws[[paste0(quantity_name, "_em")]]
  cm_values <- draws[[paste0(quantity_name, "_cm")]]

  pair_data <- tibble(
    EM = em_values[keep_rows],
    CM = cm_values[keep_rows]
  )

  # Calculate viewing limits from the full sample, not the selected points.
  axis_limits <- quantile(
    c(em_values, cm_values),
    probs = c(0.005, 0.995),
    names = FALSE
  )

  padding <- 0.04 * diff(axis_limits)
  axis_limits <- axis_limits + c(-padding, padding)

  ggplot(pair_data, aes(x = EM, y = CM)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      colour = "#AAB2BA",
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_point(
      colour = prior_color,
      alpha = 0.22,
      size = 0.7
    ) +
    coord_fixed(
      ratio = 1,
      xlim = axis_limits,
      ylim = axis_limits,
      expand = FALSE
    ) +
    scale_x_continuous(n.breaks = 4) +
    scale_y_continuous(n.breaks = 4) +
    labs(
      title = unname(pair_titles[[quantity_name]]),
      x = "EM (MMDs)",
      y = "CM (MMDs)"
    ) +
    calibration_theme() +
    theme(plot.title = element_text(size = 10, hjust = 0.5))
})

prior_em_cm_plot <- wrap_plots(pair_plots, nrow = 1)

print(prior_em_cm_plot)


# STEP 7: Save vector PDFs for LaTeX and PNG previews.
# Output names match the figure references in the manuscript.

save_calibration_plot <- function(plot, name, width, height) {

  ggsave(
    filename = file.path(figure_dir, paste0(name, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )

  ggsave(
    filename = file.path(figure_dir, paste0(name, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

save_calibration_plot(
  prior_calibration_plot,
  name = "01_prior_calibration",
  width = 7.0,
  height = 4.3
)

save_calibration_plot(
  prior_em_cm_plot,
  name = "01_prior_em_cm",
  width = 7.2,
  height = 2.8
)

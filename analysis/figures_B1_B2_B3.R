# Redraw Figures B1, B2 and B3 from the completed historical analysis.
# No Stan model is refitted.
#
# B1: Prior, study-specific posteriors and new-study predictive intervals.
# B2: Prior versus original historical new-study predictive densities.
# B3: Joint EM/CM comparisons before and after the historical update.
#
# Run from the project directory containing results_R_readable/.
# Requires MigraineEmaxSim >= 0.5.3, tidyverse, posterior, cmdstanr and patchwork.
# Install any missing plotting packages once, outside this script.

library(MigraineEmaxSim)
library(tidyverse)
library(posterior)
library(cmdstanr)
library(patchwork)


# STEP 1: Load the saved calibration and historical fit.
# The package uses fit$save_object(), which saves the draws in this RDS.

output_dir <- "results_R_readable"
figure_dir <- file.path(output_dir, "figures")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

calibration <- readRDS(file.path(output_dir, "01_calibration.rds"))
meta <- readRDS(file.path(output_dir, "stan_meta", "fit.rds"))

# Alternative if the saved fit is unavailable: read one specific run's CSVs.
# Edit the timestamp to identify your completed historical run.
# Never combine CSV files from different runs.
#
# chain_files <- list.files(
#   file.path(output_dir, "stan_meta"),
#   pattern = "^historical_meta-202610161019-[1-4]-.*\\.csv$",
#   full.names = TRUE
# )
# stopifnot(length(chain_files) == 4L)
# meta <- cmdstanr::as_cmdstan_fit(chain_files)

arms <- read_csv(
  file.path(output_dir, "analysis_arms.csv"),
  show_col_types = FALSE
)

# This reproduces the ordering used by historical_stan_data().
# The study index follows first appearance of each NCT identifier.
studies <- arms |>
  filter(.data$title != "HOPE") |>
  distinct(nct, .keep_all = TRUE) |>
  select(nct, title)

print(studies)
print(calibration$heterogeneity_fraction)


# STEP 2: Extract ordinary numeric draws in the requested column order.
# Rebuilding a base matrix prevents draws_matrix classes from leaking into
# clinical columns, which caused earlier binding and plotting errors.

extract_eta <- function(variable_names) {

  x <- posterior::as_draws_matrix(meta$draws(variables = variable_names))

  x <- matrix(
    as.numeric(x),
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = dimnames(x)
  )

  x[, variable_names, drop = FALSE]
}

prior_draws <- as_numeric_draws(calibration$draws)

eta_new <- extract_eta(paste0("eta_new[", 1:6, "]"))

# Use the original predictive draws, not the Gaussian transfer approximation.
historical_draws <- clinical_from_eta(eta_new) |>
  as_numeric_draws()


# STEP 3: Set labels and generic colors shared by all three figures.

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

stage_levels <- c("Initial prior", "Historical prediction")

stage_colors <- c(
  "Initial prior" = "tomato2",
  "Historical prediction" = "steelblue"
)

figure_theme <- function() {

  theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E8EBEE", linewidth = 0.3),
      strip.text = element_text(size = 10, colour = "#333B43"),
      axis.text = element_text(colour = "#333B43"),
      axis.title = element_text(colour = "#333B43"),
      legend.title = element_blank(),
      legend.position = "bottom",
      panel.spacing = grid::unit(0.8, "cm"),
      plot.margin = margin(8, 8, 8, 8)
    )
}


# FIGURE B1: Median and central 90% interval for every source.
# Study intervals describe clinical mean parameters under the fitted model.
# They are not raw arm estimates or intervals for individual patient outcomes.

study_summaries <- map_dfr(seq_len(nrow(studies)), function(study_index) {

  variable_names <- paste0(
    "eta_study[", study_index, ",", 1:6, "]"
  )

  eta <- extract_eta(variable_names)

  clinical_quantiles(clinical_from_eta(eta)) |>
    mutate(source = studies$title[[study_index]])
})

interval_data <- bind_rows(
  clinical_quantiles(prior_draws) |>
    mutate(source = "Initial prior"),

  study_summaries,

  clinical_quantiles(historical_draws) |>
    mutate(source = "New study")
) |>
  filter(.data$quantity %in% quantity_order) |>
  mutate(
    quantity = factor(quantity, levels = quantity_order),
    source = factor(
      source,
      levels = rev(c("Initial prior", studies$title, "New study"))
    )
  )

# Lighten steelblue toward white for the individual historical trials.
# The new-study prediction retains the full steelblue color.
lighten_color <- function(color, white_fraction) {

  original <- grDevices::col2rgb(color)[, 1]

  lighter <- original + white_fraction * (255 - original)

  grDevices::rgb(
    red = lighter[1],
    green = lighter[2],
    blue = lighter[3],
    maxColorValue = 255
  )
}

trial_colors <- map_chr(
  seq(0.25, 0.45, length.out = nrow(studies)),
  function(white_fraction) {
    lighten_color("steelblue", white_fraction)
  }
)

source_colors <- c(
  "Initial prior" = stage_colors[["Initial prior"]],
  setNames(trial_colors, studies$title),
  "New study" = stage_colors[["Historical prediction"]]
)

B1 <- ggplot(interval_data, aes(x = q50, y = source, colour = source)) +
  geom_segment(
    aes(x = q05, xend = q95, yend = source),
    linewidth = 0.8,
    lineend = "round"
  ) +
  geom_point(size = 2) +
  facet_wrap(
    vars(quantity),
    ncol = 3,
    scales = "free_x",
    labeller = as_labeller(panel_labels)
  ) +
  scale_colour_manual(values = source_colors) +
  scale_x_continuous(n.breaks = 4) +
  labs(x = "Median and central 90% interval (MMDs)", y = NULL) +
  figure_theme() +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 8)
  )


# FIGURE B2: Compare the prior and original historical predictive densities.
# All draws contribute to density estimation. Only the displayed range is
# restricted to keep long tails from obscuring the main distributions.

draw_sets <- list(
  "Initial prior" = prior_draws,
  "Historical prediction" = historical_draws
)

density_data <- map_dfr(quantity_order, function(quantity_name) {

  combined_values <- c(
    prior_draws[[quantity_name]],
    historical_draws[[quantity_name]]
  )

  display_limits <- quantile(
    combined_values,
    probs = c(0.0025, 0.9975),
    names = FALSE
  )

  imap_dfr(draw_sets, function(stage_draws, stage_name) {

    curve <- density(
      stage_draws[[quantity_name]],
      from = display_limits[1],
      to = display_limits[2],
      n = 512
    )

    tibble(
      quantity = quantity_name,
      stage = stage_name,
      value = curve$x,
      density = curve$y
    )
  })
}) |>
  mutate(
    quantity = factor(quantity, levels = quantity_order),
    stage = factor(stage, levels = stage_levels)
  )

B2 <- ggplot(
  density_data,
  aes(x = value, y = density, colour = stage, fill = stage)
) +
  geom_area(position = "identity", alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(
    vars(quantity),
    ncol = 3,
    scales = "free",
    labeller = as_labeller(panel_labels)
  ) +
  scale_colour_manual(values = stage_colors) +
  scale_fill_manual(values = stage_colors) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "MMDs", y = "Density") +
  figure_theme() +
  theme(panel.grid.major.x = element_blank())


# FIGURE B3: Compare matched EM/CM draws at each stage.
# Each point uses EM and CM values from the same draw row.
# Stages are separate distributions, with no pairing between prior and posterior.
# Use the same 1,800 selected rows across the three panels within each stage.

plot_rows <- map(draw_sets, function(stage_draws) {

  unique(round(seq(
    from = 1,
    to = nrow(stage_draws),
    length.out = min(1800L, nrow(stage_draws))
  )))
})

pair_titles <- c(
  baseline = "Baseline",
  placebo = "Placebo change",
  active = "Active change"
)

pair_plots <- map(names(pair_titles), function(quantity_name) {

  pair_data <- imap_dfr(draw_sets, function(stage_draws, stage_name) {

    # Extract vectors before building the tibble to avoid column-name masking.
    em_values <- stage_draws[[paste0(quantity_name, "_em")]]
    cm_values <- stage_draws[[paste0(quantity_name, "_cm")]]
    keep <- plot_rows[[stage_name]]

    tibble(
      stage = stage_name,
      EM = em_values[keep],
      CM = cm_values[keep]
    )
  }) |>
    mutate(stage = factor(stage, levels = stage_levels))

  # Viewing limits use the full marginal draws pooled across both stages.
  all_values <- unlist(map(draw_sets, function(stage_draws) {
    c(
      stage_draws[[paste0(quantity_name, "_em")]],
      stage_draws[[paste0(quantity_name, "_cm")]]
    )
  }), use.names = FALSE)

  limits <- quantile(all_values, c(0.005, 0.995), names = FALSE)
  padding <- 0.04 * diff(limits)
  limits <- limits + c(-padding, padding)

  ggplot(pair_data, aes(x = EM, y = CM, colour = stage)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      colour = "#AAB2BA",
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_point(alpha = 0.22, size = 0.7) +
    scale_colour_manual(values = stage_colors, drop = FALSE) +
    coord_fixed(ratio = 1, xlim = limits, ylim = limits, expand = FALSE) +
    scale_x_continuous(n.breaks = 4) +
    scale_y_continuous(n.breaks = 4) +
    labs(
      title = unname(pair_titles[[quantity_name]]),
      x = "EM (MMDs)",
      y = "CM (MMDs)"
    ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    figure_theme() +
    theme(plot.title = element_text(size = 10, hjust = 0.5))
})

B3 <- wrap_plots(pair_plots, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")


# STEP 4: Display the figures and save PDF and PNG versions.
# The filenames match the existing manuscript figure references.

print(B1)
print(B2)
print(B3)

save_figure <- function(plot, name, width, height) {

  ggsave(
    file.path(figure_dir, paste0(name, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )

  ggsave(
    file.path(figure_dir, paste0(name, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

save_figure(B1, "02_historical_trials", width = 7.2, height = 4.8)
save_figure(B2, "02_meta_update", width = 7.0, height = 4.3)
save_figure(B3, "02_meta_em_cm", width = 7.2, height = 3.0)

# Numerical summaries still use every draw. The viewing limits and selected
# scatter points affect only the figures, not the stored analysis or its PoS.

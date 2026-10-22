# Redraw Figures C1 and C2 using tomato2 and steelblue.
# No historical or HOPE model is refitted.
#
# C1: Original historical prediction versus Gaussian transfer.
# C2: Separate lognormal ED50 prior versus HOPE posterior.
#
# Run from the project directory containing results_R_readable/.
# Requires MigraineEmaxSim >= 0.5.3 and tidyverse.

library(MigraineEmaxSim)
library(tidyverse)


# STEP 1: Read the saved inputs from the completed analysis.

output_dir <- "results_R_readable"
figure_dir <- file.path(output_dir, "figures")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

predictive <- readRDS(file.path(output_dir, "02_predictive_prior.rds"))
ed50_prior <- readRDS(file.path(output_dir, "03_ed50_prior.rds"))

posterior <- read_csv(
  file.path(output_dir, "emax_posterior.csv"),
  show_col_types = FALSE
)

# Normalize older one-column matrix containers without changing draw pairing.
historical_draws <- as_numeric_draws(predictive$clinical_draws)


# STEP 2: Generate draws from the actual Gaussian prior used for HOPE.
# Use the stored mean and covariance, including any recorded covariance
# multiplier and diagonal regularization. Do not add either adjustment again.
# This is only a descriptive transfer comparison, not a new model fit.

set.seed(20261016)
n_transfer <- 500000L

standard_normal <- matrix(rnorm(n_transfer * 6), ncol = 6)

# R's chol() returns an upper-triangular factor. Multiplication on the right
# gives row-vector draws with the intended covariance matrix.
gaussian_eta <- standard_normal %*% chol(predictive$cov)

gaussian_eta <- sweep(
  gaussian_eta,
  MARGIN = 2,
  STATS = predictive$mean,
  FUN = "+"
)

gaussian_draws <- clinical_from_eta(gaussian_eta) |>
  as_numeric_draws()

# These R random draws need not equal the earlier Python plotting draws.
# They target the same saved Gaussian distribution.


# STEP 3: Choose panel labels and colors.
# "Active change" is the requested shorter panel title. By default it still
# refers to the placebo-adjusted Emax quantity plotted in the original C1.
# Set effect_type to "total" to plot A = placebo + Emax instead.
# That setting changes the plotted quantity and should be stated in the caption.

effect_type <- "placebo_adjusted"

stopifnot(effect_type %in% c("placebo_adjusted", "total"))

effect_prefix <- if (effect_type == "total") "active" else "emax"

quantity_order <- c(
  "baseline_em", "placebo_em", paste0(effect_prefix, "_em"),
  "baseline_cm", "placebo_cm", paste0(effect_prefix, "_cm")
)

panel_labels <- setNames(
  c(
    "Baseline: EM", "Placebo change: EM", "Active change: EM",
    "Baseline: CM", "Placebo change: CM", "Active change: CM"
  ),
  quantity_order
)

transfer_colors <- c(
  "Historical prediction" = "steelblue",
  "Gaussian transfer" = "tomato2"
)

potency_colors <- c(
  "Separate prior" = "tomato2",
  "After HOPE" = "steelblue"
)

figure_theme <- function() {

  theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(
        colour = "#E8EBEE",
        linewidth = 0.3
      ),
      strip.text = element_text(size = 10, colour = "#333B43"),
      axis.text = element_text(colour = "#333B43"),
      axis.title = element_text(colour = "#333B43"),
      legend.title = element_blank(),
      legend.position = "bottom",
      panel.spacing = grid::unit(0.8, "cm"),
      plot.margin = margin(8, 8, 8, 8)
    )
}


# FIGURE C1: Original prediction and Gaussian transfer on the clinical scale.
# All draws contribute to each density. The displayed range covers the pooled
# 0.25th to 99.75th percentiles so that extreme tails do not dominate the axes.

draw_sets <- list(
  "Historical prediction" = historical_draws,
  "Gaussian transfer" = gaussian_draws
)

density_data <- map_dfr(quantity_order, function(quantity_name) {

  combined_values <- c(
    historical_draws[[quantity_name]],
    gaussian_draws[[quantity_name]]
  )

  limits <- quantile(
    combined_values,
    probs = c(0.0025, 0.9975),
    names = FALSE
  )

  imap_dfr(draw_sets, function(stage_draws, stage_name) {

    curve <- density(
      stage_draws[[quantity_name]],
      from = limits[1],
      to = limits[2],
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
    stage = factor(stage, levels = names(transfer_colors))
  )

C1 <- ggplot(
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
  scale_colour_manual(values = transfer_colors) +
  scale_fill_manual(values = transfer_colors) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "MMDs", y = "Density") +
  figure_theme()


# FIGURE C2: Compare the separate ED50 prior with the saved HOPE posterior.
# The prior density is evaluated analytically from its saved log parameters.
# Smooth the posterior on log dose, then transform back to the mg scale.
# This prevents the density curve from extending to negative potency values.

ed50 <- posterior$ed50

stopifnot(all(is.finite(ed50)), all(ed50 > 0))

ed50_grid <- seq(0.01, 100, length.out = 1000)

prior_density <- tibble(
  stage = "Separate prior",
  ed50 = ed50_grid,
  density = dlnorm(
    ed50_grid,
    meanlog = ed50_prior$logmean,
    sdlog = ed50_prior$logsd
  )
)

log_density <- density(
  log(ed50),
  from = log(min(ed50_grid)),
  to = log(max(ed50_grid)),
  n = 2048
)

# If Z = log(ED50), then f_ED50(x) = f_Z(log(x)) / x.
posterior_density <- tibble(
  stage = "After HOPE",
  ed50 = ed50_grid,
  density = approx(
    x = log_density$x,
    y = log_density$y,
    xout = log(ed50_grid)
  )$y / ed50_grid
)

potency_data <- bind_rows(prior_density, posterior_density) |>
  mutate(stage = factor(stage, levels = names(potency_colors)))

C2 <- ggplot(
  potency_data,
  aes(x = ed50, y = density, colour = stage, fill = stage)
) +
  geom_area(position = "identity", alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = potency_colors) +
  scale_fill_manual(values = potency_colors) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, NA)) +
  labs(x = expression(ED[50]~"(mg)"), y = "Density") +
  figure_theme()


# STEP 4: Display the figures and save PDF and PNG versions.
# Filenames match the manuscript's existing LaTeX figure references.

print(C1)
print(C2)

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

save_figure(C1, "transfer_comparison", width = 7.0, height = 4.3)
save_figure(C2, "03_ed50_update", width = 5.0, height = 3.0)

# Suggested C1 caption with the default setting:
# Original historical predictive distributions and their Gaussian approximation.
# Active change panels show the maximum placebo-adjusted effect, Emax.
# Curves are descriptive density estimates on the clinical scale.
#
# Suggested C2 caption:
# Separately elicited lognormal ED50 prior and the HOPE-updated posterior.
# The posterior density is smoothed on log dose and transformed back to mg.

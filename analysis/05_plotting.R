library(tidyverse)


# Change this path if your analysis outputs are stored elsewhere.
output_dir <- "results_R_readable"

figure_dir <- file.path(output_dir, "figures")

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Read the saved HOPE-updated posterior.
# Each row contains one joint posterior draw.
posterior <- read_csv(
  file.path(output_dir, "emax_posterior.csv"),
  show_col_types = FALSE
)


# Calculate the dose-response curve for one population.
# Negative effects indicate reductions in MMDs relative to placebo.
summarise_emax_curve <- function(draws, population) {

  doses <- seq(0, 300, by = 2)

  emax <- draws[[paste0("emax_", population)]]

  ed50 <- draws$ed50


  # At each dose, calculate the effect for every posterior draw.
  # Matching Emax and ED50 values remain paired within each row.
  map_dfr(doses, function(dose) {

    effect <- emax * dose / (dose + ed50)


    # Summarise the posterior effect distribution at this dose.
    interval <- quantile(
      effect,
      probs = c(0.05, 0.50, 0.95),
      names = FALSE
    )

    tibble(
      population = population,
      dose = dose,
      lower = interval[1],
      median = interval[2],
      upper = interval[3]
    )
  })
}


# Combine the population-specific summaries for plotting.
curve_summary <- bind_rows(
  summarise_emax_curve(posterior, "em"),
  summarise_emax_curve(posterior, "cm")
) |>
  mutate(
    population = factor(
      population,
      levels = c("em", "cm"),
      labels = c("Episodic migraine", "Chronic migraine")
    )
  )


# Use generic colors for the two populations.
population_colors <- c(
  "Episodic migraine" = "steelblue",
  "Chronic migraine" = "tomato2"
)


# Plot posterior medians and pointwise central 90% intervals.
figure_2 <- ggplot(
  curve_summary,
  aes(x = dose)
) +

  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper,
      fill = population
    ),
    alpha = 0.17,
    colour = NA
  ) +

  geom_hline(
    yintercept = -2,
    colour = "#929CA5",
    linetype = "dashed",
    linewidth = 0.4
  ) +

  geom_line(
    aes(
      y = median,
      colour = population
    ),
    linewidth = 0.9
  ) +

  facet_wrap(
    vars(population),
    nrow = 1
  ) +

  scale_colour_manual(values = population_colors) +

  scale_fill_manual(values = population_colors) +

  scale_x_continuous(
    breaks = c(0, 100, 200, 300),
    expand = expansion(mult = c(0, 0.02))
  ) +

  scale_y_continuous(
    breaks = seq(-5, 0, by = 1)
  ) +

  coord_cartesian(
    xlim = c(0, 300),
    ylim = c(-5.3, 0.15)
  ) +

  labs(
    x = "Dose (mg)",
    y = "Active - placebo (in MMDs)"
  ) +

  theme_minimal(base_size = 10) +

  theme(
    legend.position = "none",

    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),

    panel.grid.major.y = element_line(
      colour = "#E8EBEE",
      linewidth = 0.3
    ),

    strip.text = element_text(
      size = 11,
      colour = "#333B43"
    ),

    axis.title = element_text(colour = "#333B43"),
    axis.text = element_text(colour = "#333B43"),

    panel.spacing = grid::unit(1, "cm"),

    plot.margin = margin(8, 10, 8, 8)
  )


# Display the figure.
print(figure_2)


# Save a vector PDF for LaTeX and a PNG for previewing.
ggsave(
  filename = file.path(figure_dir, "figure_2_emax.pdf"),
  plot = figure_2,
  width = 6.8,
  height = 2.7,
  units = "in",
  bg = "white"
)

ggsave(
  filename = file.path(figure_dir, "figure_2_emax.png"),
  plot = figure_2,
  width = 6.8,
  height = 2.7,
  units = "in",
  dpi = 300,
  bg = "white"
)



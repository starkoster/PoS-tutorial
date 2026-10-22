\name{plot_calibration}
\alias{plot_calibration}
\alias{analysis_theme}
\alias{save_analysis_plot}
\alias{plot_distributions}
\alias{plot_em_cm}
\alias{plot_meta_intervals}
\alias{plot_emax_curves}
\alias{plot_ed50_update}
\alias{plot_programs}
\alias{plot_derisking}
\alias{plot_selection}
\title{Simple figures for the four analysis stages}
\description{Create clinical prior/posterior plots, historical trial intervals,
Emax curves and program comparisons using a restrained blue/orange style.}
\usage{
analysis_theme()
save_analysis_plot(plot, name, directory, width = 10, height = 6)
plot_calibration(calibration)
plot_distributions(draw_sets, quantities = default_targets()$quantity,
                   title = "Clinical parameter distributions")
plot_em_cm(draw_sets)
plot_meta_intervals(prior_draws, fit, historical)
plot_emax_curves(draws, doses = seq(0, 750, 5), adjusted = TRUE)
plot_ed50_update(draws, prior)
plot_programs(summary)
plot_derisking(summary)
plot_selection(summary)
}
\arguments{
\item{plot}{A ggplot object.}
\item{name}{Output filename stem.}
\item{directory}{Figure directory.}
\item{width, height}{Figure dimensions in inches.}
\item{calibration}{Calibration result.}
\item{draw_sets}{Named list of clinical draw tables.}
\item{quantities}{Clinical columns to plot.}
\item{title}{Plot title.}
\item{prior_draws}{Initial new-study clinical draws.}
\item{fit}{Historical CmdStanMCMC fit.}
\item{historical}{Historical subset of the prepared arm table.}
\item{draws}{Joint clinical Emax posterior rows.}
\item{doses}{Dose grid.}
\item{adjusted}{Plot active minus placebo rather than total change from baseline.}
\item{prior}{ED50 prior specification.}
\item{summary}{Output of summarize_programs().}
}
\value{A ggplot object or theme. The saving helper writes PDF and PNG and
invisibly returns the plot.}
\details{Clinical intervals are central 90 percent intervals. Curve bands are
pointwise. Historical subgroup intervals are model-derived when summaries were
mixed; they are not independent observed EM/CM estimates.}

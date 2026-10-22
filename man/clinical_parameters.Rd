\name{clinical_from_eta}
\alias{clinical_from_eta}
\alias{as_numeric_draws}
\alias{eta_from_reference}
\alias{clinical_quantiles}
\alias{mean_change}
\title{Clinical parameters and Emax mean changes}
\description{Convert latent parameters to constrained clinical quantities and
evaluate the Emax mean. Baselines are positive; changes are negative for benefit.}
\usage{
clinical_from_eta(eta)
as_numeric_draws(draws)
eta_from_reference(reference)
clinical_quantiles(draws)
mean_change(draw, dose, population = c("EM", "CM"))
}
\arguments{
\item{eta}{A length-six numeric vector or a matrix with six columns.}
\item{reference}{Named vector or list containing baseline_em, baseline_cm,
placebo_em, placebo_cm, emax_em and emax_cm at a coherent joint reference point.}
\item{draws}{A data frame with one numeric clinical quantity per column.
Numeric one-column matrices from older saved outputs are also accepted.}
\item{draw}{One joint posterior row, or a table of posterior rows.}
\item{dose}{Dose or vector of doses. With multiple posterior rows, use a scalar dose.}
\item{population}{EM or CM.}
}
\value{The transformation returns a tibble of eight clinical quantities; its
inverse returns six latent values. Quantile summaries return a tidy table.
The mean function returns changes from baseline. The numeric conversion returns
a tibble containing ordinary double vectors, with values and row order retained.}
\details{The mapping enforces ordered baselines below 30, beneficial drug effects,
and stronger underlying CM effects. These restrictions apply to mean parameters,
not to normally distributed observations. Emax is a placebo-adjusted effect;
the total active change equals placebo plus Emax at the asymptote.
The numeric conversion rejects list columns and matrices with more than one
column; each quantity must supply exactly one number per posterior row.}

\name{summarize_programs}
\alias{summarize_programs}
\alias{probability_summary}
\title{Summarize program success and fixed-design screening}
\description{Report probabilities with explicit denominators, conditional
screening gains, selected doses and expected participant counts.}
\usage{
probability_summary(success)
summarize_programs(simulation)
}
\arguments{
\item{success}{Logical event indicators.}
\item{simulation}{Output of simulate_programs().}
}
\value{A probability row, or a list containing probability, screening,
Phase 2 Go, selection and participant tables.}
\details{MCSE measures simulation noise conditional on the available posterior.
If no programs advance, conditional probabilities are NA. Joint probabilities
are calculated from joint event indicators, not marginal products.}

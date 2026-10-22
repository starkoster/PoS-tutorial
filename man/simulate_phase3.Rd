\name{simulate_phase3}
\alias{simulate_phase3}
\alias{phase3_success}
\title{Separate EM and CM Phase 3 trials}
\description{Simulate normal arm means with shared placebo comparisons and
evaluate beneficial statistical and observed-effect success.}
\usage{
simulate_phase3(draw, doses, n_per_arm, sigma = 5, alpha = 0.05,
                multiplicity = c("none", "bonferroni"))
phase3_success(phase3, threshold = -2)
}
\arguments{
\item{draw}{One joint posterior row.}
\item{doses}{Distinct nonnegative doses including placebo and active treatment.}
\item{n_per_arm}{Participants per arm in each population-specific trial.}
\item{sigma}{Fixed patient SD.}
\item{alpha}{Two-sided significance level.}
\item{multiplicity}{Nominal comparisons or Bonferroni across active doses within a trial.}
\item{phase3}{Phase 3 arm-comparison table.}
\item{threshold}{Observed negative TPP threshold; zero gives statistical success alone.}
}
\value{A tibble of comparisons, or one row of EM, CM, joint and common-dose
success indicators.}
\details{EM and CM sampling errors are independent conditional on the shared
truth. Active doses share one placebo mean within a trial. A qualifying result
must be beneficial and significant in addition to meeting the observed TPP.
Joint success can occur at different doses unless common-dose success is used.}

\name{read_analysis_arms}
\alias{read_analysis_arms}
\title{Read and prepare supplied trial arms}
\description{Read arm summaries, retain original input fields and apply the
documented preparation for mixed populations and baseline information.}
\usage{
read_analysis_arms(
  path = system.file("extdata", "supplied_arms.csv", package = "MigraineEmaxSim"),
  correct_hope_inputs = TRUE
)
}
\arguments{
\item{path}{CSV path with the supplied arm-summary columns.}
\item{correct_hope_inputs}{Use published HOPE arm-specific CM counts and placebo
baseline SD 4.6. FALSE deliberately retains the supplied HOPE values.}
}
\value{A tibble with original input fields, EM fractions, baseline SEs and
an indicator specifying whether to use the baseline observation.}
\details{LY pooled-arm baselines are excluded from the subgroup likelihood;
their change summaries remain included. Raw HOPE composition and baseline
SD are retained in input_cm and input_baseline_sd. See the extdata README.}

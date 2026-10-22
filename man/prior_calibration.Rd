\name{calibrate_priors}
\alias{calibrate_priors}
\alias{default_targets}
\alias{validate_targets}
\alias{simulate_prior}
\title{Calibrate clinical quantiles to latent hyperpriors}
\description{Find latent centers and scales whose new-study prediction
approximately reproduces elicited clinical quantiles.}
\usage{
default_targets()
validate_targets(targets)
calibrate_priors(
  targets = default_targets(),
  heterogeneity_fraction = c(4 / 13, 4 / 13, rep(0.36, 4)),
  lkj_shape = 2, n_train = 20000L, n_validate = 100000L,
  seed = 4210, maxit = 3500L, n_starts = 3L
)
simulate_prior(prior, n = 50000L, seed = 4202)
}
\arguments{
\item{targets}{Table with quantity, q05, q50 and q95 for the six default quantities.}
\item{heterogeneity_fraction}{Six fractions between zero and one, allocating
total latent variance to study heterogeneity rather than mean uncertainty.}
\item{lkj_shape}{Positive LKJ shape parameter.}
\item{n_train, n_validate}{Prior simulation sizes for optimization and independent validation.}
\item{seed}{Random seed.}
\item{maxit}{Maximum Nelder--Mead iterations per start.}
\item{n_starts}{Number of starting points.}
\item{prior}{List with latent mean centers m, mean SDs s, half-normal scale
parameters q, and lkj_shape.}
\item{n}{Number of new-study draws.}
}
\value{Calibration returns prior settings, target and achieved quantiles,
clinical draws and the optimizer result. Simulation returns clinical draws.}
\details{The targets concern underlying new-study mean parameters, not
individual outcomes. Common random numbers are reused during optimization;
independent draws validate the result. The heterogeneity allocation is a
separate assumption. Calibration is approximate and requires clinical review.}

\name{fit_meta}
\alias{fit_meta}
\alias{historical_stan_data}
\alias{predictive_prior_from_meta}
\title{Historical meta-analysis and new-study prediction}
\description{Fit a constrained hierarchical Stan model and transfer its joint
new-study predictive distribution to the drug-specific analysis.}
\usage{
historical_stan_data(arms, prior, active_weight = 0.5)
fit_meta(arms, prior, output_dir = "results_R/stan_meta",
         active_weight = 0.5, seed = 4301, ...)
predictive_prior_from_meta(fit, covariance_multiplier = 1)
}
\arguments{
\item{arms}{Prepared arm table from read_analysis_arms().}
\item{prior}{Calibrated latent hyperpriors.}
\item{active_weight}{Power weight between zero and one for active change observations.}
\item{output_dir}{Writable directory for compiled model, Stan CSVs and fit.rds.}
\item{seed}{Stan random seed.}
\item{...}{Sampler settings: chains (4), warmup (4000), sampling (4000),
adapt_delta (0.995) and max_treedepth (13).}
\item{fit}{Historical CmdStanMCMC fit.}
\item{covariance_multiplier}{Positive multiplier of the empirical latent
predictive covariance; the default is one.}
}
\value{Stan data, a CmdStanMCMC fit, or a predictive-prior list with mean,
covariance and clinical predictive draws, respectively.}
\details{HOPE is excluded. The hierarchy includes half-normal heterogeneity
SDs and an LKJ correlation prior. The new-study prediction includes
heterogeneity. Its transfer is a single Gaussian moment approximation;
a diagonal ridge of 1e-8 is added. Native sampler warnings remain enabled.}

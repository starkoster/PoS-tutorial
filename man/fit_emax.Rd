\name{fit_emax}
\alias{fit_emax}
\alias{make_ed50_prior}
\alias{hope_stan_data}
\alias{emax_draws}
\title{Drug-specific Emax update with a separate potency prior}
\description{Use HOPE to update the historical predictive prior and an
independently elicited lognormal ED50 prior in Stan.}
\usage{
make_ed50_prior(q05 = sqrt(40), q50 = 20, q95 = sqrt(4000))
hope_stan_data(arms, predictive_prior, ed50_prior)
fit_emax(arms, predictive_prior, ed50_prior = make_ed50_prior(),
         output_dir = "results_R/stan_hope", seed = 4401, ...)
emax_draws(fit)
}
\arguments{
\item{q05, q50, q95}{Coherent positive lognormal quantiles: q50 must equal
the geometric mean of q05 and q95.}
\item{arms}{Prepared arm table.}
\item{predictive_prior}{Joint latent predictive-prior list from the meta-analysis.}
\item{ed50_prior}{Lognormal prior list returned by make_ed50_prior().}
\item{output_dir}{Writable directory for the model, Stan CSVs and saved fit.}
\item{seed}{Stan seed.}
\item{...}{Sampler settings as for fit_meta().}
\item{fit}{HOPE CmdStanMCMC fit.}
}
\value{A prior specification, Stan data, a CmdStanMCMC fit or a tibble of
joint clinical posterior draws with ED50. All nine draw columns are ordinary
numeric vectors, suitable for CSV export and plotting.}
\details{ED50 is shared between EM and CM. Arm-specific population weights
combine their mean curves. Whole posterior rows retain potency--effect dependence.}

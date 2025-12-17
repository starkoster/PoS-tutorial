#' Simulate trial with continuous outcome
#'
#' @param n            Sample size, provided as a scalar or vector
#' @param delta_prior  Prior mean, either two scalars (mean, std.dev, for normal distributed) or 'rdist' type function
#' @param sigma_prior  Prior std.dev, either scalar (for known variance) or 'rdist' type function
#' @param arms         Number of trial arms
#' @param nsim         Number of trial simulations
#'
#' @return 2-column matrix with simulated true and observed effect sizes
#'
#' @examples
#' # define priors for effect size and sigma
#' dprior = function(x) rnorm(x, 2,2)      # for normal prior, this could also be specified as 'dprior = c(2,2)'
#' sprior = function(x) rgamma(x, 6.5,1)   # if variance is known, then 'sprior = 6.5' is also valid
#' sim_continuous(n=100, delta_prior=dprior, sigma_prior=sprior)


sim_continuous <- function(n, delta_prior, sigma_prior, arms = 2, nsim = 1e3){

  if(is.function(delta_prior)){ # input delta_prior is a distribution
    dtru = delta_prior(nsim)
  } else{
    dtru = delta_prior[1]+delta_prior[2]*rnorm(nsim, 0, 1)
  }

  if(is.function(sigma_prior)){ # input sigma_prior is a distribution, assume equal variances in all arms
    sigma = sigma_prior(nsim)
    tau = sqrt(arms*sigma^2/n)
  } else{
    if(length(sigma_prior>1)){ # unequal variance, and/or sample size in each arm
      tau = sqrt(sum(sigma_prior^2/n))
    } else if(length(n)>1){ # unequal sample size in each arm
      tau = sqrt(sum(sigma_prior^2/n))
    } else{ # assume arms with equal variance and sample size
      tau = sqrt(arms*sigma_prior^2/n)
    }
  }

  dobs = dtru+tau*rnorm(nsim, 0, 1)

  out = cbind(dtru, dobs)
  colnames(out) = c("True", "Obs")

  return(out)
}




#' Calculate assurance
#'
#' @param n            Sample size, provided as a scalar or vector
#' @param sigma        Observed std.dev.
#' @param sigma_prior  Prior std.dev
#' @param mu_prior     Prior mean
#' @param alpha        Significance threshold
#'
#' @return scalar or vector (depending on n) of assurance probabilites
#'
#' @examples
#' n = 1:1000
#' assurance(n, sigma=6.5, sigma_prior = 2, mu_prior=2)
#'


assurance <- function(n,sigma,sigma_prior,mu_prior=0,alpha=0.025){

  if(length(n)>1){
    out = unlist(lapply(n, assurance, sigma=sigma, sigma_prior=sigma_prior, mu_prior=mu_prior, alpha=alpha))
  } else{
    tau <- sigma*sqrt(2/n)
    out = pnorm((qnorm(alpha)*tau+mu_prior)/sqrt(sigma_prior^2+tau^2))
  }

  return(out)
}


#' Calculate power for continuous endpoint
#'
#' @param n         Sample size, provided as a scalar or vector
#' @param delta     Mean effect
#' @param sigma     Std.dev. of mean effect
#' @param alpha     Significance threshold
#' @param beta      Power
#'
#' @return scalar or vector (depending on n) of power probabilites
#'
#' @examples
#' n = 1:1000
#' pwr(n, delta = 2, sigma=6.5)
#'
#'

pwr <- function(n,delta,sigma, alpha=0.05, beta=0.80){               #function that returns power

  if(length(n)>1){
    out = unlist(lapply(n, pwr, delta=delta, sigma=sigma, alpha=alpha, beta=beta))
  } else{
    out = pnorm(qnorm(alpha)+(sqrt(n)*delta)/(sqrt(2)*sigma))
  }

  return(out)
}





#' Find Gamma parameters from input mean/variance, or find mean/variance from input Gamma parameters
#'
#' @param EX      Mean of variable
#' @param VX      Variance of variable
#' @param shape   Gamma shape parameter
#' @param rate    Gamma rate parameter (=1/scale)
#'
#' @return either mean/variance or shape/rate parameters
#'
#' @examples
#'
#' repar_gamma(EX = 6.5, VX = .5)
#'
#' repar_gamma(shape = 84.5, rate = 13)
#'
repar_gamma <- function(EX=NULL, VX=NULL, shape=NULL, rate=NULL){
  if(is.null(EX) & is.null(VX)){
    EX = shape/rate
    VX = shape/rate^2
    return(c(EX=EX, VX=VX))

  } else if(is.null(shape) & is.null(rate)){
    shape = EX^2/VX
    rate = EX/VX
    return(c(shape=shape, rate=rate))
  }
}


#' Calculate std.dev. from a confidence interval
#'
#' @param ci        Confidence interval limits, given as vector of length 2
#' @param n         Sample size
#' @param ci_level  Significance level of confidence interval
#'
#' @return standard deviation, scaled with sample size
#'
#' @examples
#' # 95% confidence interval of [0,6], from 100 subjects:
#' ci2sd(c(0,6), n=100)
#' # standard devation is 15.3064
#'

ci2sd <- function(ci, n=1, ci_level = 0.95){
  za = qnorm(0.5*(1+ci_level))
  out = sqrt(n)*abs(diff(ci))/(2*za)

  return(out)
}

#' Approximate standard error from data range
#'
#' @param lo    Lowest observation
#' @param hi    highest observation
#' @param smax  max standard devation to look for
#' @param eps   Precision value for
#'
#' @return standard deviation, scaled with sample size
#'
#' @examples
#' # 95% confidence interval of [0,6], from 100 subjects:
#' ci2sd(c(0,6), n=100)
#' # standard devation is 15.3064
#'
range2sd <- function(lo, up, smax = 10, eps = 1e-4){

  out = optimize(\(sd) (abs(up-lo)/2-qnorm(1-eps, sd=sd))^2, c(eps,smax))$minimum
  return(out)
}

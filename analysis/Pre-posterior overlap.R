preposterior_overlap <- function(n){

    sim_pwr <- function(d, n){
      sig_prior = 6.5
      d_prior = c(d,0) # true effect fixed at d, second argument => prior std.dev.  = 0

      out = sim_continuous(n=n, delta_prior=d_prior, sigma_prior=sig_prior, arms = 2, nsim = 1e3)
      tau = sqrt(2*sig_prior^2/n)
      Za  = qnorm(1-0.025)

      out = mean(out[,2] > tau*Za)
      return(out)
    }

    set.seed(1234)
    ddt = 0.01 # precision for integration
    ds = seq(-2,8, ddt) # evaluation range
    pw = map_dbl(ds, sim_pwr, n = n) # simulate pwr for varying true effect sizes
    # Pre posteriors:
    y1 = dnorm(ds,2,2)*pw/(sum(dnorm(ds,2,2)*pw*ddt)) # P(true delta | success)
    y2 = dnorm(ds,2,2)*(1-pw)/(sum(dnorm(ds,2,2)*(1-pw)*ddt)) # P(true delta | failure)

    # Overlap probability
    idx1 = which(y1<y2)
    idx2 = which(y2<=y1)
    y3 = c(y1[idx1],y2[idx2])

    out = sum(y3)*ddt
    return(out)
}

ns  = c(10,25,50,100, 150, 200, 400, 600, 800, 1000)
res = map_dbl(ns, preposterior_overlap)

plot(ns, res, type='l')
spout = spline(ns, res, xout = seq(1,1000))
lines(spout$x, spout$y, col='tomato2')

diff(spout$y) %>% plot(type='l')

f1 <- function(x, y){

  n = length(y)

  out = numeric(n)
  for(i in 1:(n-1)){
    out[i] = y[i+1]+y[i]
  }

  out = cbind(x=x, fx =out)
  out = data.frame(out)

  return(out)
}

f2 <- function(x, y){

  n = length(y)

  out = numeric(n)
  for(i in 2:(n-1)){
    out[i] = y[i+1]+y[i-1]-2*y[i]
  }

  out = cbind(x=x, fx =out)
  out = data.frame(out)
  return(out)
}

f1(spout$x, spout$y) %>% plot
f2(spout$x, spout$y) %>% plot

f1d = f1(spout$x, spout$y)$fx
f2d = f2(spout$x, spout$y)$fx

abs(f2d)/(1+abs(f1d)^2)^(3/2) %>% plot



library(PoSTutorial)
library(tidyverse)
clean_up()

# Specify priors ----
  ## Effect prior ----
  dprior = function(x) rnorm(x, 2, 2) # can also be specified as dprior = c(2,2), i.e. c(mu, sigma)

  ## Observed std.dev. prior ----
  gprior = 6.5 # either fixed (known variance)
  gprior = function(x) rgamma(x, 6.5, rate=1) # or use prior
  # curve(dgamma(x, 6.5, rate=1), from=0, to = 15)

# Simulate trials ----

  # Simulate continuous endpoint with unknown sigma (using prior distribution)
  set.seed(1234)
  sim1 = sim_continuous(n=100, delta_prior=dprior, sigma_prior=gprior, arms = 2, nsim = 5e3) %>% data.frame

  # Simulate continuous endpoint with known sigma
  set.seed(1234)
  sim2 = sim_continuous(n=100, delta_prior=dprior, sigma_prior=6.5, arms = 2, nsim = 5e3) %>% data.frame

## Plot and compare simulations ----
  # Compare simulated trials
  cat("Correlation between known/unknown variance:",round(cor(sim1$Obs, sim2$Obs),3))
  par(mfrow=c(1,1), mar=c(4,4,2,0), oma=c(0,0,0,0))
  plot(sim1$Obs, sim2$Obs, pch=16, cex=.5,
       col=add_alpha('tomato2',.5), bty='n',
       xlab=expression("Unknown"~sigma~"~"~Gamma~"(6.5,1)"),
       ylab=expression("Known"~sigma),
       main=expression("Observed effects for known/unknown"~sigma))
  abline(0,1, lty=3)

  # Compared true vs observed effects
  par(mfrow=c(1,1), mar=c(4,4,0,0), oma=c(0,0,0,0))
  plot(NA, bty='n', xlim=c(-5,10), ylim=c(-5,10), xlab="True effect",ylab="Observed effect")
  sim1 %>% points(pch=16, col=add_alpha('tomato2',.25), cex=.5)
  sim2 %>% points(pch=16, col=add_alpha('steelblue',.25), cex=.5)
  legend("top", c("Known variance", "Unknown variance"), fill=c("tomato2", "steelblue"), bty='n')

# Evaluate trial PoS ----
  tau = sqrt(2*6.5^2/100)
  Za  = qnorm(1-0.025)
  lrv = 2
  tv  = 2.5

  obs = sim1$Obs
  cat("Probability of significance (better than ",round(tau*Za,3),"): ",   round(mean(obs>tau*Za),3),
      "\nProbability of lrv success (better than ",lrv,"):      ", round(mean(obs>max(tau*Za, lrv)),3),
      "\nProbability of lrv success (better than ",tv,"):    ", round(mean(obs>max(tau*Za, tv )),3),
      sep="")

  cat("P(significance | true effect>0):", mean(sim1$True>0 & sim1$Obs>tau*Za),
      "\nP(significance | true effect<0):", mean(sim1$True<0 & sim1$Obs>tau*Za))

  cat("P(significance | true effect >lrv):", mean(sim1$True>lrv & sim1$Obs>tau*Za),
      "\nP(significance | true effect >tv): ", mean(sim1$True>tv & sim1$Obs>tau*Za))


# Simulate Pre-posteriors ----
  sim_pwr <- function(d){
    n   = 100
    sig_prior = 6.5
    d_prior = c(d,0) # true effect fixed at d, second argument => prior std.dev.  = 0

    out = sim_continuous(n=n, delta_prior=d_prior, sigma_prior=sig_prior, arms = 2, nsim = 5e3)
    tau = sqrt(2*sig_prior^2/n)
    Za  = qnorm(1-0.025)

    out = mean(out[,2] > tau*Za)
    return(out)
  }

  set.seed(1234)
  ddt = 0.01 # precision for integration
  ds = seq(-2,8, ddt) # evaluation range
  pw = map_dbl(ds, sim_pwr) # simulate pwr for varying true effect sizes
  # Pre posteriors:
    y1 = dnorm(ds,2,2)*pw/(sum(dnorm(ds,2,2)*pw*ddt)) # P(true delta | success)
    y2 = dnorm(ds,2,2)*(1-pw)/(sum(dnorm(ds,2,2)*(1-pw)*ddt)) # P(true delta | failure)

  # Plotting pre-posteriors
  par(mfrow=c(1,1), mar=c(4,4,0,0), oma=c(0,0,0,0))
  plot(NA, bty='n', xlim=c(-4,14), ylim=c(0, max(y1,y2)), xlab="True effect",ylab="Observed effect")
  lines(ds, y1, lwd=2, col=add_alpha('tomato2',.75))
  lines(ds, y2, lwd=2, col=add_alpha('steelblue',.75))
  legend("topright", c("P(d | succes)", "P(d | failure)"), fill=c("tomato2","steelblue"), bty='n')


  # Overlap probability
  idx1 = which(y1<y2)
  idx2 = which(y2<=y1)
  y3 = c(y1[idx1],y2[idx2])
  polygon(c(ds,rev(ds)), c(y3, rep(0,length(ds))), border=NA, col=add_alpha('magenta',.25))

  cat("Pre-posterior overlap is", round(sum(y3)*ddt,3))







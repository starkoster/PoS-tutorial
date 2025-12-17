# clean_up()
library(tidyverse)
library(readr)
library(readxl)

res <- read_excel("data/Migraine trial results.xlsx") # Load data from spreadsheet

res = res %>%
  mutate(nActive  = ifelse(is.na(nActive), pts/2, nActive),    # use half of study population if #pts in each arm is not available
         nPlacebo = ifelse(is.na(nPlacebo), pts/2, nPlacebo),  # use half of study population if #pts in each arm is not available
         nEst     = nActive+nPlacebo,                          # number of pts for est
         se = (hi-lo)/2*qnorm(1-0.05/2),
         sd = sqrt(nEst)*se) %>%
  mutate(id2 = 1:n()) #unique identifier, if using trial number, several studies have multiple rows



cor.test(res$baseline, res$est)
par(mar=c(4,4,2,2), oma = c(0,0,0,0))
plot(res$baseline, res$est, cex=(res$nEst/200), pch=16, col='tomato2', bty='n', xlim=c(5,20), ylim=c(-4,2),
     xlab="MMD at baseline",
     ylab=expression(Delta~"MMD"),
     main="MMD at baseline vs treatment effect")
abline(h=0, lty=1, col='gray60')
legend("topright", paste0("# pts=", c(250,500,750)), pch=16, col='tomato2', pt.cex=c(250,500,750)/200, y.intersp = 2, x.intersp = 1.5, bty='n', xpd=TRUE)



# glimpse(res)
# res %>% distinct(id2, est) %>% knitr::kable()

# Set priors for meta analysis
mu0 = 0
sd0 = 2 # qualitatively guessed expected effect
f <- function(x, prob = 0.975, target=8) (qcauchy(prob, 0, x)-target)^2 # find Cauchy scale parameter such that Prob(x < 8) = 0.95 (specify here as 0.975 due to half-cauchy prior!), where 8 is 4x "expected effect"
tau0 = optimize(f, c(0,3))$minimum
pcauchy(c(2,4,8), 0, tau0)


library(brms)
M0 = brm(data = res, family = gaussian,
         est | se(se) ~ 1 + (1 | id2),
         prior = c(prior(normal(0, 2), class = Intercept),  # input mu0 and sd0
                   prior(cauchy(0, 0.630), class = sd)),    # input tau0
         iter = 2000, warmup = 1000, cores = 4, chains = 4,
         seed = 1234)
print(M0)
# tidybayes::summarise_draws(M0)


# Forest plot of meta analysis
par(mar=c(2,0,2,0), oma = c(0,0,0,0))
cols = c("steelblue", "magenta","tomato2")
ns = nrow(res)
plot(NA, bty='n', xlim=c(-5,2), ylim=c(-1,ns), ylab="", yaxt="n", xlab="Effect size")
abline(v=0, lty=1, col='gray60')
mtext("Mean effects and 95% CI's", line=.5)

# Add observed estimates and CIs
points(res$est, 1:ns, pch=16, col=cols[1])
segments(x0 = res$lo, x1 = res$hi, y0=1:ns, y1=1:ns, col=cols[1])
text(rep(-5, ns), 1:ns, labels = res$name, cex=0.5, adj = 0, xpd=TRUE)
text(rep(-4.5, ns), 1:ns, labels = paste0(res$dose,"mg"), cex=0.5, adj = 0, xpd=TRUE)
text(rep(2, ns+1), 1:(ns+1), labels = c(res$nEst,"n"), cex=0.5)
abline(h=which(res$name != lag(res$name))-.5, lty=3, col='gray75')


# Add meta analysis outcome, ie prior estimate
tmp = fixef(M0)
points(tmp[1],0,  pch=16, col=cols[2])
segments(x0 = tmp[3], x1 = tmp[4], y0=0, y1=0, lwd=2, col=cols[2])
text(rep(2, ns+1), 1:(ns+1), labels = c(res$nEst,"n"), cex=0.5, )
text(-5, 0, labels = substitute(paste(bold('PIII prior'))), cex=0.75, adj = 0, xpd=TRUE)

abline(h=.5, lty=3, col='gray75')

# Prior predictive probabilities
round(100*pnorm(c(0,-1,-2,-3), mean=tmp[1], sd = tmp[2]),2)
round(qnorm(c(.25,.5,.75), mean=tmp[1], sd=tmp[2]),2)

# Add prior predictive outcome: see https://opensource.nibr.com/bamdd/src/02ab_meta_analysis_trtdiff.html
tmp = predict(object=M0, newdata=tibble(id="New study", se=10e-19), allow_new_levels=TRUE, sample_new_levels = "gaussian")
points(tmp[1],-1,  pch=16, col=cols[3])
segments(x0 = tmp[3], x1 = tmp[4], y0=-1, y1=-1, lwd=2, col=cols[3])
text(-5, -1, labels = substitute(paste(bold('PIII predictive'))), cex=0.75, adj =0, xpd=TRUE)

# Prior predictive probabilities
round(100*pnorm(c(0,-1,-2,-3), mean=tmp[1], sd = tmp[2]),2)
round(qnorm(c(.25,.5,.75), mean=tmp[1], sd=tmp[2]),2)



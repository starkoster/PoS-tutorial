# clean_up()
library(PoSTutorial)
load("data/mmd.Rda")


# Correlation between baseline MMD and outcome ----
par(mar=c(4,4,2,2), oma = c(0,0,0,0))
mmd %>%
  select(md_mean, est) %>%
  plot(cex=log(mmd$ntotal/250), pch=16, col=add_alpha('tomato2',.75), bty='n', xlim=c(5,20), ylim=c(-4,2),
     xlab="MMD at baseline",
     ylab=expression(Delta~"MMD"),
     main="MMD at baseline vs treatment effect")
abline(h=0, lty=1, col='gray60')
legend("topright", paste0("# pts=", c(250,500,750)), pch=16, col='tomato2', pt.cex=log(c(250,500,750)/200), y.intersp = 2, x.intersp = 1.5, bty='n', xpd=TRUE)
cor.test(mmd$md_mean, mmd$est)
M0  = lm(est~md_mean, data=mmd, weights = ntotal)
x   = seq(5,20,1)
tmp = predict(M0, newdata = data.frame(md_mean = x),se.fit = TRUE)
tmp$hi = tmp$fit+1.96*tmp$se.fit
tmp$lo = tmp$fit-1.96*tmp$se.fit
polygon(c(x,rev(x)), c(tmp$hi, rev(tmp$lo)), border=NA, col=add_alpha('black',.2))
coef(M0) %>% abline(col=add_alpha('black',.25))
summary(M0)


# Meta analysis ----
glimpse(mmd)

# Set priors for meta analysis
mu0 = 0
sd0 = 2 # qualitatively guessed expected effect
f <- function(x, prob = 0.975, target=8) (qcauchy(prob, 0, x)-target)^2 # find Cauchy scale parameter such that Prob(x < 8) = 0.95 (specify here as 0.975 due to half-cauchy prior!), where 8 is 4x "expected effect"
tau0 = optimize(f, c(0,3))$minimum
# pcauchy(c(2,4,8), 0, tau0)

mmd = mmd %>%
  mutate(se = abs(lower-upper)/2*qnorm(1-0.05/2),
         id = 1:n()) %>%
  relocate(id)

library(brms)
M1 = brm(data = mmd, family = gaussian,
         est | se(se) ~ 1 + (1 | id),
         prior = c(prior(normal(0, 2), class = Intercept),  # input mu0 and sd0
                   prior(cauchy(0, 0.630), class = sd)),    # input tau0
         iter = 2000, warmup = 1000, cores = 4, chains = 4,
         seed = 1234)
print(M1)
# tidybayes::summarise_draws(M1)

remove_all <- function(string, x){
  for(s in x){
    string = str_replace_all(string, s,"")
  }
  return(string)
}


plt = mmd %>%
  mutate(lbls = tolower(mmd$group) %>% remove_all(c("placebo"," vs ", "double-blind phase: ", "period: ", "treatment ","/.", '[[:digit:]]+ ', " mg", "mg "))) %>%
  arrange(lbls, dose)

# Forest plot of meta analysis
par(mar=c(2,0,2,0), oma = c(0,0,0,0))
cols = c("steelblue", "magenta","tomato2")
ns = nrow(plt)

plot(NA, bty='n', xlim=c(-5,2), ylim=c(-1,ns), ylab="", yaxt="n", xlab="Effect size")
abline(v=0, lty=1, col='gray60')
mtext("Mean effects and 95% CI's", line=.5)

# Add observed estimates and CIs
points(plt$est, ns:1, pch=16, col=cols[1])
segments(x0 = plt$lower, x1 = plt$upper, y0=ns:1, y1=ns:1, col=cols[1])


text(rep(-4.75, ns), ns:1, labels = plt$lbls, cex=0.5, adj = 0, xpd=TRUE)
text(rep(-5, ns), ns:1, labels = paste0(plt$dose,"mg"), cex=0.5, adj = 0, xpd=TRUE)

text(rep(2, ns+1), (ns+1):1, labels = c("n", plt$n), cex=0.5)
abline(h=2+ns-which(plt$nct != lag(plt$nct))-.5, lty=3, col='gray75')


# Add meta analysis outcome, ie prior estimate
tmp = fixef(M1)
points(tmp[1],0,  pch=16, col=cols[2])
segments(x0 = tmp[3], x1 = tmp[4], y0=0, y1=0, lwd=2, col=cols[2])
text(-5, 0, labels = substitute(paste(bold('PIII prior'))), cex=0.75, adj = 0, xpd=TRUE)

abline(h=.5, lty=3, col='gray75')

# Prior predictive probabilities
round(100*pnorm(c(0,-1,-2,-3), mean=tmp[1], sd = tmp[2]),2)
round(qnorm(c(.25,.5,.75), mean=tmp[1], sd=tmp[2]),2)

# Add prior predictive outcome: see https://opensource.nibr.com/bamdd/src/02ab_meta_analysis_trtdiff.html
tmp = predict(object=M1, newdata=tibble(id="New study", se=10e-19), allow_new_levels=TRUE, sample_new_levels = "gaussian")
points(tmp[1],-1,  pch=16, col=cols[3])
segments(x0 = tmp[3], x1 = tmp[4], y0=-1, y1=-1, lwd=2, col=cols[3])
text(-5, -1, labels = substitute(paste(bold('PIII predictive'))), cex=0.75, adj =0, xpd=TRUE)

# Prior predictive probabilities
round(100*pnorm(c(0,-1,-2,-3), mean=tmp[1], sd = tmp[2]),2)
round(qnorm(c(.25,.5,.75), mean=tmp[1], sd=tmp[2]),2)


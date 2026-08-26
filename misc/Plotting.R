library(tidyverse)
prior_pred %>% tibble %>% select(1:8)

lib = "/Users/jsxo/Library/CloudStorage/OneDrive-NovoNordisk/Projects/Migraine PoS Tutorial/beamer/figures/"

pdf(file = paste0(lib,"prior_predictive.pdf"), width = 8, height=4)
par(mfrow=c(1,2), mar=c(4,4,1,1))
plot(NA, xlim=c(0,30), ylim=c(0,.1), bty='n', xlab="Baseline MMD", ylab="pdf", main="")
lines(density(prior_pred$baseline_em), col=jsxo::add_alpha('steelblue',.75), lwd=2)
lines(density(prior_pred$baseline_cm), col=jsxo::add_alpha('tomato2',.75), lwd=2)
legend("topright", c("EM","CM"), fill=c("steelblue","tomato2"), bty='n')

plot(NA, xlim=c(0,10), ylim=c(0,1.5), bty='n', xlab="MMD reduction", ylab="pdf", main="")
lines(density(prior_pred$pbo_red_em), col=jsxo::add_alpha('steelblue',.75), lwd=2)
lines(density(prior_pred$pbo_red_cm), col=jsxo::add_alpha('tomato2',.75), lwd=2)
# legend("topright", c("EM","CM"), fill=c("steelblue","tomato2"), bty='n')

# plot(NA, xlim=c(0,30), ylim=c(0,1.5), bty='n', xlab="Placebo reduction", ylab="pdf", main="Prior predictive")
lines(density(prior_pred$active_red_em), col=jsxo::add_alpha('steelblue',.75), lwd=2, lty=3)
lines(density(prior_pred$active_red_cm), col=jsxo::add_alpha('tomato2',.75), lwd=2, lty=3)
legend("topright", c("EM placebo","CM placebo","EM active","CM active"), col=c("steelblue","tomato2"), lty=c(1,1,3,3), bty='n', lwd=2)
dev.off()






f <- function(d,idx = NULL){
  if(is.null(idx)){
    idx = sample(1:nrow(draws), size = 1)
  }
  post = draws[idx,]
  em = post$e0_em+post$emax_em*emax_g(d,post$ed50)
  cm = post$e0_cm+post$emax_cm*emax_g(d,post$ed50)
  out = tibble(idx = idx, dose = d, em = em, cm=cm)
  return(out)
}

ds = seq(0,300,1)
nsim = 10e3
set.seed(1234)
lapply(sample(1:nrow(draws), size=nsim, replace=TRUE), f, d=ds) %>%
  bind_rows %>%
  group_by(dose) %>%
  summarize(em05 = quantile(em, probs = 0.05),
            em25 = quantile(em, probs = 0.25),
            em50 = quantile(em, probs = 0.50),
            em75 = quantile(em, probs = 0.75),
            em95 = quantile(em, probs = 0.95),
            cm05 = quantile(cm, probs = 0.05),
            cm25 = quantile(cm, probs = 0.25),
            cm50 = quantile(cm, probs = 0.50),
            cm75 = quantile(cm, probs = 0.75),
            cm95 = quantile(cm, probs = 0.95)) -> tmp

polygonit <- function(x, y1,y2, col='gray60'){
  polygon(c(x,rev(x)), c(y1,rev(y2)), border=NA, col=col)
}

pdf(file = paste0(lib,"emax_posterior.pdf"), width = 6, height=6)
plot(NA, xlim=c(0,300), ylim=c(-10,0), bty='n', xlab="Dose",ylab="Change from baseline")
x = tmp$dose
polygonit(tmp$dose, tmp$em05, tmp$em95, col=jsxo::add_alpha('steelblue',.5))
polygonit(tmp$dose, tmp$em25, tmp$em75, col=jsxo::add_alpha('steelblue',.5))
lines(tmp$dose, tmp$em50, col=jsxo::add_alpha('steelblue',.95), lwd=2)
polygonit(tmp$dose, tmp$cm05, tmp$cm95, col=jsxo::add_alpha('tomato2',.5))
polygonit(tmp$dose, tmp$cm25, tmp$cm75, col=jsxo::add_alpha('tomato2',.5))
lines(tmp$dose, tmp$cm50, col=jsxo::add_alpha('tomato2',.95), lwd=2)

legend("topright", c("EM","CM"), fill=c("steelblue","tomato2"), bty='n')

dev.off()







library(readr)
meta_sum <- read_csv("analysis/output/historical_meta_summary.csv") %>% tibble
emax_sum <- read_csv("analysis/output/historical_emax_summary.csv") %>% tibble

meta_sum %>% slice(c(1:2,7:10))
emax_sum %>% slice(c(1:2,7:10))


pdf(file = paste0(lib,"posterior_forest.pdf"), width = 6, height=6)

par(mar=c(4,5,1,1))
lbls = c("Baseline EM","Baseline CM", "Placebo EM", "Placebo CM", "Drug EM", "Drug CM")
plot(NA, xlim=c(0,30), ylim=c(6.5,1), bty='n', xlab="MMD scale (baseline/change from baseline)", ylab="", yaxt='n')
axis(side=2, at=1:6, labels = lbls, las=1, cex.axis=.75)

off = 0.1
tmp = meta_sum %>% slice(c(1:2,7:10))
points(tmp$median, 1:6-off, col='steelblue', pch=16)
segments(tmp$q5,1:6-off,tmp$q95,1:6-off, lend='butt', lwd=2, col='steelblue')

tmp = emax_sum %>% slice(c(1:2,8:11))
points(tmp$median, 1:6+off, col='tomato2', pch=16)
segments(tmp$q5,1:6+off,tmp$q95,1:6+off, lend='butt', lwd=2, col='tomato2')

legend("bottomright", c("Meta analysis prior", "HOPE Emax posterior"), fill=c("steelblue", "tomato2"), bty='n')
dev.off()









f <- function(d,idx = NULL){
  if(is.null(idx)){
    idx = sample(1:nrow(draws), size = 1)
  }
  post = draws[idx,]
  em = -post$emax_em*emax_g(d,post$ed50)
  cm = -post$emax_cm*emax_g(d,post$ed50)
  out = tibble(idx = idx, dose = d, em = em, cm=cm)
  return(out)
}

ds = seq(0,300,1)
nsim = 10e3
set.seed(1234)
lapply(sample(1:nrow(draws), size=nsim, replace=TRUE), f, d=ds) %>%
  bind_rows %>%
  group_by(dose) %>%
  summarize(em05 = quantile(em, probs = 0.05),
            em25 = quantile(em, probs = 0.25),
            em50 = quantile(em, probs = 0.50),
            em75 = quantile(em, probs = 0.75),
            em95 = quantile(em, probs = 0.95),
            cm05 = quantile(cm, probs = 0.05),
            cm25 = quantile(cm, probs = 0.25),
            cm50 = quantile(cm, probs = 0.50),
            cm75 = quantile(cm, probs = 0.75),
            cm95 = quantile(cm, probs = 0.95)) -> tmp

polygonit <- function(x, y1,y2, col='gray60'){
  polygon(c(x,rev(x)), c(y1,rev(y2)), border=NA, col=col)
}

pdf(file = paste0(lib,"effect_posterior.pdf"), width = 8, height=6)
plot(NA, xlim=c(0,300), ylim=c(0,5), bty='n', xlab="Dose",ylab="Drug effect")
x = tmp$dose
polygonit(tmp$dose, tmp$cm05, tmp$cm95, col=jsxo::add_alpha('tomato2',.35))
polygonit(tmp$dose, tmp$cm25, tmp$cm75, col=jsxo::add_alpha('tomato2',.35))
lines(tmp$dose, tmp$cm50, col=jsxo::add_alpha('tomato2',.99), lwd=2)

polygonit(tmp$dose, tmp$em05, tmp$em95, col=jsxo::add_alpha('steelblue',.35))
polygonit(tmp$dose, tmp$em25, tmp$em75, col=jsxo::add_alpha('steelblue',.35))
lines(tmp$dose, tmp$em50, col=jsxo::add_alpha('steelblue',.99), lwd=2)
abline(h=c(2,3,4), lty=3, col='gray50')

legend("topright", c("EM","CM"), fill=c("steelblue","tomato2"), bty='n')

dev.off()








plot_sig <- function(tmp, doses = c(20,200), sig_threshold = qnorm(0.975)*5*sqrt(2/233)){
  thrsh = 0:6
  par(mfrow=c(2,length(doses)), mar=c(3,2,2,1), oma=c(1,1,2,1), bty='n')
  for(d in doses){
    plot(tmp$true_effect, tmp$obs_effect, col=ifelse(tmp$significant, "steelblue","tomato2"), pch=ifelse(tmp$dose==d, ".",""),
         xlab="", ylab="", main=paste("Dose",d), xlim=c(0,6), ylim=c(0,6))
    if(d==doses[1]) mtext("Obs effect", side=2, line=2, cex=.75)
    mtext("True effect", side=1, line=2, cex=.75)
    abline(0,1, lty=1)
    abline(h=sig_threshold, col='gray50')
    abline(h=thrsh, v=thrsh, col='gray70',lty=3)
  }
  for(d in doses){
    plot(tmp$y_pbo,tmp$y_act, col=ifelse(tmp$significant, "steelblue","tomato2"), pch=ifelse(tmp$dose==d, ".",""),
         xlab="", ylab="", main=paste("Dose",d), xlim=c(-10,0), ylim=c(-10,0))
    if(d==doses[1]) mtext(expression(Y[active]), side=2, line=2, cex=.75)
    mtext(expression(Y[placebo]), side=1, line=2, cex=.75)
    abline(0,1, lty=1)
    abline(-sig_threshold,1, lty=1, col='gray50')
    map(thrsh, function(x) abline(-x,1, lty=3, col='gray70'))
  }
}

pdf(file = paste0(lib,"prog1_sims_em.pdf"), width = 6, height=6.5)
  prog1$raw %>% tibble %>% filter(pop=="EM") %>% plot_sig
  mtext("EM", outer=TRUE, line=0, font=2, cex=1.2)
dev.off()

pdf(file = paste0(lib,"prog1_sims_cm.pdf"), width = 6, height=6.5)
  prog1$raw %>% tibble %>% filter(pop=="CM") %>% plot_sig
  mtext("CM", outer=TRUE, line=0, font=2, cex=1.2)
dev.off()



tmp = prog1$raw

## Program 1 success table
tmp %>% group_by(sim, pop) %>%
  summarize(across(starts_with("sig_and_obs_effect_ge"),mean)) %>%
  ungroup %>%
  mutate(across(starts_with("sig_and_obs"), \(x) 1*(x>0))) %>%
  pivot_longer(cols = starts_with("sig_and_obs"), values_to = "success", names_to = "thresh") %>%
  mutate(thresh = str_replace_all(thresh, "sig_and_obs_effect_ge_","")) %>%
  pivot_wider(values_from = success, names_from = pop) %>%
  mutate(both = CM*EM) %>%
  group_by(thresh) %>%
  summarize(CM = mean(CM), EM = mean(EM), both = mean(both)) %>%
  xtable::xtable() %>% print(include.rownames=FALSE)





pdf(file = paste0(lib,"prog2_sims_em.pdf"), width = 9, height=6.5)
  prog2$phase3 %>% filter(pop=="EM") %>% rename(dose=selected_dose) %>% plot_sig(doses=c(20,100,300))
  mtext("EM", outer=TRUE, line=0, font=2, cex=1.2)
dev.off()

pdf(file = paste0(lib,"prog2_sims_cm.pdf"), width = 9, height=6.5)
  prog2$phase3 %>% filter(pop=="CM") %>% rename(dose=selected_dose) %>% plot_sig(doses=c(20,100,300))
  mtext("CM", outer=TRUE, line=0, font=2, cex=1.2)
dev.off()






## Program 2 success table by dose and threshold
prog2$phase3 %>% tibble %>% group_by(sim, selected_dose, pop) %>%
  summarize(across(starts_with("sig_and_obs_effect_ge"),mean)) %>%
  ungroup %>%
  mutate(across(starts_with("sig_and_obs"), \(x) 1*(x>0))) %>%
  pivot_longer(cols = starts_with("sig_and_obs"), values_to = "success", names_to = "thresh") %>%
  mutate(thresh = str_replace_all(thresh, "sig_and_obs_effect_ge_","")) %>%
  pivot_wider(values_from = success, names_from = pop) %>%
  mutate(both = CM*EM) %>%
  group_by(selected_dose, thresh) %>%
  summarize(CM = mean(CM), EM = mean(EM), both = mean(both)) %>%
  xtable::xtable() %>% print(include.rownames=FALSE)

## Program 2 success table by threshold
prog2$phase3 %>% tibble %>% group_by(sim, selected_dose, pop) %>%
  summarize(across(starts_with("sig_and_obs_effect_ge"),mean)) %>%
  ungroup %>%
  mutate(across(starts_with("sig_and_obs"), \(x) 1*(x>0))) %>%
  pivot_longer(cols = starts_with("sig_and_obs"), values_to = "success", names_to = "thresh") %>%
  mutate(thresh = str_replace_all(thresh, "sig_and_obs_effect_ge_","")) %>%
  pivot_wider(values_from = success, names_from = pop) %>%
  mutate(both = CM*EM) %>%
  group_by(thresh) %>%
  summarize(CM = mean(CM), EM = mean(EM), both = mean(both)) %>%
  xtable::xtable() %>% print(include.rownames=FALSE)







tmp = tibble(thresh=2:4,
             prog1=prog1$threshold_program_success[,4],
             prog2=prog2$program_success[-1,5])

tmp
layout(1)
off = .15
lwds = 40

pdf(file = paste0(lib,"prog_compare.pdf"), width = 6, height=6)
  par(mar=c(2,2.5,1,1))
  plot(NA, xlim=c(1.5,4.5), ylim=c(0,.4), bty='n', axes=FALSE)
  mtext("Probability", 2, 2.5)
  axis(1, at=2:4, labels= paste("Effect >",2:4, "MMD"))
  axis(2, at=seq(0,1,.1), labels = 100*seq(0,1,.1), las=1)
  lines(tmp$thresh-off, tmp$prog1, type='h', lwd=lwds, lend='butt', col='steelblue')
  lines(tmp$thresh+off, tmp$prog2, type='h', lwd=lwds, lend='butt', col='tomato2')
  legend("topright", c("Program 1", "Program 2"), fill=c("steelblue","tomato2"), bty='n')
dev.off()

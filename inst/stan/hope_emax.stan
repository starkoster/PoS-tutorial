functions {
  #include clinical_functions.stan
}

data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  vector<lower=0>[N] dose;
  vector<lower=0, upper=1>[N] w_em;

  int<lower=0> NB;
  vector[NB] y_b;
  vector<lower=0>[NB] se_b;
  vector<lower=0, upper=1>[NB] w_em_b;

  // Joint historical predictive prior plus an independent log-ED50 prior.
  vector[6] prior_mean;
  cov_matrix[6] prior_cov;
  real ed50_logmean;
  real<lower=0> ed50_logsd;
}

transformed data {
  matrix[6, 6] prior_L = cholesky_decompose(prior_cov);
}

parameters {
  vector[6] z_eta;
  real log_ed50;
}

transformed parameters {
  vector[6] eta = prior_mean + prior_L * z_eta;
  vector[8] clinical = clinical_from_eta(eta);
  real<lower=0> ed50 = exp(log_ed50);
  vector[N] mean_change;
  vector[NB] mean_baseline;

  // Both populations share potency, while their placebo and Emax differ.
  for (i in 1:N) {
    real fraction_of_maximum = dose[i] / (dose[i] + ed50);

    mean_change[i] = w_em[i] * (clinical[3] + clinical[5] * fraction_of_maximum)
                  + (1 - w_em[i]) * (clinical[4] + clinical[6] * fraction_of_maximum);
  }

  if (NB > 0) {
    for (i in 1:NB) {
      mean_baseline[i] = w_em_b[i] * clinical[1] + (1 - w_em_b[i]) * clinical[2];
    }
  }
}

model {
  z_eta ~ std_normal();
  log_ed50 ~ normal(ed50_logmean, ed50_logsd);

  // Reported arm-summary SEs are treated as known.
  y ~ normal(mean_change, se);
  y_b ~ normal(mean_baseline, se_b);
}

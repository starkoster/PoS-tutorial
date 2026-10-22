functions {
  #include clinical_functions.stan
}

data {
  // Change summaries, including both pure and mixed EM/CM arms.
  int<lower=1> N;
  int<lower=1> J;
  array[N] int<lower=1, upper=J> study;
  vector[N] y;
  vector<lower=0>[N] se;
  vector<lower=0, upper=1>[N] w_em;
  array[N] int<lower=0, upper=1> active;
  vector<lower=0, upper=1>[N] weight;

  // Baseline summaries retained after checking their subgroup provenance.
  int<lower=0> NB;
  array[NB] int<lower=1, upper=J> study_b;
  vector[NB] y_b;
  vector<lower=0>[NB] se_b;
  vector<lower=0, upper=1>[NB] w_em_b;

  // Calibrated latent hyperpriors.
  vector[6] m;
  vector<lower=0>[6] s;
  vector<lower=0>[6] q;
  real<lower=0> lkj_shape;
}

parameters {
  vector[6] mu_eta;
  vector<lower=0>[6] tau_eta;
  cholesky_factor_corr[6] L_Omega;
  array[J] vector[6] z_eta;
}

transformed parameters {
  matrix[6, 6] L_eta = diag_pre_multiply(tau_eta, L_Omega);
  array[J] vector[6] eta_study;
  array[J] vector[8] clinical_study;
  vector[N] mean_change;
  vector[NB] mean_baseline;

  // Noncentered study effects improve sampling near small heterogeneity SDs.
  for (j in 1:J) {
    eta_study[j] = mu_eta + L_eta * z_eta[j];
    clinical_study[j] = clinical_from_eta(eta_study[j]);
  }

  // An active arm adds its placebo-adjusted effect to the placebo change.
  for (i in 1:N) {
    vector[8] x = clinical_study[study[i]];

    mean_change[i] = w_em[i] * (x[3] + active[i] * x[5])
                  + (1 - w_em[i]) * (x[4] + active[i] * x[6]);
  }

  if (NB > 0) {
    for (i in 1:NB) {
      vector[8] x = clinical_study[study_b[i]];

      mean_baseline[i] = w_em_b[i] * x[1] + (1 - w_em_b[i]) * x[2];
    }
  }
}

model {
  // The lower bound makes each tau prior half-normal.
  mu_eta ~ normal(m, s);
  tau_eta ~ normal(0, q);
  L_Omega ~ lkj_corr_cholesky(lkj_shape);

  for (j in 1:J) {
    z_eta[j] ~ std_normal();
  }

  // Discount active cross-drug information using a power likelihood.
  for (i in 1:N) {
    target += weight[i] * normal_lpdf(y[i] | mean_change[i], se[i]);
  }

  y_b ~ normal(mean_baseline, se_b);
}

generated quantities {
  // Predict a new study, retaining the full study-level heterogeneity.
  vector[6] eta_new = multi_normal_cholesky_rng(mu_eta, L_eta);
}

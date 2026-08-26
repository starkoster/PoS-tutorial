functions {
  vector pars_from_eta(vector eta) {
    vector[12] out;

    real baseline_em;
    real baseline_cm;

    real active_red_em;
    real active_red_cm;

    real pbo_red_em;
    real pbo_red_cm;

    real drug_red_em;
    real drug_red_cm;

    baseline_em = 30 * inv_logit(eta[1]);
    baseline_cm = baseline_em + (30 - baseline_em) * inv_logit(eta[2]);

    active_red_em = baseline_em * inv_logit(eta[3]);
    active_red_cm = active_red_em + (baseline_cm - active_red_em) * inv_logit(eta[4]);

    pbo_red_em = active_red_em * inv_logit(eta[5]);
    pbo_red_cm = pbo_red_em + (active_red_cm - active_red_em) * inv_logit(eta[6]);

    drug_red_em = active_red_em - pbo_red_em;
    drug_red_cm = active_red_cm - pbo_red_cm;

    out[1] = baseline_em;
    out[2] = baseline_cm;

    out[3] = pbo_red_em;
    out[4] = pbo_red_cm;

    out[5] = drug_red_em;
    out[6] = drug_red_cm;

    out[7] = -pbo_red_em;
    out[8] = -pbo_red_cm;

    out[9] = -drug_red_em;
    out[10] = -drug_red_cm;

    out[11] = active_red_em;
    out[12] = active_red_cm;

    return out;
  }
}

data {
  int<lower=1> N;

  vector[N] y;
  vector<lower=0>[N] se;
  vector<lower=0>[N] dose;
  vector<lower=0, upper=1>[N] w_em;

  // Optional baseline likelihood for the new trial.
  // These observations are on the absolute baseline MMD scale.
  int<lower=0> N_base;
  vector[N_base] y_base;
  vector<lower=0>[N_base] se_base;
  vector<lower=0, upper=1>[N_base] w_em_base;

  vector[6] eta_prior_mu;
  cov_matrix[6] eta_prior_Sigma;

  real ed50_prior_logmean;
  real<lower=0> ed50_prior_logsd;

  int<lower=1> N_pred;
  vector<lower=0>[N_pred] dose_pred;
}

transformed data {
  matrix[6, 6] L_eta_prior;
  L_eta_prior = cholesky_decompose(eta_prior_Sigma);
}

parameters {
  vector[6] z_eta;
  real log_ed50;
}

transformed parameters {
  vector[6] eta;
  vector[12] pars;

  real baseline_em;
  real baseline_cm;

  real e0_em;
  real e0_cm;

  real emax_em;
  real emax_cm;

  real<lower=0> ed50;

  vector[N] mu;
  vector[N_base] mu_base;

  eta = eta_prior_mu + L_eta_prior * z_eta;
  pars = pars_from_eta(eta);

  baseline_em = pars[1];
  baseline_cm = pars[2];

  e0_em = pars[7];
  e0_cm = pars[8];

  emax_em = pars[9];
  emax_cm = pars[10];

  ed50 = exp(log_ed50);

  for (i in 1:N) {
    real g;
    real mu_em;
    real mu_cm;

    g = dose[i] / (dose[i] + ed50);

    mu_em = e0_em + emax_em * g;
    mu_cm = e0_cm + emax_cm * g;

    mu[i] = w_em[i] * mu_em + (1 - w_em[i]) * mu_cm;
  }

  if (N_base > 0) {
    for (i in 1:N_base) {
      mu_base[i] = w_em_base[i] * baseline_em + (1 - w_em_base[i]) * baseline_cm;
    }
  }
}

model {
  z_eta ~ std_normal();
  log_ed50 ~ normal(ed50_prior_logmean, ed50_prior_logsd);

  y ~ normal(mu, se);

  if (N_base > 0) {
    y_base ~ normal(mu_base, se_base);
  }
}

generated quantities {
  real pbo_red_em;
  real pbo_red_cm;

  real drug_red_em;
  real drug_red_cm;

  real active_red_em;
  real active_red_cm;

  vector[N] y_rep;
  vector[N_base] y_base_rep;

  vector[N_pred] curve_em;
  vector[N_pred] curve_cm;
  vector[N_pred] curve_mix_equal;

  vector[N_pred] post_mmd_em;
  vector[N_pred] post_mmd_cm;

  vector[N_pred] effect_em_vs_placebo;
  vector[N_pred] effect_cm_vs_placebo;

  pbo_red_em = pars[3];
  pbo_red_cm = pars[4];

  drug_red_em = pars[5];
  drug_red_cm = pars[6];

  active_red_em = pars[11];
  active_red_cm = pars[12];

  for (i in 1:N) {
    y_rep[i] = normal_rng(mu[i], se[i]);
  }

  if (N_base > 0) {
    for (i in 1:N_base) {
      y_base_rep[i] = normal_rng(mu_base[i], se_base[i]);
    }
  }

  for (m in 1:N_pred) {
    real g_pred;

    g_pred = dose_pred[m] / (dose_pred[m] + ed50);

    curve_em[m] = e0_em + emax_em * g_pred;
    curve_cm[m] = e0_cm + emax_cm * g_pred;

    curve_mix_equal[m] = 0.5 * curve_em[m] + 0.5 * curve_cm[m];

    post_mmd_em[m] = baseline_em + curve_em[m];
    post_mmd_cm[m] = baseline_cm + curve_cm[m];

    effect_em_vs_placebo[m] = e0_em - curve_em[m];
    effect_cm_vs_placebo[m] = e0_cm - curve_cm[m];
  }
}

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
  int<lower=1> N_change;
  int<lower=1> J;

  array[N_change] int<lower=1, upper=J> study_change;

  vector[N_change] y_change;
  vector<lower=0>[N_change] se_change;
  vector<lower=0, upper=1>[N_change] w_em;

  // 0 for placebo, 1 for historical active-reference arm.
  // Historical numerical dose values are intentionally not used.
  vector<lower=0, upper=1>[N_change] active_ref;

  // 1 for placebo rows. Active rows are normally excluded or down-weighted.
  vector<lower=0, upper=1>[N_change] weight_change;

  int<lower=0> N_base;
  array[N_base] int<lower=1, upper=J> study_base;
  vector[N_base] y_base;
  vector<lower=0>[N_base] se_base;
  vector<lower=0, upper=1>[N_base] w_em_base;

  vector[6] mu_eta_prior;
  vector<lower=0>[6] sd_mu_eta_prior;
  vector<lower=0>[6] tau_eta_scale;

  real<lower=0> lkj_shape;
}

parameters {
  vector[6] mu_eta;

  vector<lower=0>[6] tau_eta;
  cholesky_factor_corr[6] L_Omega;

  array[J] vector[6] z_eta;
}

transformed parameters {
  matrix[6, 6] L_eta;

  array[J] vector[6] eta_study;
  array[J] vector[12] pars_study;

  vector[N_change] mu_change;
  vector[N_base] mu_base;

  L_eta = diag_pre_multiply(tau_eta, L_Omega);

  for (j in 1:J) {
    eta_study[j] = mu_eta + L_eta * z_eta[j];
    pars_study[j] = pars_from_eta(eta_study[j]);
  }

  for (i in 1:N_change) {
    vector[12] ps;
    real e0_em;
    real e0_cm;
    real emax_em;
    real emax_cm;
    real mu_em;
    real mu_cm;

    ps = pars_study[study_change[i]];

    e0_em = ps[7];
    e0_cm = ps[8];
    emax_em = ps[9];
    emax_cm = ps[10];

    mu_em = e0_em + emax_em * active_ref[i];
    mu_cm = e0_cm + emax_cm * active_ref[i];

    mu_change[i] = w_em[i] * mu_em + (1 - w_em[i]) * mu_cm;
  }

  if (N_base > 0) {
    for (i in 1:N_base) {
      vector[12] ps;
      real baseline_em;
      real baseline_cm;

      ps = pars_study[study_base[i]];
      baseline_em = ps[1];
      baseline_cm = ps[2];

      mu_base[i] = w_em_base[i] * baseline_em + (1 - w_em_base[i]) * baseline_cm;
    }
  }
}

model {
  mu_eta ~ normal(mu_eta_prior, sd_mu_eta_prior);
  tau_eta ~ normal(0, tau_eta_scale);
  L_Omega ~ lkj_corr_cholesky(lkj_shape);

  for (j in 1:J) {
    z_eta[j] ~ std_normal();
  }

  for (i in 1:N_change) {
    target += weight_change[i] * normal_lpdf(y_change[i] | mu_change[i], se_change[i]);
  }

  if (N_base > 0) {
    y_base ~ normal(mu_base, se_base);
  }
}

generated quantities {
  vector[6] eta_center;
  vector[12] pars_center;

  vector[6] eta_new;
  vector[12] pars_new;

  corr_matrix[6] Omega;
  matrix[6, 6] Sigma_eta;

  real baseline_em_new;
  real baseline_cm_new;

  real e0_em_new;
  real e0_cm_new;

  real emax_em_new;
  real emax_cm_new;

  real pbo_red_em_new;
  real pbo_red_cm_new;

  real drug_red_em_new;
  real drug_red_cm_new;

  real active_red_em_new;
  real active_red_cm_new;

  vector[N_change] y_change_rep;

  eta_center = mu_eta;
  pars_center = pars_from_eta(eta_center);

  Omega = multiply_lower_tri_self_transpose(L_Omega);
  Sigma_eta = quad_form_diag(Omega, tau_eta);

  eta_new = multi_normal_cholesky_rng(mu_eta, L_eta);
  pars_new = pars_from_eta(eta_new);

  baseline_em_new = pars_new[1];
  baseline_cm_new = pars_new[2];

  pbo_red_em_new = pars_new[3];
  pbo_red_cm_new = pars_new[4];

  drug_red_em_new = pars_new[5];
  drug_red_cm_new = pars_new[6];

  e0_em_new = pars_new[7];
  e0_cm_new = pars_new[8];

  emax_em_new = pars_new[9];
  emax_cm_new = pars_new[10];

  active_red_em_new = pars_new[11];
  active_red_cm_new = pars_new[12];

  for (i in 1:N_change) {
    y_change_rep[i] = normal_rng(mu_change[i], se_change[i]);
  }
}

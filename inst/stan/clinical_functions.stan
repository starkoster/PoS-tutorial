vector clinical_from_eta(vector eta) {

  vector[6] fraction = inv_logit(eta);
  vector[8] clinical;

  // 1--2: positive baseline burdens, ordered EM < CM < 30.
  clinical[1] = 30 * fraction[1];
  clinical[2] = clinical[1] + (30 - clinical[1]) * fraction[2];

  // 7--8: negative total active changes, bounded by baseline burden.
  clinical[7] = -clinical[1] * fraction[3];
  clinical[8] = clinical[7] - (clinical[2] + clinical[7]) * fraction[4];

  // 3--4: negative placebo changes, smaller than total active reductions.
  clinical[3] = clinical[7] * fraction[5];
  clinical[4] = clinical[3] + (clinical[8] - clinical[7]) * fraction[6];

  // 5--6: maximum placebo-adjusted effects = active minus placebo.
  clinical[5] = clinical[7] - clinical[3];
  clinical[6] = clinical[8] - clinical[4];

  return clinical;
}

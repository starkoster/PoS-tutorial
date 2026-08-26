# MigraineEmaxSim

Bayesian historical meta-analysis and ED50 Emax simulation tools for migraine development programs.

## Data columns

The included `data/hist_trials.csv` has one row per trial arm and includes:

- `nct`, `title`: trial identifiers;
- `CM`: proportion chronic migraine, so `w_em = 1 - CM`;
- `arm`: `placebo` or `active`;
- `dose`: historical max/reference dose, except HOPE where the actual dose is used;
- `n`: arm size;
- `base_mmd`: baseline monthly migraine days;
- `base_mmd_sd`: baseline MMD standard deviation;
- `chg`: change from baseline;
- `chg_se`: standard error for `chg`;
- `chg_sd`: inferred standard deviation for change.

Historical numerical doses are **not** used as dose-response information because they come from different drugs. Historical active arms can optionally be used as discounted active-reference evidence.

## Subjective priors

The only subjective priors in the historical meta-analysis are the six `eta` priors:

```r
priors <- default_eta_priors()
eta_prior_table(priors)
```

They are elicited through clinical medians on the monthly-migraine-day scale, not directly on the abstract Stan scale. The default clinical medians are:

```r
list(
  baseline_em   = 8,
  baseline_cm   = 18,
  active_red_em = 3,
  active_red_cm = 5,
  pbo_red_em    = 1,
  pbo_red_cm    = 2
)
```

These imply median placebo-adjusted drug reductions of 2 MMD in EM and 3 MMD in CM. The model transformation guarantees:

- `baseline_cm > baseline_em`;
- `pbo_red_cm > pbo_red_em`;
- `drug_red_cm > drug_red_em`;
- placebo and active reductions cannot exceed baseline MMD;
- baseline MMD stays within the 0 to 30 days-per-month range.

Baseline MMD observations enter the Stan likelihood as:

```text
base_mmd_i ~ Normal(w_em_i * baseline_em_j + (1 - w_em_i) * baseline_cm_j,
                    base_mmd_sd_i / sqrt(n_i))
```

This makes the observed baseline MMD informative for the maximum allowed placebo and active reductions and for the CM > EM baseline ordering.

## Prior predictive calibration

Run:

```r
priors <- default_eta_priors()
prior_pred <- simulate_eta_prior(priors, n = 200000, seed = 123)
prior_predictive_summary(prior_pred, threshold = 15)
```

If the prior puts too much probability on reductions above 15 MMD, tighten either `sd_mu_eta_prior` or `tau_eta_scale`, or reduce the clinical medians.

For a weaker sensitivity analysis:

```r
priors <- weak_eta_priors()
prior_predictive_summary(simulate_eta_prior(priors), threshold = 15)
```

## Full analysis

```r
library(MigraineEmaxSim)
res <- run_full_analysis(root = ".", n_program_sims = 5000)
```

This writes CSV summaries to `analysis/output`, including:

- `meta_analysis_prior_specification.csv`;
- `prior_predictive_checks.csv`;
- `historical_meta_summary.csv`;
- `historical_emax_summary.csv`;
- development-program simulation summaries.

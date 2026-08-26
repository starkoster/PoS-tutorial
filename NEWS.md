# MigraineEmaxSim 0.3.1

* Made `default_eta_priors()` specify placebo-adjusted drug reductions (`drug_red_em`, `drug_red_cm`) directly in the default prior medians, then compute total active reductions internally.
* Added `normalize_prior_medians()` to validate clinical prior medians and provide clearer messages when `active_red - pbo_red` does not satisfy CM > EM drug reduction.
* This avoids ambiguity between total active/reference reduction and placebo-adjusted drug reduction.

# MigraineEmaxSim 0.3.0

* Added support for `base_mmd` and `base_mmd_sd` columns.
* Baseline MMD is now included as a likelihood contribution in both `historical_meta.stan` and `historical_emax.stan`.
* Baseline SDs are converted to standard errors using `base_mmd_sd / sqrt(n)`.
* Added `eta_prior_table()` and `weak_eta_priors()` to make the meta-analysis priors explicit and easier to calibrate clinically.
* `run_full_analysis()` now writes `meta_analysis_prior_specification.csv`.

# MigraineEmaxSim 0.2.1

- Fixed use of `stats::tapply()` in program simulation summaries. `tapply()` is in base R, not the stats namespace.
- Updated `extract_emax_draws()` to coerce posterior draws to a plain data frame before subsetting, avoiding the harmless `draws_df` metadata warning.

# MigraineEmaxSim 0.2.0

- Initial package build for historical meta-analysis, HOPE Emax fitting, and development program simulations.

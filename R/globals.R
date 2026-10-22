# Declare names used in tidy evaluation for package-check tools.
# This does not change evaluation or attach packages to the user's session.
utils::globalVariables(c(
  "CM", "EM", "active_mean", "base_mmd_sd", "baseline_se", "both",
  "both_same_dose", "cm_count", "count", "delta", "density", "direct",
  "dose", "ed50", "estimand", "fraction_of_go", "gated_from_start",
  "given_go", "go", "lower", "median", "n", "population", "probability",
  "program", "q05", "q50", "q95", "quantile", "quantity", "selected_dose",
  "significant", "source", "stage", "strategy", "success", "target",
  "threshold", "title", "upper", "use_baseline", "value"
))

source(file.path(if (dir.exists("analysis")) "analysis" else ".", "_setup.R"))

result <- run_full_analysis(
  root = pkg_root,
  csv_path = file.path(pkg_root, "data", "hist_trials.csv"),
  output_dir = out_dir,
  n_program_sims = 5000,
  use_historical_active = TRUE,
  active_borrow_weight = 0.50,
  ed50_median = 20,
  ed50_lower = 10,
  ed50_upper = 100
)

print(result$prior_checks)
print(result$prog1$significance_only)
print(result$prog1$threshold_program_success)
print(result$prog2$phase2_summary)
print(result$prog2$selection_summary)
print(result$prog2$program_success)



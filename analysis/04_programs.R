library(MigraineEmaxSim)

# Check the version actually loaded in this R session, before doing any work.
if (utils::compareVersion(as.character(getNamespaceVersion("MigraineEmaxSim")),
                          "0.5.3") < 0) {
  stop("Install MigraineEmaxSim >= 0.5.3 and restart R before running this script.",
       call. = FALSE)
}

output_dir <- "results_R_readable"

# This also accepts existing posterior tables with one-column matrix quantities.
draws <- readRDS(file.path(output_dir, "03_emax_draws.rds")) |>
  MigraineEmaxSim::as_numeric_draws()

# STEP 1: inspect the prospective Phase 2 and Phase 3 designs.
# Edit the named entries here if a different program is being considered.
design <- default_program_design()
print(design)

# STEP 2: fix the MCP candidate shapes from current knowledge, before future data.
ed50_candidates <- stats::quantile(draws$ed50, c(0.10, 0.50, 0.90), names = FALSE)

mcp <- make_mcp_test(
  doses = design$p2_doses,
  ed50_candidates = ed50_candidates,
  n_per_arm = design$p2_n,
  sigma = design$sigma,
  alpha = design$p2_alpha
)

# STEP 3: simulate whole programs. Inside this function, the main loop calls
# simulate_phase2() and then simulate_phase3() with the SAME posterior truth.
# For a quick review run, change n_sim to 1000L before the full simulation.
programs <- simulate_programs(
  draws = draws,
  mcp = mcp,
  n_sim = 50000L,
  seed = 4601,
  design = design
)

# STEP 4: report success from program start and conditional on Phase 2 Go.
# The paired Program 1 screening analysis keeps its Phase 3 design unchanged.
summary <- summarize_programs(programs)

print(summary$probabilities)
print(summary$derisking)

saveRDS(programs, file.path(output_dir, "04_program_records.rds"))

readr::write_csv(summary$probabilities, file.path(output_dir, "program_probabilities.csv"))
readr::write_csv(summary$derisking, file.path(output_dir, "fixed_design_screening.csv"))
readr::write_csv(summary$selection, file.path(output_dir, "dose_selection.csv"))
readr::write_csv(summary$participants, file.path(output_dir, "expected_participants.csv"))
readr::write_csv(programs$phase2, file.path(output_dir, "phase2_decisions.csv"))

save_analysis_plot(
  plot_programs(summary), "04_program_comparison", file.path(output_dir, "figures")
)

save_analysis_plot(
  plot_derisking(summary), "04_fixed_design_screening", file.path(output_dir, "figures")
)

save_analysis_plot(
  plot_selection(summary), "04_dose_selection", file.path(output_dir, "figures"),
  width = 6, height = 4
)


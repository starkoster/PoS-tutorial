# Optional line-by-line example after Stage 3. This is separate from run_all.R.
library(MigraineEmaxSim)

draws <- readRDS("results_R_readable/03_emax_draws.rds")
design <- default_program_design()

mcp <- make_mcp_test(
  doses = design$p2_doses,
  ed50_candidates = stats::quantile(draws$ed50, c(0.1, 0.5, 0.9), names = FALSE),
  n_per_arm = design$p2_n,
  sigma = design$sigma,
  alpha = design$p2_alpha
)

# One joint posterior draw is the unknown truth throughout this program.
set.seed(101)
truth <- draws[sample.int(nrow(draws), 1), ]

# Program 1: go directly to the fixed Phase 3 design.
direct_phase3 <- simulate_phase3(
  truth, design$p1_doses, design$p1_n, sigma = design$sigma
)

print(phase3_success(direct_phase3, threshold = -2))

# Program 2, first step: simulate the mixed-population Phase 2 trial.
phase2 <- simulate_phase2(truth, mcp, design$p2_em, design$p2_go_delta)

print(phase2$arms)
print(phase2$go)

# Program 2, second step: simulate selected-dose Phase 3 only after Go.
if (phase2$go) {

  selected_phase3 <- simulate_phase3(
    truth,
    doses = c(0, phase2$selected_dose),
    n_per_arm = design$p2_program_p3_n,
    sigma = design$sigma
  )

  print(phase3_success(selected_phase3, threshold = -2))
}

# Program 1's secondary analysis applies this Go decision to its SAME Phase 3.
direct_success <- phase3_success(direct_phase3, threshold = -2)
screened_success <- direct_success |>
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x & phase2$go))

print(screened_success)

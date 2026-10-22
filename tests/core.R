library(MigraineEmaxSim)

# Clinical constraints and the inverse transformation must agree numerically.
set.seed(801)
eta <- matrix(rnorm(6000), ncol = 6)
clinical <- clinical_from_eta(eta)

stopifnot(
  all(clinical$baseline_cm > clinical$baseline_em),
  all(clinical$baseline_cm < 30),
  all(clinical$placebo_cm < clinical$placebo_em),
  all(clinical$placebo_em < 0),
  all(clinical$emax_cm < clinical$emax_em),
  all(clinical$emax_em < 0),
  all(clinical$baseline_em + clinical$active_em > 0),
  all(clinical$baseline_cm + clinical$active_cm > 0)
)

for (row in seq_len(100)) {
  recovered <- eta_from_reference(as.list(clinical[row, ]))
  stopifnot(max(abs(recovered - eta[row, ])) < 1e-8)
}

# A lognormal elicitation must reproduce all three requested quantiles.
ed50 <- make_ed50_prior()

stopifnot(
  abs(qlnorm(0.05, ed50$logmean, ed50$logsd) - ed50$q05) < 1e-8,
  abs(qlnorm(0.50, ed50$logmean, ed50$logsd) - ed50$q50) < 1e-8,
  abs(qlnorm(0.95, ed50$logmean, ed50$logsd) - ed50$q95) < 1e-8,
  inherits(try(make_ed50_prior(10, 20, 100), silent = TRUE), "try-error")
)

# Correct only documented HOPE inputs and omit the pooled LY baselines.
arms <- read_analysis_arms()
hope <- dplyr::filter(arms, title == "HOPE")
placebo <- dplyr::filter(hope, dose == 0)

stopifnot(
  placebo$base_mmd_sd == 4.6,
  placebo$input_baseline_sd == 4.3,
  abs(placebo$CM - 62 / 94) < 1e-12,
  all(!arms$use_baseline[arms$title == "LY3451838"]),
  sum(arms$use_baseline[arms$title != "HOPE"]) == 6
)

# Package-based MCP contrasts must point towards negative effects.
design <- default_program_design()
mcp <- make_mcp_test(design$p2_doses, c(8, 20, 50))

decreasing <- -6 * mcp$doses / (mcp$doses + 20)

test_response <- function(response) {
  DoseFinding::MCTtest(
    dose = mcp$doses,
    resp = response,
    models = mcp$models,
    S = mcp$covariance,
    type = "general",
    df = Inf,
    alpha = mcp$alpha,
    alternative = "one.sided",
    contMat = mcp$contrasts,
    critV = mcp$critical_value,
    pVal = FALSE
  )
}

benefit <- test_response(decreasing)
harm <- test_response(-decreasing)

stopifnot(
  max(benefit$tStat) > mcp$critical_value,
  max(harm$tStat) < 0
)

# Validate the DoseFinding result against the known-variance contrast statistic.
manual <- as.vector(decreasing %*% mcp$contrasts) /
  sqrt(diag(t(mcp$contrasts) %*% mcp$covariance %*% mcp$contrasts))

stopifnot(max(abs(benefit$tStat - manual)) < 1e-10)

# Under a constant dose response, the MCP rejection rate should be about 5%.
null <- tibble::tibble(
  placebo_em = -2, placebo_cm = -3,
  emax_em = 0, emax_cm = 0, ed50 = 20
)

set.seed(802)
rejections <- replicate(3000, simulate_phase2(null, mcp)$significant)
stopifnot(abs(mean(rejections) - 0.05) < 0.02)

# Shared placebo means induce approximately 0.5 correlation between comparisons.
# One beneficial two-sided 5% test rejects with probability 2.5% under the null.
set.seed(803)
trials <- replicate(
  3000,
  dplyr::filter(simulate_phase3(null, c(0, 20, 200), 233), population == "EM"),
  simplify = FALSE
)

deltas <- do.call(rbind, lapply(trials, function(x) x$delta))
significant <- vapply(trials, function(x) x$significant[1], logical(1))

stopifnot(
  abs(cor(deltas)[1, 2] - 0.5) < 0.07,
  abs(mean(significant) - 0.025) < 0.015
)

# Pairing is essential: screened Program 1 must use the SAME Phase 3 outcomes.
truth <- tibble::tibble(
  placebo_em = -1, placebo_cm = -3,
  emax_em = -2, emax_cm = -4, ed50 = 20
)

programs <- simulate_programs(truth, mcp, n_sim = 150, seed = 804)

direct <- dplyr::filter(programs$events, strategy == "P1_direct")
screened <- dplyr::filter(programs$events, strategy == "P1_screened")

paired <- dplyr::left_join(
  direct, screened,
  by = c("program", "threshold"),
  suffix = c("_direct", "_screened")
)

stopifnot(
  all(paired$both_screened == (paired$both_direct & paired$go_direct)),
  all(paired$EM_screened == (paired$EM_direct & paired$go_direct)),
  all(paired$CM_screened == (paired$CM_direct & paired$go_direct)),
  all(programs$phase3_selected$program %in% programs$phase2$program[programs$phase2$go])
)

# If no programs advance, retain their denominators and return NA given-Go PoS.
never_go <- mcp
never_go$critical_value <- Inf

stopped <- simulate_programs(truth, never_go, n_sim = 10, seed = 805)
summary <- summarize_programs(stopped)

stopifnot(
  nrow(stopped$phase3_selected) == 0,
  !any(dplyr::filter(stopped$events, strategy == "P2_selected")$both),
  all(is.na(dplyr::filter(summary$probabilities, estimand == "given_go")$probability))
)

message("Core statistical regression checks passed.")

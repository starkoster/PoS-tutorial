#' Simulate Program 1: straight to Phase 3
#'
#' Program 1 has two Phase 3 trials, one EM and one CM, each with placebo,
#' 20 mg, and 200 mg, 233 subjects per arm, and subject SD 5 by default.
#'
#' @param draws Posterior Emax draws.
#' @param n_sim Number of simulations.
#' @param seed Random seed.
#' @param sigma Subject-level SD.
#' @param doses Dose vector including placebo.
#' @param n_per_arm Subjects per arm.
#' @param thresholds Positive placebo-adjusted MMD reduction thresholds.
#' @param alpha_p3 Two-sided alpha for active-placebo tests.
#' @return A list of raw and summarized simulation results.
#' @export
simulate_program1_direct_phase3 <- function(draws,
                                            n_sim = 5000,
                                            seed = 101,
                                            sigma = 5,
                                            doses = c(0, 20, 200),
                                            n_per_arm = 233,
                                            thresholds = c(2, 3, 4),
                                            alpha_p3 = 0.05) {
  required <- c("e0_em", "emax_em", "e0_cm", "emax_cm", "ed50")
  missing <- setdiff(required, names(draws))
  if (length(missing) > 0) {
    stop("Missing columns in draws: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  set.seed(seed)
  draw_id <- sample(seq_len(nrow(draws)), n_sim, replace = TRUE)
  rows <- vector("list", n_sim * 2 * (length(doses) - 1))
  idx <- 1

  for (s in seq_len(n_sim)) {
    draw <- draws[draw_id[s], , drop = FALSE]

    for (pop in c("EM", "CM")) {
      mu <- sapply(doses, function(d) mean_emax_pop(draw, pop, d))
      ybar <- stats::rnorm(length(doses), mean = mu, sd = sigma / sqrt(n_per_arm))
      y_pbo <- ybar[doses == 0]

      for (dose in doses[doses != 0]) {
        y_act <- ybar[doses == dose]
        tst <- active_vs_placebo_test(
          y_pbo = y_pbo,
          y_act = y_act,
          n_pbo = n_per_arm,
          n_act = n_per_arm,
          sigma = sigma,
          alpha_two_sided = alpha_p3
        )

        rows[[idx]] <- data.frame(
          sim = s,
          pop = pop,
          dose = dose,
          y_pbo = y_pbo,
          y_act = y_act,
          obs_effect = tst$obs_effect,
          true_effect = true_effect_pop(draw, pop, dose),
          p_two = tst$p_two,
          significant = tst$significant
        )
        idx <- idx + 1
      }
    }
  }

  simdat <- do.call(rbind, rows)

  for (thr in thresholds) {
    lab <- gsub("\\.", "_", as.character(thr))
    simdat[[paste0("obs_effect_ge_", lab)]] <- simdat$obs_effect >= thr
    simdat[[paste0("true_effect_ge_", lab)]] <- simdat$true_effect >= thr
    simdat[[paste0("sig_and_obs_effect_ge_", lab)]] <- simdat$significant & simdat$obs_effect >= thr
  }

  metric_cols <- c(
    "significant",
    grep("^obs_effect_ge_", names(simdat), value = TRUE),
    grep("^true_effect_ge_", names(simdat), value = TRUE),
    grep("^sig_and_obs_effect_ge_", names(simdat), value = TRUE)
  )

  by_pop_dose <- stats::aggregate(
    simdat[, metric_cols, drop = FALSE],
    by = list(pop = simdat$pop, dose = simdat$dose),
    FUN = mean
  )

  trial_success_list <- list()
  for (thr in thresholds) {
    lab <- gsub("\\.", "_", as.character(thr))
    col <- paste0("sig_and_obs_effect_ge_", lab)

    em_success <- base::tapply(simdat[[col]][simdat$pop == "EM"], simdat$sim[simdat$pop == "EM"], any)
    cm_success <- base::tapply(simdat[[col]][simdat$pop == "CM"], simdat$sim[simdat$pop == "CM"], any)

    trial_success_list[[lab]] <- data.frame(
      threshold = thr,
      em_trial_success = mean(em_success),
      cm_trial_success = mean(cm_success),
      program_success_both_trials = mean(em_success & cm_success)
    )
  }
  threshold_program_success <- do.call(rbind, trial_success_list)

  em_sig_any <- base::tapply(simdat$significant[simdat$pop == "EM"], simdat$sim[simdat$pop == "EM"], any)
  cm_sig_any <- base::tapply(simdat$significant[simdat$pop == "CM"], simdat$sim[simdat$pop == "CM"], any)

  significance_only <- data.frame(
    em_trial_any_dose_significant = mean(em_sig_any),
    cm_trial_any_dose_significant = mean(cm_sig_any),
    program_both_trials_any_dose_significant = mean(em_sig_any & cm_sig_any)
  )

  list(
    raw = simdat,
    by_pop_dose = by_pop_dose,
    significance_only = significance_only,
    threshold_program_success = threshold_program_success
  )
}

#' Simulate Program 2: Phase 2, then selected-dose Phase 3
#'
#' Phase 2 is a 50/50 EM/CM mixed trial with placebo, 20 mg, 100 mg, and
#' 300 mg, 140 subjects per arm. Go requires a significant Emax-shaped
#' dose-response test and at least one dose with observed effect >= 2 MMD.
#' If Go, the smallest such dose is selected for two Phase 3 trials, one EM
#' and one CM, each with placebo and selected dose and 174 subjects per arm.
#'
#' @param draws Posterior Emax draws.
#' @param n_sim Number of simulations.
#' @param seed Random seed.
#' @param sigma Subject-level SD.
#' @param p2_doses Phase 2 doses.
#' @param p2_n_per_arm Phase 2 subjects per arm.
#' @param p2_prop_em Phase 2 EM proportion.
#' @param p2_fixed_split Whether to force exactly `p2_prop_em` per arm.
#' @param p2_effect_go Required observed effect for dose selection.
#' @param p2_alpha_dr Phase 2 dose-response alpha.
#' @param p2_ed50_candidates Candidate ED50 values for dose-response test.
#' @param p2_n_null Null simulations for the MCP critical value.
#' @param p3_n_per_arm Phase 3 subjects per arm.
#' @param thresholds Positive placebo-adjusted MMD reduction thresholds.
#' @param alpha_p3 Two-sided Phase 3 alpha.
#' @return A list of raw and summarized simulation results.
#' @export
simulate_program2_phase2_then_phase3 <- function(draws,
                                                 n_sim = 5000,
                                                 seed = 202,
                                                 sigma = 5,
                                                 p2_doses = c(0, 20, 100, 300),
                                                 p2_n_per_arm = 140,
                                                 p2_prop_em = 0.50,
                                                 p2_fixed_split = TRUE,
                                                 p2_effect_go = 2,
                                                 p2_alpha_dr = 0.05,
                                                 p2_ed50_candidates = NULL,
                                                 p2_n_null = 200000,
                                                 p3_n_per_arm = 174,
                                                 thresholds = c(2, 3, 4),
                                                 alpha_p3 = 0.05) {
  required <- c("e0_em", "emax_em", "e0_cm", "emax_cm", "ed50")
  missing <- setdiff(required, names(draws))
  if (length(missing) > 0) {
    stop("Missing columns in draws: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (is.null(p2_ed50_candidates)) {
    p2_ed50_candidates <- as.numeric(stats::quantile(draws$ed50, probs = c(0.10, 0.50, 0.90)))
  }

  p2_test <- make_emax_mcp_test(
    doses = p2_doses,
    n_per_arm = p2_n_per_arm,
    sigma = sigma,
    ed50_candidates = p2_ed50_candidates,
    alpha = p2_alpha_dr,
    n_null = p2_n_null,
    seed = seed + 999
  )

  set.seed(seed)
  draw_id <- sample(seq_len(nrow(draws)), n_sim, replace = TRUE)
  phase2_rows <- vector("list", n_sim)
  phase3_rows <- vector("list", n_sim * 2)
  idx_p3 <- 1
  active_doses <- p2_doses[p2_doses != 0]

  for (s in seq_len(n_sim)) {
    draw <- draws[draw_id[s], , drop = FALSE]
    ybar_p2 <- numeric(length(p2_doses))

    for (a in seq_along(p2_doses)) {
      dose <- p2_doses[a]

      if (p2_fixed_split) {
        n_em <- round(p2_n_per_arm * p2_prop_em)
      } else {
        n_em <- stats::rbinom(1, size = p2_n_per_arm, prob = p2_prop_em)
      }
      n_cm <- p2_n_per_arm - n_em

      mu_em <- mean_emax_pop(draw, "EM", dose)
      mu_cm <- mean_emax_pop(draw, "CM", dose)
      mu_mix <- (n_em * mu_em + n_cm * mu_cm) / p2_n_per_arm

      ybar_p2[a] <- stats::rnorm(1, mean = mu_mix, sd = sigma / sqrt(p2_n_per_arm))
    }

    p2_dr <- p2_test(ybar_p2)
    y_pbo <- ybar_p2[p2_doses == 0]
    p2_effects <- sapply(active_doses, function(dose) y_pbo - ybar_p2[p2_doses == dose])
    names(p2_effects) <- paste0("effect_", active_doses)

    doses_meeting_go <- active_doses[p2_effects >= p2_effect_go]
    any_dose_effect_go <- length(doses_meeting_go) > 0
    go <- p2_dr$significant && any_dose_effect_go
    selected_dose <- if (go) min(doses_meeting_go) else NA_real_

    phase2_rows[[s]] <- data.frame(
      sim = s,
      dose_response_significant = p2_dr$significant,
      dose_response_p_adj = p2_dr$p_adj,
      dose_response_tmax = p2_dr$tmax,
      best_ed50_candidate = p2_dr$best_ed50_candidate,
      any_dose_effect_ge_2 = any_dose_effect_go,
      go = go,
      selected_dose = selected_dose,
      effect_20 = if ("effect_20" %in% names(p2_effects)) p2_effects[["effect_20"]] else NA_real_,
      effect_100 = if ("effect_100" %in% names(p2_effects)) p2_effects[["effect_100"]] else NA_real_,
      effect_300 = if ("effect_300" %in% names(p2_effects)) p2_effects[["effect_300"]] else NA_real_
    )

    if (go) {
      for (pop in c("EM", "CM")) {
        mu_pbo <- mean_emax_pop(draw, pop, 0)
        mu_act <- mean_emax_pop(draw, pop, selected_dose)

        y_pbo_p3 <- stats::rnorm(1, mean = mu_pbo, sd = sigma / sqrt(p3_n_per_arm))
        y_act_p3 <- stats::rnorm(1, mean = mu_act, sd = sigma / sqrt(p3_n_per_arm))

        tst <- active_vs_placebo_test(
          y_pbo = y_pbo_p3,
          y_act = y_act_p3,
          n_pbo = p3_n_per_arm,
          n_act = p3_n_per_arm,
          sigma = sigma,
          alpha_two_sided = alpha_p3
        )

        phase3_rows[[idx_p3]] <- data.frame(
          sim = s,
          pop = pop,
          selected_dose = selected_dose,
          y_pbo = y_pbo_p3,
          y_act = y_act_p3,
          obs_effect = tst$obs_effect,
          true_effect = true_effect_pop(draw, pop, selected_dose),
          p_two = tst$p_two,
          significant = tst$significant
        )
        idx_p3 <- idx_p3 + 1
      }
    }
  }

  phase2 <- do.call(rbind, phase2_rows)
  phase3 <- if (idx_p3 > 1) do.call(rbind, phase3_rows[seq_len(idx_p3 - 1)]) else data.frame()

  if (nrow(phase3) > 0) {
    for (thr in thresholds) {
      lab <- gsub("\\.", "_", as.character(thr))
      phase3[[paste0("obs_effect_ge_", lab)]] <- phase3$obs_effect >= thr
      phase3[[paste0("true_effect_ge_", lab)]] <- phase3$true_effect >= thr
      phase3[[paste0("sig_and_obs_effect_ge_", lab)]] <- phase3$significant & phase3$obs_effect >= thr
    }
  }

  phase2_summary <- data.frame(
    metric = c("dose_response_significant", "any_dose_effect_ge_2", "phase2_go"),
    probability = c(
      mean(phase2$dose_response_significant),
      mean(phase2$any_dose_effect_ge_2),
      mean(phase2$go)
    )
  )

  selected_label <- ifelse(phase2$go, as.character(phase2$selected_dose), "No Go")
  selection_summary <- as.data.frame(prop.table(table(selected_label)))
  names(selection_summary) <- c("selection", "probability")

  if (nrow(phase3) > 0) {
    metric_cols <- c(
      "significant",
      grep("^obs_effect_ge_", names(phase3), value = TRUE),
      grep("^true_effect_ge_", names(phase3), value = TRUE),
      grep("^sig_and_obs_effect_ge_", names(phase3), value = TRUE)
    )

    phase3_by_pop_conditional_on_go <- stats::aggregate(
      phase3[, metric_cols, drop = FALSE],
      by = list(pop = phase3$pop),
      FUN = mean
    )
  } else {
    phase3_by_pop_conditional_on_go <- data.frame()
  }

  go_vec <- phase2$go
  program_rows <- list()

  em_sig <- rep(FALSE, n_sim)
  cm_sig <- rep(FALSE, n_sim)

  if (nrow(phase3) > 0) {
    em_rows <- phase3$pop == "EM"
    cm_rows <- phase3$pop == "CM"
    em_sig[phase3$sim[em_rows]] <- phase3$significant[em_rows]
    cm_sig[phase3$sim[cm_rows]] <- phase3$significant[cm_rows]
  }

  program_rows[["sig_only"]] <- data.frame(
    endpoint = "significance_only",
    threshold = NA_real_,
    p_phase2_go = mean(go_vec),
    p_p3_both_success_given_go = ifelse(sum(go_vec) > 0, mean((em_sig & cm_sig)[go_vec]), NA_real_),
    p_overall_program_success = mean(go_vec & em_sig & cm_sig)
  )

  for (thr in thresholds) {
    lab <- gsub("\\.", "_", as.character(thr))
    col <- paste0("sig_and_obs_effect_ge_", lab)
    em_success <- rep(FALSE, n_sim)
    cm_success <- rep(FALSE, n_sim)

    if (nrow(phase3) > 0) {
      em_rows <- phase3$pop == "EM"
      cm_rows <- phase3$pop == "CM"
      em_success[phase3$sim[em_rows]] <- phase3[[col]][em_rows]
      cm_success[phase3$sim[cm_rows]] <- phase3[[col]][cm_rows]
    }

    program_rows[[paste0("thr_", lab)]] <- data.frame(
      endpoint = "significant_and_observed_effect_threshold",
      threshold = thr,
      p_phase2_go = mean(go_vec),
      p_p3_both_success_given_go = ifelse(sum(go_vec) > 0, mean((em_success & cm_success)[go_vec]), NA_real_),
      p_overall_program_success = mean(go_vec & em_success & cm_success)
    )
  }

  program_success <- do.call(rbind, program_rows)

  list(
    phase2 = phase2,
    phase3 = phase3,
    phase2_summary = phase2_summary,
    selection_summary = selection_summary,
    phase3_by_pop_conditional_on_go = phase3_by_pop_conditional_on_go,
    program_success = program_success
  )
}

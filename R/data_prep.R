#' Find the first existing file among candidate paths
#'
#' @param paths Character vector of candidate paths.
#' @return A single path.
#' @keywords internal
first_existing_file <- function(paths) {
  paths <- paths[nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) {
    stop("None of the candidate files exist: ", paste(paths, collapse = ", "), call. = FALSE)
  }
  normalizePath(hit[[1]], mustWork = TRUE)
}

#' Load the included migraine trial CSV
#'
#' The default CSV is copied into both `data/hist_trials.csv` and
#' `inst/extdata/hist_trials.csv`. The `inst/extdata` copy is available after
#' installation via `system.file()`.
#'
#' @param path Optional path to a CSV file. If `NULL`, the package default is used.
#' @return A cleaned `data.frame` with a `w_em` column and, when available,
#' baseline aliases `base_y`, `base_se`, and `base_sd` added.
#' @export
load_trial_csv <- function(path = NULL) {
  if (is.null(path)) {
    path <- first_existing_file(c(
      system.file("extdata", "hist_trials.csv", package = "MigraineEmaxSim", mustWork = FALSE),
      file.path("data", "hist_trials.csv"),
      file.path("inst", "extdata", "hist_trials.csv")
    ))
  }

  dat <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  clean_trial_data(dat)
}

#' Clean and validate trial data
#'
#' Expected columns are `nct`, `title`, `CM`, `arm`, `dose`, `n`, `chg`, and
#' `chg_se`. The input CSV has `CM` as the proportion of chronic migraine
#' patients; this function adds `w_em = 1 - CM`.
#'
#' Baseline columns are optional but strongly recommended. If `base_mmd` and
#' `base_mmd_sd` are present, the function creates `base_y = base_mmd`,
#' `base_sd = base_mmd_sd`, and `base_se = base_mmd_sd / sqrt(n)`. If a
#' baseline standard error is already supplied as `base_mmd_se` or `base_se`,
#' that value is used instead.
#'
#' @param dat Trial-level arm data.
#' @return A cleaned `data.frame`.
#' @export
clean_trial_data <- function(dat) {
  names(dat) <- trimws(names(dat))

  drop_col <- names(dat) %in% c("", "X", "X1") |
    grepl("^Unnamed", names(dat), ignore.case = TRUE)
  dat <- dat[, !drop_col, drop = FALSE]

  required <- c("nct", "title", "CM", "arm", "dose", "n", "chg", "chg_se")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  dat$nct <- as.character(dat$nct)
  dat$title <- as.character(dat$title)
  dat$arm <- tolower(trimws(as.character(dat$arm)))

  numeric_cols <- intersect(
    c(
      "CM", "dose", "n", "chg", "chg_se", "chg_sd",
      "base_mmd", "base_mmd_sd", "base_mmd_se",
      "base_y", "base_sd", "base_se"
    ),
    names(dat)
  )
  for (nm in numeric_cols) {
    dat[[nm]] <- as.numeric(dat[[nm]])
  }

  if (!all(dat$arm %in% c("placebo", "active"))) {
    stop("Column 'arm' must contain only 'placebo' or 'active'.", call. = FALSE)
  }
  if (any(!is.finite(dat$CM) | dat$CM < 0 | dat$CM > 1)) {
    stop("Column 'CM' must be a proportion between 0 and 1.", call. = FALSE)
  }
  if (any(!is.finite(dat$n) | dat$n <= 0)) {
    stop("Column 'n' must contain positive arm sizes.", call. = FALSE)
  }
  if (any(!is.finite(dat$chg))) {
    stop("Column 'chg' must contain finite change-from-baseline values.", call. = FALSE)
  }

  if (any(!is.finite(dat$chg_se) | dat$chg_se <= 0)) {
    if ("chg_sd" %in% names(dat)) {
      idx <- !is.finite(dat$chg_se) | dat$chg_se <= 0
      dat$chg_se[idx] <- dat$chg_sd[idx] / sqrt(dat$n[idx])
    }
  }
  if (any(!is.finite(dat$chg_se) | dat$chg_se <= 0)) {
    stop("Column 'chg_se' must contain positive standard errors.", call. = FALSE)
  }

  # Standardize baseline naming. The Stan models use baseline means and
  # standard errors as an additional likelihood contribution.
  if (!"base_y" %in% names(dat) && "base_mmd" %in% names(dat)) {
    dat$base_y <- dat$base_mmd
  }
  if (!"base_sd" %in% names(dat) && "base_mmd_sd" %in% names(dat)) {
    dat$base_sd <- dat$base_mmd_sd
  }
  if (!"base_se" %in% names(dat)) {
    if ("base_mmd_se" %in% names(dat)) {
      dat$base_se <- dat$base_mmd_se
    } else if ("base_sd" %in% names(dat)) {
      dat$base_se <- dat$base_sd / sqrt(dat$n)
    }
  }

  if ("base_y" %in% names(dat)) {
    if (any(is.finite(dat$base_y) & (dat$base_y < 0 | dat$base_y > 30))) {
      warning("Some baseline MMD means are outside the 0 to 30 range.")
    }
  }
  if ("base_se" %in% names(dat)) {
    bad_base_se <- is.finite(dat$base_y) & (!is.finite(dat$base_se) | dat$base_se <= 0)
    if (any(bad_base_se)) {
      warning("Some baseline MMD rows have missing/non-positive baseline SE and will be ignored in baseline likelihood.")
    }
  }

  dat$w_em <- 1 - dat$CM
  dat
}

#' Split trial data into historical and HOPE data
#'
#' @param dat Clean or raw trial data.
#' @param hope_title Title identifying the new trial. Defaults to `"HOPE"`.
#' @param hope_nct Optional NCT identifier identifying the new trial.
#' @return A list with `all`, `historical`, and `hope` data frames.
#' @export
prepare_trial_data <- function(dat, hope_title = "HOPE", hope_nct = NULL) {
  dat <- clean_trial_data(dat)

  if (!is.null(hope_nct)) {
    is_hope <- dat$nct == hope_nct
  } else {
    is_hope <- tolower(dat$title) == tolower(hope_title)
  }

  if (!any(is_hope)) {
    stop("No rows identified as the HOPE/new trial.", call. = FALSE)
  }
  if (!any(!is_hope)) {
    stop("No historical rows found after excluding HOPE/new trial.", call. = FALSE)
  }

  list(
    all = dat,
    historical = dat[!is_hope, , drop = FALSE],
    hope = dat[is_hope, , drop = FALSE]
  )
}

#' Build historical change data for the meta-analysis
#'
#' The historical numerical `dose` column is not used. By default historical
#' active rows are excluded because those active doses are from other drugs.
#' If `use_historical_active = TRUE`, active rows are treated only as a
#' reference active effect through `active_ref = 1`; they are still not treated
#' as dose-response information.
#'
#' @param historical Historical arm-level data.
#' @param use_historical_active Whether to include historical active rows.
#' @param active_borrow_weight Likelihood weight for included historical active rows.
#' @return A data frame ready for `historical_meta_stan_data()`.
#' @export
make_historical_change_data <- function(historical,
                                        use_historical_active = FALSE,
                                        active_borrow_weight = 0) {
  historical <- clean_trial_data(historical)

  if (!isTRUE(use_historical_active)) {
    active_borrow_weight <- 0
  }
  if (!is.finite(active_borrow_weight) || active_borrow_weight < 0 || active_borrow_weight > 1) {
    stop("active_borrow_weight must be between 0 and 1.", call. = FALSE)
  }

  out <- historical
  out$active_ref <- ifelse(out$arm == "active", 1, 0)
  out$weight <- ifelse(out$active_ref == 1, active_borrow_weight, 1)
  out <- out[out$weight > 0, , drop = FALSE]

  if (nrow(out) == 0) {
    stop("No historical rows remain after applying active-row borrowing settings.", call. = FALSE)
  }
  if (!any(out$active_ref == 0)) {
    stop("Historical meta-analysis must include at least one placebo row.", call. = FALSE)
  }

  nct_levels <- unique(historical$nct)
  out$study_id <- match(out$nct, nct_levels)
  out <- out[order(out$study_id, out$active_ref, out$CM), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Build baseline MMD data
#'
#' Baseline MMD means and SDs are used as direct evidence about the baseline
#' monthly migraine day parameters. For rows with `base_mmd_sd`, the standard
#' error is computed as `base_mmd_sd / sqrt(n)`. These baseline rows are used
#' regardless of whether the corresponding active change row is borrowed,
#' because baseline disease burden is not a drug-specific dose-response effect.
#'
#' @param dat Trial data with baseline columns.
#' @return Baseline data frame, possibly empty.
#' @export
make_baseline_data <- function(dat) {
  dat <- clean_trial_data(dat)

  if (!all(c("base_y", "base_se") %in% names(dat))) {
    return(data.frame(
      nct = character(0),
      title = character(0),
      study_id = integer(0),
      base_y = numeric(0),
      base_se = numeric(0),
      w_em = numeric(0)
    ))
  }

  nct_levels <- unique(dat$nct)
  keep <- is.finite(dat$base_y) & is.finite(dat$base_se) & dat$base_se > 0
  out <- dat[keep, c("nct", "title", "base_y", "base_se", "w_em"), drop = FALSE]
  out$study_id <- match(out$nct, nct_levels)
  rownames(out) <- NULL
  out
}

#' Build HOPE Emax data
#'
#' @param hope HOPE/new-trial arm-level data.
#' @return A data frame sorted by dose and ready for `historical_emax_stan_data()`.
#' @export
make_hope_emax_data <- function(hope) {
  hope <- clean_trial_data(hope)
  out <- hope[order(hope$dose), , drop = FALSE]
  if (!any(out$dose == 0)) {
    stop("HOPE data must include a placebo arm with dose 0.", call. = FALSE)
  }
  if (length(unique(out$dose)) < 3) {
    warning("HOPE data has fewer than three unique doses; ED50 may be weakly identified.")
  }
  rownames(out) <- NULL
  out
}

#' Build Stan data for the historical meta-analysis
#'
#' @param change_data Output of `make_historical_change_data()`.
#' @param base_data Output of `make_baseline_data()` or `NULL`.
#' @param priors Output of `default_eta_priors()`.
#' @param lkj_shape LKJ shape parameter for the between-study correlation prior.
#' @return A list for Stan.
#' @export
historical_meta_stan_data <- function(change_data,
                                      base_data = NULL,
                                      priors,
                                      lkj_shape = 2) {
  if (is.null(base_data)) {
    base_data <- data.frame(
      study_id = integer(0),
      base_y = numeric(0),
      base_se = numeric(0),
      w_em = numeric(0)
    )
  }

  list(
    N_change = nrow(change_data),
    J = max(change_data$study_id),
    study_change = as.integer(change_data$study_id),
    y_change = as.numeric(change_data$chg),
    se_change = as.numeric(change_data$chg_se),
    w_em = as.numeric(change_data$w_em),
    active_ref = as.numeric(change_data$active_ref),
    weight_change = as.numeric(change_data$weight),
    N_base = nrow(base_data),
    study_base = as.integer(base_data$study_id),
    y_base = as.numeric(base_data$base_y),
    se_base = as.numeric(base_data$base_se),
    w_em_base = as.numeric(base_data$w_em),
    mu_eta_prior = as.numeric(priors$mu_eta_prior),
    sd_mu_eta_prior = as.numeric(priors$sd_mu_eta_prior),
    tau_eta_scale = as.numeric(priors$tau_eta_scale),
    lkj_shape = lkj_shape
  )
}

#' Build Stan data for the HOPE Emax model
#'
#' @param hope_data Output of `make_hope_emax_data()`.
#' @param eta_prior List with `mu` and `Sigma` from `extract_eta_prior()`.
#' @param ed50_prior List from `make_ed50_prior()`.
#' @param dose_grid Optional prediction dose grid.
#' @param base_data Optional baseline data. If `NULL`, baseline rows are
#' extracted from `hope_data` when available.
#' @return A list for Stan.
#' @export
historical_emax_stan_data <- function(hope_data,
                                      eta_prior,
                                      ed50_prior,
                                      dose_grid = NULL,
                                      base_data = NULL) {
  hope_data <- make_hope_emax_data(hope_data)

  if (is.null(base_data)) {
    base_data <- make_baseline_data(hope_data)
  }
  if (is.null(dose_grid)) {
    dose_grid <- seq(0, max(hope_data$dose) * 1.2, length.out = 100)
  }

  list(
    N = nrow(hope_data),
    y = as.numeric(hope_data$chg),
    se = as.numeric(hope_data$chg_se),
    dose = as.numeric(hope_data$dose),
    w_em = as.numeric(hope_data$w_em),
    N_base = nrow(base_data),
    y_base = as.numeric(base_data$base_y),
    se_base = as.numeric(base_data$base_se),
    w_em_base = as.numeric(base_data$w_em),
    eta_prior_mu = as.numeric(eta_prior$mu),
    eta_prior_Sigma = as.matrix(eta_prior$Sigma),
    ed50_prior_logmean = ed50_prior$logmean,
    ed50_prior_logsd = ed50_prior$logsd,
    N_pred = length(dose_grid),
    dose_pred = as.numeric(dose_grid)
  )
}

#' Backwards-compatible split helper
#'
#' @param dat Clean or raw trial data.
#' @param new_trial_title Title identifying the new trial.
#' @param new_trial_nct Optional NCT identifier identifying the new trial.
#' @return List with `historical`, `new_trial`, and `all` data frames.
#' @export
split_trial_data <- function(dat, new_trial_title = "HOPE", new_trial_nct = NULL) {
  x <- prepare_trial_data(dat, hope_title = new_trial_title, hope_nct = new_trial_nct)
  list(all = x$all, historical = x$historical, new_trial = x$hope)
}

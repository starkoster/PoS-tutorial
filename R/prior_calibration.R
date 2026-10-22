default_targets <- function() {

  # For negative effects, the fifth percentile is the largest reduction.
  tibble::tribble(
    ~quantity,       ~q05, ~q50, ~q95,
    "baseline_em",   6.0, 10.0, 14.5,
    "baseline_cm",  11.0, 15.0, 21.0,
    "placebo_em",   -2.5, -0.5, -0.1,
    "placebo_cm",   -4.0, -1.0, -0.25,
    "emax_em",      -2.5, -0.5, -0.1,
    "emax_cm",      -4.0, -1.0, -0.25
  )
}


validate_targets <- function(targets) {

  # Fix the order so clinical quantities always match the latent coordinates.
  required <- default_targets()$quantity

  stopifnot(all(c("quantity", "q05", "q50", "q95") %in% names(targets)))

  if (anyDuplicated(targets$quantity) || !setequal(targets$quantity, required)) {
    stop("Supply exactly the six quantities shown in default_targets().")
  }

  targets <- targets[match(required, targets$quantity), ]
  quantiles <- as.matrix(targets[, c("q05", "q50", "q95")])

  if (any(!is.finite(quantiles)) ||
      any(quantiles[, 1] >= quantiles[, 2] | quantiles[, 2] >= quantiles[, 3])) {
    stop("Each target must have finite q05 < q50 < q95.")
  }

  if (any(quantiles[1:2, ] <= 0 | quantiles[1:2, ] >= 30) ||
      any(quantiles[3:6, ] >= 0)) {
    stop("Baselines must be between 0 and 30; effect targets must be negative.")
  }

  if (any(quantiles[2, ] < quantiles[1, ]) ||
      any(quantiles[4, ] > quantiles[3, ]) ||
      any(quantiles[6, ] > quantiles[5, ])) {
    stop("The target quantiles conflict with the model's EM/CM ordering.")
  }

  targets
}


lkj_cholesky <- function(shape = 2, dimension = 6L) {

  # Canonical partial correlations give an exact LKJ Cholesky draw.
  # This small helper keeps calibration consistent with Stan's LKJ prior.
  correlation_factor <- diag(dimension)

  for (row in 2:dimension) {

    remaining_length <- 1

    for (column in seq_len(row - 1L)) {

      beta_shape <- shape + (dimension - column - 1) / 2
      partial_correlation <- 2 * stats::rbeta(1, beta_shape, beta_shape) - 1

      correlation_factor[row, column] <- partial_correlation * remaining_length
      remaining_length <- remaining_length * sqrt(1 - partial_correlation^2)
    }

    correlation_factor[row, row] <- remaining_length
  }

  correlation_factor
}


prior_random_numbers <- function(n, seed, lkj_shape) {

  # Generate once and reuse during optimization: the objective then has no
  # extra Monte Carlo noise when only the candidate prior changes.
  set.seed(seed)

  mean_z <- matrix(stats::rnorm(n * 6), ncol = 6)
  scale_z <- abs(matrix(stats::rnorm(n * 6), ncol = 6))
  study_z <- matrix(stats::rnorm(n * 6), ncol = 6)

  for (row in seq_len(n)) {
    study_z[row, ] <- lkj_cholesky(lkj_shape) %*% study_z[row, ]
  }

  list(mean_z = mean_z, scale_z = scale_z, study_z = study_z)
}


prior_eta_draws <- function(prior, random) {

  # A new study varies through uncertainty in mu and heterogeneity around mu.
  mean_part <- sweep(random$mean_z, 2, prior$s, "*")
  study_part <- sweep(random$scale_z * random$study_z, 2, prior$q, "*")

  sweep(mean_part + study_part, 2, prior$m, "+")
}


simulate_prior <- function(prior, n = 50000L, seed = 4202) {

  # Draw clinical mean parameters of new studies, without patient-level noise.
  random <- prior_random_numbers(n, seed, prior$lkj_shape)

  clinical_from_eta(prior_eta_draws(prior, random))
}


calibrate_priors <- function(
    targets = default_targets(),
    heterogeneity_fraction = c(4 / 13, 4 / 13, rep(0.36, 4)),
    lkj_shape = 2,
    n_train = 20000L,
    n_validate = 100000L,
    seed = 4210,
    maxit = 3500L,
    n_starts = 3L) {

  targets <- validate_targets(targets)

  stopifnot(
    length(heterogeneity_fraction) == 6,
    all(heterogeneity_fraction > 0 & heterogeneity_fraction < 1),
    lkj_shape > 0,
    n_train >= 1000,
    n_validate >= 1000,
    n_starts >= 1
  )

  # Clinical quantiles identify total latent variation. Its split between
  # mean uncertainty and heterogeneity is a separate user-specified judgment.
  h <- heterogeneity_fraction
  target_matrix <- as.matrix(targets[, c("q05", "q50", "q95")])

  reference <- stats::setNames(targets$q50, targets$quantity)
  starting_center <- eta_from_reference(reference)

  random <- prior_random_numbers(n_train, seed, lkj_shape)

  unpack <- function(parameters) {

    # Optimize log SDs so every candidate scale stays positive.
    total_sd <- exp(parameters[7:12])

    list(
      m = parameters[1:6],
      s = total_sd * sqrt(1 - h),
      q = total_sd * sqrt(h),
      lkj_shape = lkj_shape
    )
  }

  objective <- function(parameters) {

    # Avoid extreme trial values that make the numerical search unstable.
    if (any(abs(parameters[1:6]) > 12) || any(abs(parameters[7:12]) > 4)) {
      return(1e6)
    }

    candidate <- unpack(parameters)
    draws <- clinical_from_eta(prior_eta_draws(candidate, random))

    achieved <- clinical_quantiles(draws)[1:6, c("q05", "q50", "q95")]

    # Relative quantile errors receive equal weight across clinical quantities.
    sum(log(abs(as.matrix(achieved)) / abs(target_matrix))^2)
  }

  # Several nearby starts reduce dependence on one local search trajectory.
  fits <- purrr::map(seq_len(n_starts), function(start) {

    multiplier <- c(1, 0.8, 1.2)[1 + (start - 1) %% 3]
    starting_sd <- c(0.42, 0.56, 1.05, 1.45, 0.65, 1.25) * multiplier

    stats::optim(
      par = c(starting_center, log(starting_sd)),
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = maxit, reltol = 1e-9)
    )
  })

  best <- fits[[which.min(purrr::map_dbl(fits, "value"))]]
  prior <- unpack(best$par)

  # Validate with an independent sample, rather than the optimization sample.
  draws <- simulate_prior(prior, n_validate, seed + 1L)
  achieved <- clinical_quantiles(draws)

  comparison <- dplyr::left_join(
    targets,
    achieved,
    by = "quantity",
    suffix = c("_target", "_achieved")
  )

  if (best$convergence != 0) {
    warning("Calibration search reached its iteration limit; review achieved quantiles.")
  }

  list(
    prior = prior,
    targets = targets,
    comparison = comparison,
    draws = draws,
    optimizer = best,
    heterogeneity_fraction = h
  )
}

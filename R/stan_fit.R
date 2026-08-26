#' Resolve a package Stan file path
#'
#' @param file Stan file name, for example `"historical_meta.stan"`.
#' @param root Optional package root.
#' @return Existing path to the Stan file.
#' @export
migraine_stan_file <- function(file, root = NULL) {
  candidates <- character(0)
  if (!is.null(root)) {
    candidates <- c(
      file.path(root, "stan", file),
      file.path(root, "inst", "stan", file)
    )
  }

  candidates <- c(
    candidates,
    system.file("stan", file, package = "MigraineEmaxSim", mustWork = FALSE),
    file.path("stan", file),
    file.path("inst", "stan", file)
  )

  first_existing_file(candidates)
}

#' Fit the constrained historical meta-analysis
#'
#' Historical numerical dose values are ignored. Historical active rows are
#' excluded by default via `use_historical_active = FALSE`.
#'
#' @param historical Historical arm-level data.
#' @param priors Output of `default_eta_priors()`.
#' @param use_historical_active Whether to include historical active rows as
#' non-dose reference active effects.
#' @param active_borrow_weight Likelihood weight for historical active rows.
#' @param base_data Optional baseline data.
#' @param stan_file Optional Stan file path.
#' @param chains Number of chains.
#' @param parallel_chains Number of parallel chains.
#' @param iter_warmup Warmup iterations.
#' @param iter_sampling Sampling iterations.
#' @param seed Random seed.
#' @param adapt_delta Stan adapt_delta.
#' @param metric Stan metric.
#' @param ... Passed to `cmdstanr::CmdStanModel$sample()`.
#' @return A `CmdStanMCMC` fit.
#' @export
fit_historical_meta <- function(historical,
                                priors = default_eta_priors(),
                                use_historical_active = FALSE,
                                active_borrow_weight = 0,
                                base_data = NULL,
                                stan_file = NULL,
                                chains = 4,
                                parallel_chains = chains,
                                iter_warmup = 2000,
                                iter_sampling = 4000,
                                seed = 123,
                                adapt_delta = 0.97,
                                metric = "dense_e",
                                ...) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' is required to fit Stan models.", call. = FALSE)
  }

  change_data <- make_historical_change_data(
    historical = historical,
    use_historical_active = use_historical_active,
    active_borrow_weight = active_borrow_weight
  )

  if (is.null(base_data)) {
    base_data <- make_baseline_data(historical)
  }

  stan_data <- historical_meta_stan_data(
    change_data = change_data,
    base_data = base_data,
    priors = priors
  )

  if (is.null(stan_file)) {
    stan_file <- migraine_stan_file("historical_meta.stan")
  }

  mod <- cmdstanr::cmdstan_model(stan_file)
  mod$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    adapt_delta = adapt_delta,
    metric = metric,
    ...
  )
}

#' Extract a Gaussian eta prior from the historical posterior predictive draws
#'
#' @param fit_hist Historical meta-analysis fit.
#' @param variance_inflation Optional covariance inflation factor.
#' @param ridge Small diagonal ridge for numerical Cholesky stability.
#' @return A list with `mu`, `Sigma`, and raw `draws`.
#' @export
extract_eta_prior <- function(fit_hist, variance_inflation = 1, ridge = 1e-8) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required.", call. = FALSE)
  }
  if (!is.finite(variance_inflation) || variance_inflation <= 0) {
    stop("variance_inflation must be positive.", call. = FALSE)
  }

  eta_names <- paste0("eta_new[", 1:6, "]")
  mat <- posterior::as_draws_matrix(fit_hist$draws(eta_names))
  draws <- as.matrix(mat[, eta_names, drop = FALSE])
  colnames(draws) <- paste0("eta", 1:6)

  mu <- colMeans(draws)
  Sigma <- stats::cov(draws) * variance_inflation + diag(ridge, 6)

  list(mu = mu, Sigma = Sigma, draws = draws)
}

#' Fit the HOPE ED50 Emax model using the historical eta prior
#'
#' @param hope HOPE/new-trial arm-level data.
#' @param eta_prior Output of `extract_eta_prior()`.
#' @param ed50_prior Output of `make_ed50_prior()`.
#' @param dose_grid Optional prediction grid.
#' @param stan_file Optional Stan file path.
#' @param chains Number of chains.
#' @param parallel_chains Number of parallel chains.
#' @param iter_warmup Warmup iterations.
#' @param iter_sampling Sampling iterations.
#' @param seed Random seed.
#' @param adapt_delta Stan adapt_delta.
#' @param metric Stan metric.
#' @param ... Passed to `cmdstanr::CmdStanModel$sample()`.
#' @return A `CmdStanMCMC` fit.
#' @export
fit_historical_emax <- function(hope,
                                eta_prior,
                                ed50_prior,
                                dose_grid = NULL,
                                stan_file = NULL,
                                chains = 4,
                                parallel_chains = chains,
                                iter_warmup = 1000,
                                iter_sampling = 3000,
                                seed = 456,
                                adapt_delta = 0.97,
                                metric = "dense_e",
                                ...) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' is required to fit Stan models.", call. = FALSE)
  }

  hope_data <- make_hope_emax_data(hope)
  hope_base_data <- make_baseline_data(hope_data)
  stan_data <- historical_emax_stan_data(
    hope_data = hope_data,
    eta_prior = eta_prior,
    ed50_prior = ed50_prior,
    dose_grid = dose_grid,
    base_data = hope_base_data
  )

  if (is.null(stan_file)) {
    stan_file <- migraine_stan_file("historical_emax.stan")
  }

  mod <- cmdstanr::cmdstan_model(stan_file)
  mod$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    adapt_delta = adapt_delta,
    metric = metric,
    ...
  )
}

#' Extract Emax posterior draws for simulation
#'
#' @param fit_emax HOPE Emax fit.
#' @return Data frame with `e0_em`, `emax_em`, `e0_cm`, `emax_cm`, and `ed50`.
#' @export
extract_emax_draws <- function(fit_emax) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required.", call. = FALSE)
  }
  nm <- c("e0_em", "emax_em", "e0_cm", "emax_cm", "ed50")
  out <- as.data.frame(posterior::as_draws_df(fit_emax$draws(nm)))
  out[, nm, drop = FALSE]
}

#' Alias for locating a Stan file
#'
#' @param filename Stan file name.
#' @return Existing path to the Stan file.
#' @export
stan_file_path <- function(filename) {
  migraine_stan_file(filename)
}

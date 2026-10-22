fit_stan <- function(
    file,
    data,
    output_dir,
    seed,
    chains = 4L,
    warmup = 4000L,
    sampling = 4000L,
    adapt_delta = 0.995,
    max_treedepth = 13L) {

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Install cmdstanr and CmdStan as described in README.md.")
  }

  # Compile a writable copy, rather than writing into an installed package.
  model_dir <- file.path(output_dir, "model")
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  stan_dir <- system.file("stan", package = "MigraineEmaxSim", mustWork = TRUE)
  file.copy(list.files(stan_dir, full.names = TRUE), model_dir, overwrite = TRUE)

  model <- cmdstanr::cmdstan_model(file.path(model_dir, file))

  # CmdStanR's own diagnostic warnings remain enabled. There is no custom
  # diagnostic gate in the core workflow; optional review is a separate script.
  fit <- model$sample(
    data = data,
    seed = seed,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = warmup,
    iter_sampling = sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    output_dir = output_dir,
    refresh = 500
  )

  # The fit and original Stan chain files remain available for later review.
  fit$save_object(file = file.path(output_dir, "fit.rds"))

  fit
}


draw_matrix <- function(fit, variables) {

  # Explicit column order prevents parameter labels depending on extraction order.
  draws <- posterior::as_draws_matrix(fit$draws(variables = variables))

  # Return a BASE numeric matrix. Keeping draws_matrix here propagates its
  # class through the clinical transformation and creates matrix-valued columns.
  draws <- as_plain_matrix(draws)

  draws[, variables, drop = FALSE]
}

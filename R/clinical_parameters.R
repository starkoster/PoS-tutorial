# The same transformation is used in R and in inst/stan/clinical_functions.stan.
# Baseline MMDs are positive. Placebo, active and drug effects are negative.

as_plain_matrix <- function(x) {

  # as.matrix() can retain posterior's draws_matrix subclass. Reconstruct the
  # numeric matrix explicitly, preserving dimensions, labels and column order.
  x <- as.matrix(x)

  matrix(
    as.numeric(x),
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = dimnames(x)
  )
}


as_numeric_draws <- function(draws) {

  # A clinical draw table has one scalar per quantity and posterior row.
  # Older files may store each quantity as an n-by-1 posterior draws_matrix.
  stopifnot(is.data.frame(draws))

  columns <- purrr::imap(draws, function(values, quantity_name) {

    # Reject genuine multicolumn matrices and list columns rather than flatten
    # them silently, which could change the number or pairing of posterior rows.
    dimensions <- dim(values)
    scalar_column <- is.null(dimensions) ||
      identical(dimensions, c(nrow(draws), 1L))

    if (!is.numeric(values) || !scalar_column || length(values) != nrow(draws)) {
      stop("Expected one numeric value per draw in column: ", quantity_name,
           call. = FALSE)
    }

    # as.numeric() removes dimensions, draw classes and other container metadata.
    as.numeric(values)
  })

  # Rebuild the table so no class from the original draw container survives.
  tibble::as_tibble(columns, .name_repair = "check_unique")
}


clinical_from_eta <- function(eta) {

  # Each row represents one set of six unconstrained latent parameters.
  if (is.null(dim(eta))) eta <- matrix(eta, nrow = 1)

  # Strip draw-container metadata before arithmetic so every clinical column
  # becomes an ordinary numeric vector, including for Stan posterior inputs.
  eta <- as_plain_matrix(eta)
  stopifnot(ncol(eta) == 6L, all(is.finite(eta)))

  fraction <- stats::plogis(eta)

  # CM has a larger baseline burden, bounded above by 30 MMDs.
  baseline_em <- 30 * fraction[, 1]
  baseline_cm <- baseline_em + (30 - baseline_em) * fraction[, 2]

  # Total active reductions cannot exceed the corresponding baseline burden.
  active_em <- -baseline_em * fraction[, 3]
  active_cm <- active_em - (baseline_cm + active_em) * fraction[, 4]

  # Placebo reductions are smaller than the total active reductions.
  placebo_em <- active_em * fraction[, 5]
  placebo_cm <- placebo_em + (active_cm - active_em) * fraction[, 6]

  tibble::tibble(
    baseline_em,
    baseline_cm,
    placebo_em,
    placebo_cm,
    emax_em = active_em - placebo_em,
    emax_cm = active_cm - placebo_cm,
    active_em,
    active_cm
  ) |>
    as_numeric_draws()
}


eta_from_reference <- function(reference) {

  # Invert one coherent joint clinical reference point, not marginal quantiles.
  baseline_em <- reference[["baseline_em"]]
  baseline_cm <- reference[["baseline_cm"]]

  placebo_em <- reference[["placebo_em"]]
  placebo_cm <- reference[["placebo_cm"]]

  active_em <- placebo_em + reference[["emax_em"]]
  active_cm <- placebo_cm + reference[["emax_cm"]]

  fraction <- c(
    baseline_em / 30,
    (baseline_cm - baseline_em) / (30 - baseline_em),
    -active_em / baseline_em,
    (active_em - active_cm) / (baseline_cm + active_em),
    placebo_em / active_em,
    (placebo_cm - placebo_em) / (active_cm - active_em)
  )

  if (any(!is.finite(fraction) | fraction <= 0 | fraction >= 1)) {
    stop("The clinical reference point is incompatible with the constraints.")
  }

  stats::qlogis(fraction)
}


clinical_quantiles <- function(draws) {

  draws <- as_numeric_draws(draws)

  # Summarize each numeric column, retaining a tidy row for each quantity.
  purrr::imap_dfr(draws, function(values, quantity) {

    quantiles <- stats::quantile(values, c(0.05, 0.50, 0.95), names = FALSE)

    tibble::tibble(
      quantity = quantity,
      q05 = quantiles[1],
      q50 = quantiles[2],
      q95 = quantiles[3]
    )
  })
}


mean_change <- function(draw, dose, population = c("EM", "CM")) {

  # One posterior draw supplies the placebo effect, maximum effect and ED50.
  population <- match.arg(population)
  suffix <- tolower(population)

  placebo <- as.numeric(draw[[paste0("placebo_", suffix)]])
  maximum_effect <- as.numeric(draw[[paste0("emax_", suffix)]])
  ed50 <- as.numeric(draw$ed50)

  placebo + maximum_effect * dose / (dose + ed50)
}

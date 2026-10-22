library(MigraineEmaxSim)

# Reproduce posterior's actual draw-container classes without fitting Stan.
eta_values <- matrix(seq(-2, 2, length.out = 48), nrow = 8, ncol = 6)
eta_columns <- paste0("eta_new[", 1:6, "]")
colnames(eta_values) <- eta_columns

stan_draws <- posterior::as_draws_matrix(eta_values)
stopifnot(inherits(stan_draws, "draws_matrix"))

# The clinical mapping must behave identically for plain and posterior matrices.
expected <- clinical_from_eta(eta_values)
actual <- clinical_from_eta(stan_draws)

plain_columns <- function(table) {
  all(vapply(table, function(column) {
    is.double(column) && is.null(dim(column)) && !inherits(column, "draws")
  }, logical(1)))
}

stopifnot(
  plain_columns(actual),
  isTRUE(all.equal(actual, expected))
)

# A mock exposes the same draws() method used by CmdStanR. It returns a real
# posterior object, so this checks extraction without invoking a Stan sampler.
mock_fit <- list(draws = function(variables) {
  posterior::as_draws_matrix(eta_values[, variables, drop = FALSE])
})

extract <- getFromNamespace("draw_matrix", "MigraineEmaxSim")
reordered <- extract(mock_fit, rev(eta_columns))

stopifnot(
  is.matrix(reordered),
  !inherits(reordered, "draws"),
  identical(colnames(reordered), rev(eta_columns)),
  isTRUE(all.equal(unname(reordered), unname(eta_values[, rev(eta_columns)])))
)

prediction <- predictive_prior_from_meta(mock_fit)

stopifnot(
  plain_columns(prediction$clinical_draws),
  isTRUE(all.equal(prediction$mean, colMeans(eta_values))),
  isTRUE(all.equal(prediction$cov, stats::cov(eta_values) + diag(1e-8, 6)))
)

# Also protect single-row and single-column extraction from dimension dropping.
one_fit <- list(draws = function(variables) {
  posterior::as_draws_matrix(eta_values[1, variables, drop = FALSE])
})

one_column <- extract(one_fit, eta_columns[1])

stopifnot(
  identical(dim(one_column), c(1L, 1L)),
  !inherits(one_column, "draws")
)

# Construct the legacy failure case: one draw-matrix object per clinical column.
legacy <- tibble::as_tibble(lapply(expected, function(column) {
  posterior::as_draws_matrix(matrix(column, ncol = 1))
}))

stopifnot(inherits(legacy$baseline_em, "draws_matrix"))

# Mixed prior/posterior column types must bind and pivot to a numeric value column.
distribution_plot <- plot_distributions(list(Prior = expected, Historical = legacy))

stopifnot(
  is.double(distribution_plot$data$value),
  is.null(dim(distribution_plot$data$value)),
  nrow(distribution_plot$data) == 2L * 6L * nrow(expected)
)

invisible(ggplot2::ggplot_build(distribution_plot))

paired_plot <- plot_em_cm(list(Historical = legacy))
stopifnot(plain_columns(dplyr::select(paired_plot$data, EM, CM)))
invisible(ggplot2::ggplot_build(paired_plot))

# Recovered legacy summaries and Emax mean curves must preserve numerical values.
stopifnot(isTRUE(all.equal(clinical_quantiles(legacy), clinical_quantiles(expected))))

legacy <- dplyr::mutate(legacy, ed50 = 20)
expected <- dplyr::mutate(expected, ed50 = 20)

stopifnot(isTRUE(all.equal(
  mean_change(legacy, 100, "EM"),
  mean_change(expected, 100, "EM")
)))

message("Posterior draw-type regression checks passed.")

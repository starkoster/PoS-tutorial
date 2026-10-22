library(MigraineEmaxSim)

# Exercise Stage 3 extraction, CSV export and plotting using real posterior
# containers. The small mock fit requires neither CmdStan nor new sampling.
eta <- matrix(seq(-2, 2, length.out = 144), nrow = 24, ncol = 6)
colnames(eta) <- paste0("eta[", 1:6, "]")

ed50 <- exp(seq(log(8), log(40), length.out = nrow(eta)))
fit_values <- cbind(eta, ed50 = ed50)

fit <- list(draws = function(variables) {
  posterior::as_draws_matrix(fit_values[, rev(variables), drop = FALSE])
})

expected <- clinical_from_eta(eta) |>
  dplyr::mutate(ed50 = ed50)

draws <- emax_draws(fit)

plain_columns <- function(table) {
  all(vapply(table, function(column) {
    identical(class(column), "numeric") && is.null(dim(column))
  }, logical(1)))
}

# The eight clinical quantities AND potency must be numeric vectors in order.
stopifnot(
  ncol(draws) == 9L,
  nrow(draws) == nrow(eta),
  plain_columns(draws),
  identical(names(draws), names(expected)),
  isTRUE(all.equal(draws, expected, tolerance = 1e-12))
)

# Reproduce the user's exact legacy table: eight matrix-valued columns and
# one ordinary numeric ED50 column. Conversion must retain whole posterior rows.
legacy <- tibble::as_tibble(lapply(expected[1:8], function(column) {
  posterior::as_draws_matrix(matrix(column, ncol = 1))
})) |>
  dplyr::mutate(ed50 = ed50)

stopifnot(inherits(legacy$baseline_em, "draws_matrix"))

temporary_dir <- tempfile("emax-export-")
dir.create(temporary_dir)

for (input in list(draws, legacy)) {

  # Include the saved-RDS recovery path before the reported write_csv() call.
  rds_file <- file.path(temporary_dir, "03_emax_draws.rds")
  saveRDS(input, rds_file)
  recovered <- as_numeric_draws(readRDS(rds_file))

  csv_file <- file.path(temporary_dir, "emax_posterior.csv")
  readr::write_csv(recovered, csv_file)
  round_trip <- readr::read_csv(csv_file, col_types = readr::cols(.default = "d"))

  stopifnot(
    plain_columns(recovered),
    identical(names(round_trip), names(expected)),
    isTRUE(all.equal(as.matrix(round_trip), as.matrix(expected), tolerance = 1e-12))
  )

  # Keep the original input, including the legacy matrix columns, in this plot.
  plot <- plot_distributions(list(Prior = clinical_from_eta(eta), HOPE = input))
  stopifnot(is.double(plot$data$value), is.null(dim(plot$data$value)))
  invisible(ggplot2::ggplot_build(plot))

  # Exercise actual saving as well as ggplot construction.
  save_analysis_plot(plot, "03_hope_update", temporary_dir)
  stopifnot(
    file.info(file.path(temporary_dir, "03_hope_update.pdf"))$size > 0,
    file.info(file.path(temporary_dir, "03_hope_update.png"))$size > 0
  )
}

# A genuine multicolumn matrix or list is not a scalar clinical quantity.
# Refuse these inputs rather than silently flattening and changing row pairing.
invalid <- expected
invalid$baseline_em <- cbind(expected$baseline_em, expected$baseline_em)
stopifnot(inherits(try(as_numeric_draws(invalid), silent = TRUE), "try-error"))

invalid$baseline_em <- as.list(expected$baseline_em)
stopifnot(inherits(try(as_numeric_draws(invalid), silent = TRUE), "try-error"))

# Single-row tables must stay single-row tables, including matrix-valued inputs.
single <- tibble::tibble(value = matrix(3, nrow = 1, ncol = 1))
stopifnot(identical(as_numeric_draws(single)$value, 3))

unlink(temporary_dir, recursive = TRUE)
message("Emax export and plotting regression checks passed.")

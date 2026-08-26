get_pkg_root <- function() {
  if (dir.exists("R") && dir.exists("stan")) {
    return(normalizePath("."))
  }
  if (dir.exists("../R") && dir.exists("../stan")) {
    return(normalizePath(".."))
  }
  stop("Run from the package root or from the analysis directory.")
}

pkg_root <- get_pkg_root()
r_files <- list.files(file.path(pkg_root, "R"), pattern = "[.]R$", full.names = TRUE)
for (f in r_files) source(f)

out_dir <- file.path(pkg_root, "analysis", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

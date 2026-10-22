# Run this script once if the dependencies are not already installed.
# Then install MigraineEmaxSim from its source archive and restart R.

install.packages(c(
  "dplyr", "tidyr", "purrr", "tibble", "readr", "ggplot2",
  "DoseFinding", "posterior"
))

install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org")
)

# Run the following only if CmdStan is not already installed:
# cmdstanr::install_cmdstan()

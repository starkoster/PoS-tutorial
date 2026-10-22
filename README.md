# MigraineEmaxSim 0.5.3

A readable four-stage analysis:

1. Calibrate priors from elicited clinical quantiles.
2. Fit the historical meta-analysis in Stan.
3. Update the drug-specific Emax model in Stan, using a separate ED50 prior.
4. Simulate the competing programs with separate Phase 2 and Phase 3 functions.

Baselines are positive. Changes and active-minus-placebo effects are negative:
an observed `delta <= -2` means a reduction of at least two MMDs versus placebo.

## Install

Unzip the download, open R in the extracted `MigraineEmaxSim` directory, and run:

```r
# Install dependencies once. Skip this if they are already available.
source("analysis/00_install_dependencies.R")

# Install this source package. Restart R after replacing an older version.
install.packages(".", repos = NULL, type = "source")

# Install CmdStan only if it is not already installed.
# cmdstanr::install_cmdstan()
```

The separate `MigraineEmaxSim_0.5.3.tar.gz` archive can also be installed with
`install.packages("/path/to/MigraineEmaxSim_0.5.3.tar.gz", repos = NULL,
type = "source")` after installing the dependencies.

The package imports the tidyverse components it uses: dplyr, tidyr, purrr,
tibble, readr and ggplot2. Installing the full tidyverse is also fine.
Stan sampling uses cmdstanr and an installed CmdStan C++ toolchain.

## Run the four stages

From the extracted package directory, after restarting R:

```r
source("analysis/01_calibrate_prior.R")

source("analysis/02_historical_meta.R")

source("analysis/03_hope_update.R")

source("analysis/04_programs.R")
```

Or run `source("analysis/run_all.R")`. Results go into
`results_R_readable/`, keeping the previous `results_R/` folder separate.
Each script reads the saved output of the preceding stage. The scripts are
also installed under `system.file("analysis", package = "MigraineEmaxSim")`.

Start by reading the four scripts. Then read `R/phase2.R`, `R/phase3.R` and
`R/program_simulation.R` for the simulation details. The optional
`analysis/walk_through_one_program.R` shows the decisions line by line for
one posterior truth.

## Where the main choices live

| Step | Edit here | Main functions |
|---|---|---|
| Clinical quantiles and variance allocation | `analysis/01_calibrate_prior.R` | `calibrate_priors()` |
| Historical arms and active discount | `analysis/02_historical_meta.R` | `read_analysis_arms()`, `fit_meta()` |
| ED50 elicitation and HOPE update | `analysis/03_hope_update.R` | `make_ed50_prior()`, `fit_emax()` |
| Doses, sample sizes, Go rule, TPP thresholds | `R/program_simulation.R` or edit `design` in Stage 4 | `default_program_design()` |
| Phase 2 observations, MCP test and dose selection | `R/phase2.R` | `make_mcp_test()`, `simulate_phase2()` |
| Separate EM/CM Phase 3 trials | `R/phase3.R` | `simulate_phase3()`, `phase3_success()` |
| Competing strategies and paired screening | `R/program_simulation.R` | `simulate_programs()` |

## Phase 2 testing

The package uses `DoseFinding::Mods()` and `DoseFinding::MCTtest()` for the
multiple-contrast testing step of MCP-Mod. CRAN's MCPMod page states that
future development is in DoseFinding:
<https://CRAN.R-project.org/package=MCPMod>.

The Emax candidates have `maxEff = -1`, explicitly making decreasing effects
beneficial. The test uses the general-summary interface with arm-mean
covariance `diag(sigma^2 / n_per_arm)` and `df = Inf`, matching the known-SD
normal simulation model. Its one-sided alpha is explicitly 0.05.

The package computes the critical value once, then reuses it and the contrast
matrix in each simulated trial. This replaces the handwritten contrasts and
200,000-null-simulation calibration. There is no per-trial numerical
integration for p-values. Go uses the package's maximum contrast statistic
and its critical value, together with the observed clinical-effect threshold.

The simulation uses the **MCP testing step**, followed by the existing
smallest-qualifying-dose rule. It does not fit a dose-response model to each
future Phase 2 dataset or use the full `DoseFinding::MCPMod()` modeling step.

Official API documentation:

- <https://openpharma.github.io/DoseFinding/reference/Mods.html>
- <https://openpharma.github.io/DoseFinding/reference/MCTtest.html>

## Simulation assumptions

Each simulated program samples one whole row of the updated Emax posterior.
That truth is held fixed across Phase 2 and both population-specific Phase 3
trials. Sampling errors are independent between trials conditional on truth.
Active-dose comparisons within a trial share one simulated placebo mean.
Patient SD is fixed at 5; Phase 2 has fixed EM/CM arm counts.

Program 1 proceeds directly to fixed-dose Phase 3. Program 2 runs selected-dose
Phase 3 only after Phase 2 Go. The secondary screened Program 1 applies the
same gate to the **same potential fixed-design Phase 3 outcomes**. It therefore
isolates screening from changes in Phase 3 doses or sample sizes.

Both success from program start and success conditional on Go are reported.
A conditional improvement does not imply a larger from-start success
probability for the unchanged Phase 3 design. Joint success is recorded
directly; it is not the product of marginal EM/CM probabilities.

The default multidose Phase 3 benchmark remains nominal, without an across-dose
adjustment. Set `design$p3_multiplicity <- "bonferroni"` to use that explicit
alternative. Safety, development duration and financial utility are outside
these efficacy simulations.

## Priors, constraints and input corrections

Clinical quantiles concern underlying mean parameters in a new study, rather
than patient-level outcomes. Calibration allocates latent variance between
mean uncertainty and study heterogeneity using separately chosen fractions.
The full LKJ correlation prior is included in calibration simulations.
The achieved quantile table should be reviewed; numerical calibration is
approximate and does not establish clinical agreement.

The model retains beneficial-only effects and stronger CM effects through its
clinical transformation. The historical prediction includes heterogeneity;
its transfer to HOPE is a single joint Gaussian approximation on the latent
scale. HOPE shares one ED50 between EM and CM. These are scientific assumptions,
not consequences of simplifying the code.

`inst/extdata/supplied_arms.csv` preserves the supplied data. By default,
`read_analysis_arms()` applies these documented preparations:

- HOPE CM fractions use 62/94, 32/46 and 68/97 for placebo, 100 and 750 mg.
- HOPE placebo baseline SD is corrected from supplied 4.3 to published 4.6.
- Repeated LY pooled-arm baselines are excluded from the subgroup baseline
  likelihood; all change summaries remain included.

The original values remain in `input_cm` and `input_baseline_sd` in the returned
table. `correct_hope_inputs = FALSE` retains the supplied HOPE values for a
deliberate legacy comparison. See `inst/extdata/README.md` for sources.

The HOPE SD correction and package-based MCP calibration mean that the previous
reported PoS numbers are not results for this revision. Run Stages 3 and 4 to
obtain revised results. The historical likelihood itself is unchanged.

## Diagnostics and checks

The custom `check_stan_fit()` gate and its CSV exports have been removed from
the core workflow as requested. CmdStanR's sampling warnings remain enabled.
The optional `analysis/optional_diagnostics.R` gives a short manual review.
Fitting defaults now use 4,000 warmup iterations, 4,000 retained iterations per
chain and `adapt_delta = 0.995`, following the two-divergence run.

Core statistical regression checks are in `tests/core.R`. After installing:

```r
source("tests/core.R")
source("tests/plots.R")
source("tests/draw_types.R")
source("tests/emax_export.R")
```

They cover clinical constraints, data preparation, MCP direction and null
behavior, shared-placebo dependence, and paired Go/Phase 3 identities.
These checks do not require fitting Stan.
The separate plotting checks cover multiple row counts and stages, retained
EM/CM pairing, and construction of the affected figures.
The draw-type checks use actual `posterior::as_draws_matrix()` objects and a
mock extraction interface. They cover conversion to numeric columns, variable
ordering and dimensions, and combining ordinary prior columns with legacy
matrix-valued posterior columns. They do not invoke a Stan sampler.

## Recover the Stage 3 export and plots without refitting

Version 0.5.3 explicitly converts all nine Emax draw columns to ordinary numeric
vectors before saving. Distribution plots also normalize old matrix-valued
clinical columns before combining them with prior draws. The model, data and
posterior values are unchanged by this conversion.

Install this version and **restart R**, without restoring the old patch functions.
The scripts check the version actually loaded in the session, so an old loaded
namespace cannot silently bypass the update. Use the new scripts from this ZIP.

If Stage 3 already fitted successfully and failed at CSV export or plotting:

```r
# Run from the project folder containing your existing results_R_readable/.
# Supply the path to the NEW recovery script from this download.
source("/path/to/MigraineEmaxSim/analysis/recover_stage3.R")
```

Alternatively edit `output_dir` near the top of the recovery script to the
absolute path of your existing results folder. The script reads
`03_emax_draws.rds`, converts its legacy matrix columns and rewrites the RDS,
`emax_posterior.csv`, `emax_summary.csv` and all four Stage 3 figures (PDF/PNG).
If the draw table is absent, it extracts draws from `stan_hope/fit.rds` instead.
**It does not call either Stan fitting function.** Both earlier stages remain
usable; there is no need to repeat prior calibration or the meta-analysis.

For legacy runs, recovery uses the original default ED50 prior for its prior
curve. If you changed that elicitation, enter those same quantiles in the
fallback block of `recover_stage3.R`. New Stage 3 runs save the chosen prior as
`03_ed50_prior.rds`, which recovery reads automatically.

You can also repair an in-memory clinical table directly after installation:

```r
draws <- MigraineEmaxSim::as_numeric_draws(draws)
readr::write_csv(draws, file.path(output_dir, "emax_posterior.csv"))
```

`as_numeric_draws()` preserves values, column names and posterior row pairing.
It accepts numeric vectors and numeric one-column matrices, and rejects list
columns and matrices containing multiple quantities. It does not silently
flatten an invalid clinical table. Stage 4 also normalizes saved inputs.

The new `tests/emax_export.R` reproduces all eight legacy clinical matrix
columns plus ED50, then exercises extraction, an RDS/CSV round trip, the mixed
prior/posterior plot and saving PDF/PNG. It uses real posterior containers and
a mock fit, with no Stan sampling.

R and CmdStan are unavailable in the authoring environment. This source
revision has not been runtime-tested in R or compiled with Stan here.
`R CMD check`, the supplied R tests and a full Stan run remain to be performed.
No fresh posterior or program results are bundled as if those runs occurred.

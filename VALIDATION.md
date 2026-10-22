# Validation status of version 0.5.3

This patch addresses the reported Stage 3 write_csv() and bind_rows() errors.
All nine Emax columns are explicitly normalized to ordinary numeric vectors.
The shared conversion also handles legacy clinical tables before plotting;
whole posterior rows and values are retained. Core scripts check the loaded
namespace version before any fitting. A separate recovery script regenerates
Stage 3 outputs using saved draws or the saved fit, without Stan sampling.

Completed in the authoring environment:

- Reviewed numeric conversion at Stan extraction, clinical mapping, Emax
  extraction, CSV export and plotting boundaries.
- Compared with 0.5.2: Stan files, raw inputs, historical preprocessing,
  calibration, fitting likelihoods and trial-simulation rules are unchanged.
- Checked all 36 exported functions have definitions and help entries.
- Checked matching delimiters and closed strings in R and Stan sources.
  This is a structural scan, not an R parser or compiler check.
- Checked the driver still calls exactly four core stages. Recovery is optional,
  has no fitting calls and repeats the Stage 3 exports and figure calls.
- Checked duplicated analysis scripts match the installed inst/analysis copies.
- Independently checked the clinical transformation and inverse, effect ordering
  and baseline bounds numerically in Python on 10,000 latent draws.
- Checked the raw supplied arm CSV remains unchanged and preprocessing retains
  ten historical change rows and six historical baseline rows.
- Checked archive integrity and expected package, recovery and test contents.

R regression scripts supplied (NOT executed here):

- tests/core.R: clinical constraints, data preparation, MCP testing direction,
  shared-placebo dependence and paired Go/Phase 3 identities.
- tests/plots.R: multiple row counts, matched EM/CM values and related figures.
- tests/draw_types.R: actual posterior containers, matrix extraction and
  mixed-type historical plotting inputs.
- tests/emax_export.R: actual posterior containers plus a mock Emax fit;
  all nine exported columns; legacy RDS repair; CSV write/read round trips;
  the reported mixed Prior/HOPE plot; PDF/PNG saving; invalid-shape rejection.

Not completed because R and CmdStan are unavailable:

- R execution of the supplied functions, analysis scripts and regression tests.
- R CMD build or R CMD check.
- Stan compilation, posterior sampling or new program PoS calculations.

The included tar.gz is a source-installation archive assembled from the
package files, not an R CMD build product. The ZIP also contains the source
analysis directory. Identical scripts are installed under inst/analysis.

To recover the failed Stage 3 output step, install 0.5.3, restart R, and use
analysis/recover_stage3.R with your existing results_R_readable directory.
Use your original ED50 quantiles in its fallback block if they were customized.
For complete new analyses, run the four stages. Native Stan warnings and
optional manual diagnostics remain available.

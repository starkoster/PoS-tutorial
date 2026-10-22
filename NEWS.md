# MigraineEmaxSim 0.5.3

- Added as_numeric_draws(), which explicitly rebuilds a scalar numeric draw
  table and accepts legacy one-column posterior matrices without changing values.
- Guaranteed ordinary numeric vectors in every emax_draws() output column;
  normalized again at Stage 3 export and before combining distribution plots.
- Added analysis/recover_stage3.R to regenerate CSVs and figures from the saved
  clinical draws or saved Emax fit, without refitting either Stan model.
- Added loaded-namespace version checks before core scripts do any work and
  qualified critical Stage 3 calls to avoid older global patch definitions.
- Stage 4 accepts legacy saved draw tables; Stage 3 saves the ED50 elicitation
  for accurate prior-curve recovery in subsequent sessions.
- Added regression coverage for the exact Stage 3 export and plotting calls,
  including all eight matrix-valued clinical quantities and numeric ED50.
- No changes to data, Stan models, likelihoods or simulation rules. R runtime
  validation remains unavailable in the authoring environment.

# MigraineEmaxSim 0.5.2

- Fixed posterior class leakage: Stan extraction now explicitly reconstructs
  a base numeric matrix before selecting variables and transforming them.
- Clinical transformation outputs are ordinary numeric columns even when
  input is a posterior draws_matrix. Labels and posterior row order are retained.
- Distribution and calibration plots normalize legacy matrix-valued columns
  before combining stages; EM/CM plots, quantiles and mean curves also accept them.
- Added regression cases using actual posterior::as_draws_matrix objects and
  mock draw extraction, including column reordering and a mixed-type plot input.
- Included a standalone patch for existing sessions. Existing Stan fits and
  saved clinical values can be reused; the numerical model is unchanged.
- R runtime validation remains unavailable in the authoring environment.

# MigraineEmaxSim 0.5.1

- Fixed plot_em_cm(): a quantity column created inside tibble() masked the
  scalar callback variable after recycling to multiple rows. EM/CM vectors
  are now extracted before constructing the plotting table.
- Renamed local plotting callback arguments so they differ from output columns.
- Added focused plotting regression checks for one or multiple stages, one row
  or thousands of rows, correct EM/CM pairing, and plot construction.
- Included a standalone plotting patch so existing calibration objects can be
  used immediately, without reinstalling or repeating the calibration.
- No statistical model, fitting inputs, calibration calculation or simulation
  rule changed in this patch. R runtime validation remains unavailable here.

# MigraineEmaxSim 0.5.0

- Reworked around four short analysis scripts and separately named R modules.
- Added more whitespace, explicit intermediate quantities and comments inside
  every function; tidyverse tools handle tables, joins and summaries.
- Retained Stan for historical meta-analysis and the HOPE Emax update.
- Replaced custom MCP contrast/null-simulation code with DoseFinding's Mods
  and MCTtest, using decreasing effects and the known-variance summary interface.
- Kept Phase 2 and Phase 3 simulation separate, with a plainly written program
  loop and paired fixed-design screening for Program 1.
- Removed the custom diagnostic stopping gate and large diagnostic outputs
  from the core workflow. Native CmdStanR warnings remain; manual checks are
  in a separate optional script.
- Changed fitting defaults to adapt_delta 0.995 and 4,000 warmup iterations.
- Corrected HOPE placebo baseline SD to 4.6 during data preparation, retaining
  the supplied value and arm-composition inputs in separate columns.
- Added a one-program walkthrough and focused statistical regression checks.
- This is a source revision, with no newly fitted results. R/Stan runtime
  validation remains outstanding in the authoring environment.

## Interface changes from 0.4.0

- Use `analysis/run_all.R` instead of the large `run_analysis()` wrapper.
- `check_stan_fit()` is no longer called or exported.
- `make_mcp_test()` no longer takes `n_null` or stores an empirical null sample.
- Phase 2 returns the decision and maximum contrast statistic; it does not
  calculate an adjusted p-value for every simulated program.
- `simulate_prior()` focuses on the full new-study prior; optional prior-layer
  and prior-width diagnostic machinery is omitted from this streamlined package.
- Program summaries retain probabilities, denominators and MCSEs; the larger
  simulation-interval and diagnostic reporting machinery is omitted.
- `plot_selection()` now takes the result of `summarize_programs()`.
- Existing numerical results, C++/Python samplers and manuscript-generation
  tools are not included in this focused R/Stan source package.

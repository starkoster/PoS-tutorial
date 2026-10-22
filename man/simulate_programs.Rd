\name{simulate_programs}
\alias{simulate_programs}
\alias{default_program_design}
\title{Direct and staged development programs}
\description{Simulate the competing strategies and a paired fixed-design
screening analysis using a plainly written Phase 2/Phase 3 loop.}
\usage{
default_program_design()
simulate_programs(draws, mcp, n_sim = 50000L, seed = 4601,
                   design = default_program_design())
}
\arguments{
\item{draws}{Joint Emax posterior rows.}
\item{mcp}{Phase 2 test specification matching the design.}
\item{n_sim}{Number of programs.}
\item{seed}{Program-simulation seed.}
\item{design}{Named list returned by default_program_design(), optionally edited.}
}
\value{A list with event indicators, Phase 2 decisions and arm means, direct
and selected Phase 3 comparisons, generating truths and design settings.}
\details{One truth is held fixed within each program. Program 1 skips Phase 2
in its primary strategy. Program 2 runs selected-dose Phase 3 only after Go.
Screened Program 1 applies the gate to the same fixed Phase 3 outcomes,
including counterfactual outcomes for No Go programs. Stopped programs are
unsuccessful from program start. Conditional and from-start PoS differ.}

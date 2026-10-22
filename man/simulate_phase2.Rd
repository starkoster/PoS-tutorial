\name{simulate_phase2}
\alias{simulate_phase2}
\alias{make_mcp_test}
\title{Phase 2 simulation and package-based MCP testing}
\description{Generate mixed-population arm means, test for a decreasing
dose-response relationship and apply the observed-effect Go rule.}
\usage{
make_mcp_test(doses, ed50_candidates, n_per_arm = 140L,
              sigma = 5, alpha = 0.05, seed = 4501)
simulate_phase2(draw, mcp, prop_em = 0.5, go_delta = -2)
}
\arguments{
\item{doses}{Distinct nonnegative doses including one placebo and at least two active doses.}
\item{ed50_candidates}{Positive candidate ED50 values, fixed before future data.}
\item{n_per_arm}{Total participants per arm.}
\item{sigma}{Known common patient SD.}
\item{alpha}{One-sided MCP significance level.}
\item{seed}{Seed used for the package's critical-value integration.}
\item{draw}{One joint clinical posterior row including ED50.}
\item{mcp}{Design-specific test specification returned by make_mcp_test().}
\item{prop_em}{EM fraction, converted to a fixed rounded count per arm.}
\item{go_delta}{Negative observed active-minus-placebo Go threshold.}
}
\value{A reusable test specification, or a list containing Go, selected dose,
significance, maximum contrast statistic and simulated arm means.}
\details{DoseFinding::Mods() and MCTtest() supply the contrasts, critical value
and test. A negative maximum effect specifies decreasing benefit. The general
interface uses known arm-mean covariance and infinite degrees of freedom.
Go requires MCP significance and a qualifying observed difference. The smallest
qualifying dose is selected. The modeling step of MCP-Mod is not performed.}

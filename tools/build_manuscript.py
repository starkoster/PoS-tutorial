"""Build editable numerical tables and replacement Section 5 from actual results."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2];R=ROOT/'results';M=ROOT/'manuscript';M.mkdir(exist_ok=True)
LABEL={'baseline_em':'Baseline EM','baseline_cm':'Baseline CM','placebo_em':'Placebo change EM',
 'placebo_cm':'Placebo change CM','emax_em':'Maximum drug effect EM','emax_cm':'Maximum drug effect CM',
 'active_em':'Maximum active change EM','active_cm':'Maximum active change CM','ed50':'$ED_{50}$ (mg)'}
def read(name):return pd.read_csv(R/f'tables/{name}.csv')
def tab(name,headers,rows,align=None):
    align=align or 'l'+'r'*(len(headers)-1)
    text='\\begin{tabular}{@{}'+align+'@{}}\n\\toprule\n'+' & '.join(headers)+'\\\\\n\\midrule\n'
    text+='\n'.join(' & '.join(map(str,row))+'\\\\' for row in rows)+'\n\\bottomrule\n\\end{tabular}\n'
    (M/(name+'.tex')).write_text(text)
def interval(row):return f'{row.q50:.2f} ({row.q05:.2f}, {row.q95:.2f})'

c=read('calibration_quantiles');p=read('calibrated_prior')
tab('table_calibration',['Quantity','Target 5/50/95\\%','Achieved 5/50/95\\%'],
 [[LABEL[r.quantity],f'{r.q05_target:.2f}, {r.q50_target:.2f}, {r.q95_target:.2f}',f'{r.q05_achieved:.2f}, {r.q50_achieved:.2f}, {r.q95_achieved:.2f}'] for r in c.itertuples()],'lll')
tab('table_latent',['Coordinate','$m_l$','$s_l$','$q_l$','$h_l$'],
 [[f'$\\eta_{{{int(r.dimension)}}}$',f'{r.m:.4f}',f'{r.s:.4f}',f'{r.q:.4f}',f'{r.heterogeneity_fraction:.4f}'] for r in p.itertuples()])
w=read('prior_width_sensitivity')
tab('table_width',['Width factor','$\\Pr(B_{\\mathrm E}>15)$','$\\Pr(D_{\\mathrm E}\\leq-2)$','$\\Pr(D_{\\mathrm C}\\leq-2)$'],
 [[f'{r.width_multiplier:.2f}',f'{r.p_em_baseline_gt15*100:.2f}\\%',f'{r.p_em_emax_le_minus2*100:.2f}\\%',f'{r.p_cm_emax_le_minus2*100:.2f}\\%'] for r in w.itertuples()])
pair=read('prior_between_study_differences').iloc[:6]
tab('table_heterogeneity',['Quantity','5th percentile of difference','95th percentile of difference'],
 [[LABEL[r.quantity],f'{r.q05:.2f}',f'{r.q95:.2f}'] for r in pair.itertuples()])
d=read('analysis_arms')
tab('table_arms',['Study','Arm/dose (mg)','$n$','CM\\%','Baseline (SD)','Change (SE)'],
 [[r.title,('Pbo' if r.arm=='placebo' else str(r.dose)),str(r.n),f'{r.CM*100:.1f}',
   f'{r.base_mmd:.1f} ({r.base_mmd_sd:.2f})'+('$^*$' if r.title=='LY3451838' else ''),f'{r.chg:.1f} ({r.chg_se:.2f})'] for r in d.itertuples()],'llrrll')
h=read('historical_predictive_summary').set_index('quantity');g=read('hope_prior_summary').set_index('quantity');e=read('emax_posterior_summary').set_index('quantity')
tab('table_update',['Quantity','Historical prediction','Gaussian HOPE prior','HOPE posterior'],
 [[LABEL[q],interval(h.loc[q]),interval(g.loc[q]),interval(e.loc[q])] for q in h.index],'llll')
probs=read('program_probabilities');derisk=read('fixed_design_derisking')
rows=[]
for threshold in [0,-2,-3,-4]:
    for strat,label in [('P1_direct','P1 direct'),('P1_screened','P1 screened'),('P2_selected','P2 selected')]:
        sub=probs[(probs.threshold==threshold)&(probs.strategy==strat)&(probs.estimand=='from_start')].set_index('population')
        cond=probs[(probs.threshold==threshold)&(probs.strategy==strat)&(probs.estimand=='given_go')].set_index('population')
        rows.append(['Statistical' if threshold==0 else f'$\\leq {threshold}$',label,
          f'{sub.loc["EM","probability"]*100:.2f}',f'{sub.loc["CM","probability"]*100:.2f}',
          f'{sub.loc["both","probability"]*100:.2f}',f'{cond.loc["both","probability"]*100:.2f}'])
tab('table_programs',['Threshold','Strategy','EM','CM','Both','Both given Go'],rows,'llrrrr')
rows=[]
for threshold in [0,-2,-3,-4]:
    r=derisk[(derisk.threshold==threshold)&(derisk.population=='both')].iloc[0]
    rows.append(['Statistical' if threshold==0 else f'$\\leq {threshold}$',*[f'{r[col]*100:.2f}' for col in ['direct','given_go','given_no_go','gated_from_start','conditional_gain']]])
tab('table_derisking',['Threshold','$\\Pr(S)$','$\\Pr(S\\mid G)$','$\\Pr(S\\mid G^c)$','$\\Pr(G\\cap S)$','Gain (pp)'],rows,'lrrrrr')
go=read('phase2_go').iloc[0];sel=read('phase2_selection');sample=read('expected_participants')
tab('table_selection',['Selected dose (mg)','Count','Percentage among Go'],[[str(r.dose),str(r.count),f'{100*r.fraction_of_go:.2f}\\%'] for r in sel.itertuples()])
md=pd.read_csv(R/'diagnostics/meta_mcmc.csv');hd=pd.read_csv(R/'diagnostics/hope_mcmc.csv')
tab('table_diagnostics',['Model','Chains','Warmup/chain','Retained/chain','Max $\\widehat R$','Min bulk ESS'],
 [['Historical','4','15,000','30,000',f'{md.rhat.max():.4f}',f'{md.ess_bulk.min():.0f}'],
 ['HOPE','4','8,000','20,000',f'{hd.rhat.max():.4f}',f'{hd.ess_bulk.min():.0f}']],'lrrrrr')

text=r'''\section{Numerical specification and results}

The following results use newly calibrated priors and newly fitted models;
they do not reuse the previous package's posterior or PoS results. All effect
summaries are negative changes or active-minus-placebo differences. Unless
specified otherwise, parameter intervals are central 90\% intervals.

\subsection{Proposed clinical quantiles and variance calibration}

The six proposed clinical targets and their independently validated achieved
quantiles are shown in Table~\ref{tab:calibration}. They concern new-study
mean parameters, before incorporating the supplied historical data. They are
analyst proposals retained for this tutorial, not documented stakeholder
elicitation or evidence that experts have accepted the joint prior.
Structured elicitation requires recording the target population, evidence
available to experts, individual judgments and reconciliation of the resulting
distribution~\cite{elicitation}.

\begin{table}[H]\centering\small
\caption{Clinical target quantiles and independently validated calibration, in
MMDs. For negative effects, the fifth percentile is the most beneficial endpoint.}
\label{tab:calibration}\input{table_calibration.tex}\end{table}

The clinical mapping used $M=30$ MMDs. The correlation prior was
$\Omega\sim\mathrm{LKJ}(2)$, which has the identity as its modal correlation
matrix while permitting nonzero correlations~\cite{stanlkj}. The heterogeneity
variance fractions were fixed at $h_1=h_2=4/13$ and $h_3=\cdots=h_6=0.36$,
retaining the allocation in the earlier analysis. This allocation is a separate
assumption and is not identified by the six marginal clinical intervals.

Calibration used 20,000 fixed prior simulations and three starting points.
Every simulation included a full LKJ(2) matrix. Independent validation used
100,000 new draws, including new correlation matrices. The numerical reference
used bounded nonlinear least squares; the R implementation uses the same
quantile loss and Nelder--Mead optimization. Their finite-simulation candidates
need not coincide. The exact latent settings used for the reported results
are provided at full precision in the accompanying CSV and summarized below.

\begin{table}[H]\centering\small
\caption{Calibrated latent prior specification. The $q_l$ column is the
half-normal scale of $\tau_l$, not a fixed between-study SD.}
\input{table_latent.tex}\end{table}

The largest relative quantile mismatch was 10.08\%, at the median EM maximum
drug effect: the target was $-0.50$ MMD and the achieved median was $-0.55$.
The achieved EM fifth-percentile effects were approximately $-2.32$ for
placebo and $-2.31$ for maximum drug effect, against targets of $-2.50$.
Calibration is approximate. Multiple starts did not establish a global optimum
or prove that the remaining discrepancy cannot be reduced.

\fig{01_prior_calibration.pdf}{1}{Calibrated new-study clinical priors (blue)
and proposed fifth, median and ninety-fifth percentile targets (orange).
These distributions describe study mean parameters, without patient sampling
error. Displayed density curves are smoothed; numerical quantiles use the draws.}

\fig{03_prior_layers.pdf}{1}{Prior variance decomposition. Grey varies only
the reference-class mean; light blue fixes that mean at its center and varies
new-study heterogeneity; dark blue includes both. Nonlinearity means that these
clinical distributions cannot be combined by simply adding their quantiles.}

\fig{02_prior_em_cm.pdf}{1}{Dependence implied by the calibrated clinical
mapping. Top: independent latent study deviations. Bottom: integration over
the full LKJ(2) prior. All remaining prior settings are the same. Columns show
baseline, placebo change and total active change; darker hexagons contain more
draws. The diagonal marks equality. Similar integrated patterns do not establish
posterior or PoS insensitivity to particular nonzero correlation matrices.}

The strong EM/CM dependence persists with identity latent correlation because
the mapping shares baselines and effects. The prior's modal matrix and this
visual similarity are the proposed rationale for regularizing additional
latent correlations; neither is a mathematical demonstration of decorrelation.

To interpret the heterogeneity scales independently of uncertainty in the
mean, two new studies were simulated under each common draw of
$(\boldsymbol\mu,\boldsymbol\tau,\Omega)$. Table~\ref{tab:heterogeneity} reports
the distribution of their clinical differences. Shared mean uncertainty is
not independently redrawn between the pair. These are implications to assess
with clinicians, rather than additional elicited targets.

\begin{table}[H]\centering\small
\caption{Central 90\% intervals for differences between two comparable
studies, conditional on sharing hyperparameters and then integrated over their
prior. Values are differences in underlying mean MMD parameters.}
\label{tab:heterogeneity}\input{table_heterogeneity.tex}\end{table}

Multiplying both $s_l$ and $q_l$ by a common factor changes initial-prior
efficacy probabilities (Table~\ref{tab:width}). These are probabilities about
asymptotic true effects, not trial PoS. Under the calibrated settings,
$\Pr(B_{\mathrm E}>15)=4.34\%$. If the intended EM eligibility definition
excludes such baseline means, the support or targets require revision before
using the prior as a clinically accepted specification.

\begin{table}[H]\centering\small
\caption{Initial-prior width sensitivity with latent centers held fixed.
No historical or HOPE refitting is performed in this check.}
\label{tab:width}\input{table_width.tex}\end{table}

\subsection{Historical and same-drug input summaries}

The supplied data include 462 participants from selected CONQUER rows
(NCT03559257), 591 from DELIVER rows (NCT04418765), 38 from LY3451838 rows
(NCT04498910), and 237 from HOPE (NCT05133323). These counts describe the
included rows, not necessarily every trial arm or the complete analysis
populations. Historical active doses were 120, 300 and 1500 mg, respectively;
their numerical values were not used as a pooled dose-response predictor.
HOPE studied Lu AG09222, the investigational drug being updated.

\begin{table}[H]\centering\small
\caption{Numerical arm inputs. Baseline entries are mean (SD); change entries
are mean (SE). $^*$Repeated LY baseline summaries are shown for auditability
but omitted from the baseline likelihood. All LY change rows remain included.}
\input{table_arms.tex}\end{table}

Two explicit input decisions distinguish this fit from the earlier one.
HOPE uses the published arm-specific CM counts: 62/94 on placebo, 32/46 at
100 mg and 68/97 at 750 mg~\cite{hope}. This replaces the approximate common
70\% composition. The supplied LY baseline mean and SD are identical in its
EM and CM rows within each treatment arm. Without verified subgroup baseline
provenance, these baseline terms are excluded; their repeated values are
preserved in the original CSV. No unreported subgroup baselines were imputed.

The historical change likelihood used all ten supplied rows. Active change
terms had power 0.50, placebo terms power one. The six retained CONQUER and
DELIVER baseline rows had power one. This discounts active cross-drug
information while retaining disease-burden information. The selected active
responses are only proxies for maximum efficacy. In addition to the normal
arm-summary approximation, the unreported within-arm baseline/change
covariance remains a limitation.

\fig{12_historical_inputs.pdf}{1}{Historical inputs with normal-approximation
95\% intervals. Mixed summaries are labeled as mixed. Baseline observations
excluded from the likelihood are omitted from the baseline panel. Historical
active changes remain negative, matching the model and subsequent figures.}

\subsection{Historical update and predictive transfer}

Historical fitting moves the initial prior toward larger placebo and active
changes, while retaining uncertainty about a new study. The new-study
predictive medians were $-0.92$ and $-2.38$ MMDs for placebo change, and
$-1.60$ and $-3.05$ MMDs for maximum drug effects in EM and CM. Predictive
median baselines were 9.47 and 18.40 MMDs. These are predictions for a
new study; they are not the reference-class means or fixed historical study
estimates.

\fig{04_meta_marginals.pdf}{1}{Initial new-study prior (grey) and historical
posterior new-study prediction (blue) in each of the six elicited clinical
dimensions. Prediction retains between-study heterogeneity.}

\fig{05_meta_studies.pdf}{1}{Historical study posterior medians and central
90\% intervals alongside initial and posterior new-study predictive intervals.
The LY baseline distributions are model-derived because its ambiguous baseline
rows were excluded. DELIVER subgroup effects are inferred from mixed data and
the hierarchical model; they are not separately observed EM and CM estimates.}

\fig{06_meta_em_cm.pdf}{1}{EM/CM cross-checks for baseline, placebo change
and total active change, before and after historical evidence. Orange points
are model-derived historical study posterior medians. They are not observed
subgroup pairs for trials reporting mixed summaries. The lower row is a
new-study prediction rather than the distribution of the three historical means.}

\fig{13_historical_ppc.pdf}{1}{Historical arm checks. Orange points are
observed summaries; dark blue gives central 90\% intervals for their underlying
model means; light blue gives central 90\% observational-replication intervals
using the original reported SEs. Active likelihood discounting affects the fit,
not the observational SE used for these replications.}

The empirical predictive latent draws were summarized by their mean and
covariance and transferred through a single Gaussian approximation, with
covariance multiplier one and diagonal ridge $10^{-8}$. Figure~\ref{fig:transfer}
shows the approximation on the clinical scale. It changes some shapes, most
visibly the EM baseline distribution. The achieved quantile differences are
provided in \texttt{gaussian\_transfer\_check.csv}; they should not be hidden
by describing the Gaussian transfer as exact Bayesian propagation.

\begin{figure}[H]\centering
\includegraphics[width=\linewidth]{07_transfer_check.pdf}
\caption{Empirical historical new-study prediction (blue) and the Gaussian
latent approximation used as the HOPE prior (orange). Both are transformed
through the same constrained clinical mapping.}\label{fig:transfer}\end{figure}

\begin{table}[H]\centering\scriptsize
\caption{Clinical parameter medians and central 90\% intervals, in MMDs.
Historical prediction and transferred HOPE prior are distinguished explicitly.}
\resizebox{\linewidth}{!}{\input{table_update.tex}}\end{table}

\subsection{HOPE Emax update and dose-response results}

The separately elicited potency prior was
\[
 \log ED_{50}\sim N(\log20,0.69994^2).
\]
Its median was 20 mg and central 90\% interval 6.32--63.25 mg. These are
the actual implemented quantiles; 10 and 100 mg are not substituted as if
they were the fifth and ninety-fifth percentiles of this same prior.

After HOPE, median maximum drug effects were $-1.69$ MMDs in EM
(90\% interval $-3.07$ to $-0.67$) and $-3.23$ MMDs in CM
($-4.77$ to $-2.00$). Median placebo changes were $-1.28$ and $-4.62$ MMDs.
The posterior $ED_{50}$ median was 19.74 mg, with 90\% interval
6.60--59.81 mg and mean 24.85 mg. Potency remains relatively weakly informed
by two active doses, as seen in the prior/posterior comparison.

\fig{08_hope_update.pdf}{1}{Transferred clinical prior (grey) and HOPE
posterior (blue). Population-specific changes reflect the joint model and
historical borrowing as well as the mixed-population HOPE observations.}

\fig{09_emax_adjusted.pdf}{1}{HOPE posterior active-minus-placebo dose
responses. Lines are pointwise medians and bands pointwise central 90\%
intervals, rather than simultaneous credible bands. The orange line marks
$-2$ MMDs. Negative effects indicate improvement.}

\fig{10_emax_change.pdf}{1}{Posterior total mean change from baseline,
including placebo response. These curves differ from placebo-adjusted effects
and approach the maximum active changes at high dose.}

\fig{11_ed50_update.pdf}{0.65}{Separately specified lognormal $ED_{50}$ prior
and HOPE posterior. The similarity of these marginal densities indicates
limited updating of potency under this specification. The joint posterior,
including dependence between potency and maximum effects, is retained for
simulation.}

\subsection{Development-program settings and simulation outcomes}

\begin{table}[H]\centering\small
\caption{Prospective designs; Phase~3 arm sizes apply separately to EM and CM.}
\begin{tabular}{@{}p{3.5cm}p{4.5cm}p{5.5cm}@{}}\toprule
Strategy/stage & Doses (mg) & Allocation\\\midrule
Program 1, direct Phase 3 & 0, 20, 200 & 233 per arm in each population\\
Common Phase 2 & 0, 20, 100, 300 & 140 per arm; fixed 70 EM and 70 CM\\
Program 2, selected Phase 3 & 0, selected dose & 174 per arm in each population\\
Program 1, secondary screening & Same Phase 2, then fixed 0, 20, 200 & Phase 3 unchanged at 233 per arm\\
\bottomrule\end{tabular}\end{table}

The 200 mg higher Program~1 dose resolves the earlier discrepancy between
the executed configuration and an older 100 mg slide label. Patient SD was
fixed at 5 MMDs. The Phase~2 Go threshold was an observed difference
$\leq-2$ MMDs together with a one-sided MCP rejection at 0.05. The three
candidate $ED_{50}$ values were fixed at posterior quantiles 8.46, 19.74 and
47.08 mg. With 200,000 null simulations, the maximum-statistic critical value
was 1.7849. An independent 100,000-null check gave rejection rate 0.05127
(Monte Carlo SE 0.00070), compatible with the nominal level within simulation
variation.

Phase~3 tests were beneficial two-sided normal tests at nominal 0.05 without
multi-dose adjustment. The package offers Bonferroni adjustment as an explicit
alternative. TPP thresholds were observed differences of $\leq-2$, $\leq-3$
or $\leq-4$ MMDs. At the specified sample sizes and SD, a $-2$ MMD observed
difference already exceeds the nominal significance boundary, about $-0.91$
for 233 per arm and $-1.05$ for 174 per arm. Consequently, the nonzero TPP
scenarios are largely determined by their observed-effect thresholds.

Across 50,000 programs, Phase~2 Go occurred in 35,821 (71.64\%). The Go
probability's conditional-on-posterior Monte Carlo SE was 0.20 percentage
points. Dose selection among advancing programs is reported below.

\begin{table}[H]\centering\small
\caption{Smallest qualifying Phase~2 dose among programs receiving Go.}
\input{table_selection.tex}\end{table}

\fig{16_phase2_selection.pdf}{0.65}{Selected doses among advancing programs.
Selection is based on observed Phase~2 effects following the MCP Go requirement,
without refitting a future Emax model or changing the generating truth.}

\begin{table}[H]\centering\small
\caption{Efficacy success percentages. EM, CM and Both use all 50,000 programs
as denominator. The final column conditions on the 35,821 common Go programs.
The conditional column for direct Program~1 is a screening diagnostic.}
\input{table_programs.tex}\end{table}

\fig{14_program_comparison.pdf}{1}{Success from program start for direct
Program~1, screened Program~1 and selected-dose Program~2. Error bars are
95\% Wilson intervals for simulation noise conditional on the retained
posterior. ``Stat.'' uses statistical success without an additional effect
threshold. Complete strategies differ in dose, sample size and gating.}

At the $-2$ MMD target, direct Program~1 had joint success 27.82\%
(conditional-on-posterior simulation interval 27.43--28.21\%). Program~2
had joint success from program start 16.57\% (16.25--16.90\%) and conditional
joint success given Go 23.13\%. This strategy comparison does not isolate
de-risking: Program~2 selects a dose, often below 200 mg, uses fewer Phase~3
participants per arm and requires a Phase~2 Go event. EM is the limiting
population under these assumptions; the mixed-population Go rule does not
require its EM effect to meet the Phase~3 target.

\begin{table}[H]\centering\small
\caption{Paired fixed-design screening analysis for success in both EM and CM.
Probabilities and gains are percentages and percentage points, respectively.}
\input{table_derisking.tex}\end{table}

With the Phase~3 design held fixed, the $-2$ MMD joint success probability
increased from 27.82\% overall to 35.37\% conditional on Go, an enrichment
of 7.55 percentage points (paired Monte Carlo SE 0.112 percentage points).
Among No Go programs, potential joint Phase~3 success was 8.74\%.
Success from the start of the screened Program~1 was 25.34\%
(simulation interval 24.96--25.73\%). The identity
$0.71642\times0.353731=0.25342$ makes the denominators explicit. A gate
removes some programs that would have succeeded, so the conditional gain
does not imply higher unconditional gated success than the unchanged direct
design.

\fig{15_fixed_design_derisking.pdf}{1}{Screening isolated using shared
potential Phase~3 outcomes. Grey is direct success; blue conditions the same
outcomes on Phase~2 Go; orange requires both Go and Phase~3 success from
program start. No increase in conditional success is imposed by the code.}

Joint success may occur at different doses in EM and CM. For the $-2$ MMD
scenario, direct Program~1 common-dose success was 27.66\%, compared with
27.82\% when different doses were allowed. The corresponding screened
from-start figures were 25.20\% and 25.34\%. Program~2 uses one selected
dose in both populations, so its two definitions coincide.

Expected participant counts were 1398 for direct Program~1,
1561.56 for screened Program~1 and 1058.63 for Program~2, including the
560-participant Phase~2 trial where applicable. These counts are design
consequences rather than a cost-effectiveness or development-time analysis.

\fig{17_phase3_joint_effects.pdf}{1}{Observed active-minus-placebo effects
in the EM and CM Phase~3 trials. The first two panels show potential fixed-design
Program~1 effects at each dose; the last shows the selected-dose trials actually
run after Go in Program~2. Orange lines mark the $-2$ MMD thresholds.}

\subsection{Computation, reproducibility and remaining scientific checks}

The supplied numerical execution used independent C++/Python code for the
same likelihoods, priors and constrained mapping as the R/Stan package.
The historical sampler combined Gaussian Gibbs updates for the latent mean,
elliptical-slice study updates, Metropolis updates for log heterogeneity SDs
and canonical partial correlations, and noncentered interweaving for the SDs.
Its partial-correlation priors imply the stated LKJ(2) distribution. The HOPE
sampler used joint elliptical-slice updates under its Gaussian latent and
log-potency prior. Warmup adaptation was discarded.

\begin{table}[H]\centering\small
\caption{Reference-execution diagnostics over monitored model parameters.
Rank-normalized split/folded $\widehat R$ and bulk ESS follow the definitions
used for modern MCMC assessment~\cite{rhat}. These are not Stan diagnostics.}
\input{table_diagnostics.tex}\end{table}

Prior-only simulation checked that the historical sampler recovered the
specified mean hyperpriors, half-normal scale moments, latent new-study
variances and LKJ pair moments. The LKJ pair variance was approximately
0.111, against the analytical value $1/9$. Diagnostics, individual chains,
trace plots, original inputs, analysis-ready inputs, calibrated settings and
program records are included. The numerical reference is independently
reproducible from its source.

R and CmdStan were unavailable in the execution environment. The R functions
and Stan models have therefore not been run or validated by \texttt{R CMD check}
here. The package includes executable R regression checks for sign constraints,
the inverse transformation, LKJ moments, null behavior, shared placebo
dependence and paired screening identities. An R/Stan rerun, including its
sampler diagnostics, remains necessary to validate that implementation.
The two optimizers and random-number generators can yield different
finite-simulation calibration candidates; the supplied full-precision prior
CSV identifies the candidate used in every reported result.

The principal scientific issues remain acceptance of the quantile targets
and variance allocation, the EM baseline upper tail, beneficial-only and
CM-dominance restrictions, cross-mechanism active-effect borrowing and the
maximum-dose proxy, normal arm-summary likelihoods, the Gaussian transfer,
shared potency, and fixed future variance. Prior-width plots are calibration
checks, not sensitivity analyses of fitted PoS. Assessing these assumptions'
decision consequences requires targeted refitting and new simulations.
The historical selection and subgroup baseline provenance also require a
documented evidence audit before a clinical or regulatory use.

\begin{thebibliography}{9}
\bibitem{hope} Ashina M et al. A monoclonal antibody to PACAP for migraine
prevention. \emph{New England Journal of Medicine}, 2024.
\href{https://doi.org/10.1056/NEJMoa2314577}{doi:10.1056/NEJMoa2314577}.
The arm-specific CM counts are from Table 1.
\bibitem{stanlkj} Stan Development Team. Correlation matrix distributions.
\url{https://mc-stan.org/docs/functions-reference/correlation_matrix_distributions.html}.
Accessed 15 October 2026.
\bibitem{elicitation} O'Hagan A. Expert knowledge elicitation: subjective but
scientific. Author's article summary and elicitation discussion.
\url{https://tonyohagan.co.uk/academic/abs/ElicSubjSci.html}.
\bibitem{rhat} Vehtari A, Gelman A, Simpson D, Carpenter B, B\"urkner PC.
Rank-normalization, folding, and localization: an improved $\widehat R$ for
assessing convergence of MCMC. \emph{Bayesian Analysis} 16:667--718, 2021.
\href{https://doi.org/10.1214/20-BA1221}{doi:10.1214/20-BA1221}.
\end{thebibliography}
'''
(M/'section5.tex').write_text(text)
print('Tables and Section 5 written from the current numerical outputs.')

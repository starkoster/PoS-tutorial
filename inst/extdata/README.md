# Arm-summary inputs

`supplied_arms.csv` is the unmodified supplied input. `clinical_targets.csv`
contains the proposed clinical quantiles also returned by `default_targets()`.

The preparation function retains original HOPE fields and, by default,
substitutes the published arm-specific CM fractions and placebo baseline SD.
HOPE Table 1 gives placebo baseline 16.9 (SD 4.6), whereas the supplied CSV
uses SD 4.3. The latter is the paper's overall baseline SD. This revision
corrects the arm-specific likelihood input; old fitted results are not reused.

HOPE: Ashina et al., New England Journal of Medicine 2024;391:800–809.
<https://doi.org/10.1056/NEJMoa2314577>. The drug is Lu AG09222, the same
drug whose Emax curve is updated, not a historical cross-compound reference.

The repeated LY3451838 baseline summaries refer to pooled treatment arms,
not the separately reported EM/CM subgroups. They are excluded from the
baseline likelihood. All four LY change rows remain, using the supplied
month 3 subgroup changes rather than the month 1 primary endpoint.

LY3451838 report: Johnson et al., Cephalalgia 2025;45(8).
<https://doi.org/10.1177/03331024251368757>.

Historical reference trials:

- CONQUER, NCT03559257: <https://clinicaltrials.gov/study/NCT03559257>
- DELIVER, NCT04418765: <https://clinicaltrials.gov/study/NCT04418765>
- LY3451838, NCT04498910: <https://clinicaltrials.gov/study/NCT04498910>
- Same-drug update HOPE, NCT05133323: <https://clinicaltrials.gov/study/NCT05133323>

Other supplied summary values remain as given. Mechanism, route, prior-failure
definitions and assessment windows are not completely harmonized. The code
does not imply an independent reconstruction of every estimator from raw data.

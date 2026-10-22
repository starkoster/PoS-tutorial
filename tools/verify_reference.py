"""Checks of supplied numerical outputs. Run after run_reference.py."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT=Path(__file__).resolve().parents[2];R=ROOT/'results'
checks=[]
def check(name,value):
    checks.append({'check':name,'passed':bool(value)})
    if not value:raise AssertionError(name)
def clinical_support(x):
    return np.all((x.baseline_em>0)&(x.baseline_cm>x.baseline_em)&(x.baseline_cm<30)&
     (x.placebo_cm<x.placebo_em)&(x.placebo_em<0)&(x.emax_cm<x.emax_em)&(x.emax_em<0)&
     (x.baseline_em+x.active_em>0)&(x.baseline_cm+x.active_cm>0))
for name in ['prior_clinical','historical_predictive_clinical','hope_prior_clinical','emax_posterior']:
    x=pd.read_csv(R/f'draws/{name}.csv.gz');check(f'{name}: constraints',clinical_support(x))
for model in ['meta','hope']:
    d=pd.read_csv(R/f'diagnostics/{model}_mcmc.csv')
    check(f'{model}: rank Rhat <= 1.01',d.rhat.max()<1.01)
    check(f'{model}: bulk and tail ESS >= 400',min(d.ess_bulk.min(),d.ess_tail.min())>=400)
cal=pd.read_csv(R/'tables/calibration_quantiles.csv');p=pd.read_csv(R/'draws/prior_clinical.csv.gz')
for _,row in cal.iterrows():
    q=p[row.quantity].quantile([.05,.5,.95]).to_numpy()
    check(row.quantity+': calibration quantiles match validation draws',np.max(abs(q-row[['q05_achieved','q50_achieved','q95_achieved']].to_numpy(dtype=float)))<1e-10)
mcp=pd.read_csv(R/'tables/mcp_calibration.csv').iloc[0]
check('MCP independent null check within 1.96 MCSE of nominal',abs(mcp.independent_null_rejection-.05)<1.96*mcp.null_mcse)
p2=pd.read_csv(R/'draws/phase2_records.csv.gz');p3=pd.read_csv(R/'draws/phase3_records.csv.gz');ev=pd.read_csv(R/'draws/program_events.csv.gz')
qual=p2[['delta_20','delta_100','delta_300']].to_numpy()<=-2
check('Phase 2 Go rule',np.all(p2.go==(p2.significant&qual.any(axis=1))))
expected=np.min(np.where(qual,np.array([20,100,300]),np.inf),axis=1)
check('Smallest qualifying dose selected',np.all(p2.loc[p2.go,'selected_dose']==expected[p2.go]))
actual_p2=p3[p3.strategy=='P2_selected'];check('Selected Phase 3 only after Go',set(actual_p2.program)==set(p2.loc[p2.go,'program']))
for threshold in [0,-2,-3,-4]:
    direct=ev[(ev.threshold==threshold)&(ev.strategy=='P1_direct')].set_index('program')
    screened=ev[(ev.threshold==threshold)&(ev.strategy=='P1_screened')].set_index('program')
    for pop in ['EM','CM','both','both_same_dose']:
        check(f'Paired gate identity {threshold} {pop}',np.all(screened[pop]==(direct[pop]&direct.go)))
for strategy in ev.strategy.unique():
    x=ev[ev.strategy==strategy]
    counts=x.groupby('threshold').both.sum().reindex([-4,-3,-2,0]);check(strategy+': threshold monotonicity',np.all(np.diff(counts)>=0))
    check(strategy+': joint event is EM and CM',np.all(x.both==(x.EM&x.CM)))
    check(strategy+': common-dose implies joint',np.all(~x.both_same_dose|x.both))
for strategy,n in [('P1_fixed',233),('P2_selected',174)]:
    x=p3[p3.strategy==strategy]
    check(strategy+': normal p values',np.max(abs(x.p_two-2*norm.cdf(-abs(x.delta/(5*np.sqrt(2/n))))))<1e-12)
res=p3[(p3.strategy=='P1_fixed')&(p3.population=='EM')].copy();res['noise']=res.delta-res.true_delta
noise=res.pivot(index='program',columns='dose',values='noise');rho=noise.corr().iloc[0,1]
check('Shared placebo noise correlation near 0.5',abs(rho-.5)<.02)
prior=pd.read_csv(R/'diagnostics/prior_only_validation.csv')
for row in prior.itertuples():
    if row.check=='mu_mean':scale=pd.read_csv(R/'tables/calibrated_prior.csv').s[row.dimension-1]
    elif row.check=='omega_pair_mean':scale=1
    else:scale=abs(row.expected)
    check(f'Prior-only {row.check} {row.dimension}',abs(row.observed-row.expected)/scale<.1)
pd.DataFrame(checks).to_csv(R/'diagnostics/reference_checks.csv',index=False)
print(json.dumps({'checks':len(checks),'passed':sum(x['passed'] for x in checks),'R_execution':'not available'},indent=2))

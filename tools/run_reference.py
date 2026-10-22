"""Run the SAME statistical model as the R/Stan package using a separate sampler.

This is a numerical reference execution, not an execution of R or CmdStan.
Requires Python 3, numpy, scipy, pandas, matplotlib and a C++17 compiler.
Run from the bundle root: python MigraineEmaxSim/tools/run_reference.py
All quantities use negative changes. The supplied input is never overwritten.
"""
from pathlib import Path
import os, sys, json, subprocess, hashlib, platform
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
from scipy.special import expit, logit, ndtri
from scipy.optimize import least_squares
from scipy.stats import norm, rankdata
from scipy.fft import rfft, irfft, next_fast_len

PKG = Path(__file__).resolve().parents[1]
ROOT = PKG.parent
OUT = ROOT/'results'
for folder in ['tables','figures','draws','diagnostics']:
    (OUT/folder).mkdir(parents=True,exist_ok=True)
NAMES=['baseline_em','baseline_cm','placebo_em','placebo_cm','emax_em','emax_cm','active_em','active_cm']
TARGET=np.array([[6,10,14.5],[11,15,21],[-2.5,-.5,-.1],[-4,-1,-.25],[-2.5,-.5,-.1],[-4,-1,-.25]])
H=np.array([4/13,4/13,.36,.36,.36,.36])

def clinical(eta):
    u=expit(np.atleast_2d(eta));b=30*u[:,0];c=b+(30-b)*u[:,1]
    a=-b*u[:,2];ac=a-(c+a)*u[:,3]
    p=a*u[:,4];pc=p+(ac-a)*u[:,5]
    return np.column_stack([b,c,p,pc,a-p,ac-pc,a,ac])

def noise(n,seed,identity=False):
    r=np.random.default_rng(seed);a=r.normal(size=(n,6));h=abs(r.normal(size=(n,6)));z=r.normal(size=(n,6))
    L=np.zeros((n,6,6));L[:,0,0]=1
    for i in range(1,6):
        remain=np.ones(n)
        for j in range(i):
            alpha=2+(4-j)/2
            rho=2*r.beta(alpha,alpha,size=n)-1
            L[:,i,j]=rho*remain;remain*=np.sqrt(1-rho*rho)
        L[:,i,i]=remain
    if not identity:z=np.einsum('nij,nj->ni',L,z)
    return a,h,z,L

def quantiles(x):return np.quantile(x,[.05,.5,.95],axis=0).T
def table(x,names=NAMES):
    q=quantiles(x)
    return pd.DataFrame({'quantity':names,'mean':np.mean(x,axis=0),'q05':q[:,0],'q50':q[:,1],'q95':q[:,2]})
def unpack(x):
    t=np.exp(x[6:]);return x[:6],t*np.sqrt(1-H),t*np.sqrt(H)
def pred(x,eps,layer='new_study',factor=1):
    m,s,q=unpack(x);a,h,z,_=eps
    e=np.broadcast_to(m,a.shape).copy()
    if layer in ['new_study','mean']:e+=factor*s*a
    if layer in ['new_study','heterogeneity']:e+=factor*q*h*z
    return clinical(e)

def calibration():
    train=noise(20000,4210)
    m0=logit([1/3,1/4,.1,1/14,.5,.5])
    def fun(x):return np.log(abs(quantiles(pred(x,train))[:6])/abs(TARGET)).ravel()
    runs=[]
    for f in [1,.8,1.2]:
        initial=np.r_[m0,np.log(np.array([.42,.56,1.05,1.45,.65,1.25])*f)]
        runs.append(least_squares(fun,initial,bounds=(np.r_[np.full(6,-9),np.full(6,-4)],np.r_[np.full(6,5),np.full(6,1)]),
            diff_step=.003,max_nfev=400,ftol=1e-8,xtol=1e-8,gtol=1e-8))
    fit=min(runs,key=lambda r:sum(r.fun*r.fun));x=fit.x;m,s,q=unpack(x)
    valid=noise(100000,4211);draws=pred(x,valid);ach=quantiles(draws)[:6]
    comparison=pd.DataFrame({'quantity':NAMES[:6]})
    for j,p in enumerate(['q05','q50','q95']):
        comparison[p+'_target']=TARGET[:,j];comparison[p+'_achieved']=ach[:,j]
    comparison['max_relative_error']=np.max(abs(ach/TARGET-1),axis=1)
    comparison.to_csv(OUT/'tables/calibration_quantiles.csv',index=False)
    pd.DataFrame({'dimension':np.arange(1,7),'m':m,'s':s,'q':q,'heterogeneity_fraction':H}).to_csv(OUT/'tables/calibrated_prior.csv',index=False)
    pd.DataFrame({'quantity':NAMES[:6],'q05':TARGET[:,0],'q50':TARGET[:,1],'q95':TARGET[:,2]}).to_csv(PKG/'inst/extdata/clinical_targets.csv',index=False)
    pd.DataFrame(draws,columns=NAMES).to_csv(OUT/'draws/prior_clinical.csv.gz',index=False)
    pd.DataFrame([{'start':i+1,'loss':sum(r.fun*r.fun),'success':bool(r.success),'message':r.message} for i,r in enumerate(runs)]).to_csv(OUT/'diagnostics/calibration_optimization.csv',index=False)
    sensitivity=[]
    for layer in ['mean','heterogeneity','new_study']:
        t=table(pred(x,valid,layer));t['layer']=layer;sensitivity.append(t)
    pd.concat(sensitivity).to_csv(OUT/'tables/prior_layers.csv',index=False)
    sensitivity=[]
    for f in [.75,1,1.25]:
        p=pred(x,valid,factor=f)
        sensitivity.append({'width_multiplier':f,'p_em_baseline_gt15':np.mean(p[:,0]>15),
          'p_em_emax_le_minus2':np.mean(p[:,4]<=-2),'p_cm_emax_le_minus2':np.mean(p[:,5]<=-2)})
    pd.DataFrame(sensitivity).to_csv(OUT/'tables/prior_width_sensitivity.csv',index=False)
    # A shared-hyperparameter pair isolates differences between studies.
    r=np.random.default_rng(4212);z1=r.normal(size=(100000,6));z2=r.normal(size=(100000,6))
    L=valid[3];tau=q*valid[1];mu=m+s*valid[0]
    c1=clinical(mu+tau*np.einsum('nij,nj->ni',L,z1));c2=clinical(mu+tau*np.einsum('nij,nj->ni',L,z2))
    table(c1-c2).to_csv(OUT/'tables/prior_between_study_differences.csv',index=False)
    omega=L@L.transpose(0,2,1);off=omega[:,np.triu_indices(6,1)[0],np.triu_indices(6,1)[1]]
    assert abs(off.mean())<.002 and abs(off.var()-1/9)<.003
    assert np.linalg.eigvalsh(omega[:200]).min()>0
    pd.DataFrame([{'check':'LKJ pair mean','value':off.mean(),'expected':0},
                  {'check':'LKJ pair variance','value':off.var(),'expected':1/9}]).to_csv(OUT/'diagnostics/lkj_checks.csv',index=False)
    print('Calibration complete; largest relative quantile mismatch:',comparison.max_relative_error.max(),flush=True)
    return x,draws,valid

def arms():
    d=pd.read_csv(PKG/'inst/extdata/supplied_arms.csv');d['input_cm']=d.CM
    for dose,count in [(0,62),(100,32),(750,68)]:
        idx=(d.title=='HOPE')&(d.dose==dose);d.loc[idx,'CM']=count/d.loc[idx,'n']
    d['w_em']=1-d.CM;d['baseline_se']=d.base_mmd_sd/np.sqrt(d.n)
    d['use_baseline']=d.title!='LY3451838'
    d['baseline_note']=np.where(d.use_baseline,'Supplied arm baseline','Omitted: repeated LY baseline needs verification')
    d.to_csv(OUT/'tables/analysis_arms.csv',index=False)
    return d

def run_chains(mode,input_path,warm=15000,keep=30000):
    exe=ROOT/'reference_sampler'
    if not exe.exists():
        subprocess.run(['g++','-O3','-std=c++17',str(PKG/'tools/reference_sampler.cpp'),'-o',str(exe)],check=True)
    paths=[OUT/f'draws/{mode}_chain_{c+1}.csv.gz' for c in range(4)]
    def one(c):
        temp=OUT/f'draws/{mode}_chain_{c+1}.csv'
        subprocess.run([str(exe),mode,str(input_path),str(temp),str(4301+107*c if mode=='meta' else 4401+109*c),str(warm),str(keep),'1'],check=True)
        df=pd.read_csv(temp);df.to_csv(paths[c],index=False,compression='gzip');temp.unlink()
        return df
    with ThreadPoolExecutor(max_workers=4) as pool:dfs=list(pool.map(one,range(4)))
    return np.stack([df.to_numpy() for df in dfs]),dfs[0].columns.tolist()

def split_chains(x):
    n=x.shape[1]//2;return np.concatenate([x[:,:n],x[:,-n:]],axis=0)
def rhat(x):
    n=x.shape[1];W=np.var(x,axis=1,ddof=1).mean(axis=0);B=n*np.var(x.mean(axis=1),axis=0,ddof=1)
    return np.sqrt(((n-1)/n*W+B/n)/W)
def ess(x):
    # Geyer initial positive, monotone paired sequence, with between-chain V+.
    m,n,p=x.shape;center=x-x.mean(axis=1,keepdims=True);size=next_fast_len(2*n)
    f=rfft(center,n=size,axis=1);ac=irfft(f*f.conj(),n=size,axis=1)[:,:n,:]/n
    w=np.var(x,axis=1,ddof=1).mean(axis=0);b=n*np.var(x.mean(axis=1),axis=0,ddof=1)
    vp=(n-1)/n*w+b/n;rho=1-(w[None,:]-ac.mean(axis=0))/vp[None,:];rho[0]=1
    out=[]
    for k in range(p):
        pairs=rho[:2*(n//2),k].reshape(-1,2).sum(axis=1)
        idx=np.flatnonzero(pairs<0);end=idx[0] if len(idx) else len(pairs)
        pairs=np.minimum.accumulate(pairs[:end]);tau=max(1,-1+2*pairs.sum())
        out.append(m*n/tau)
    return np.array(out)
def diagnostics(x,names):
    s=split_chains(x);flat=s.reshape(-1,s.shape[-1]);z=np.empty_like(flat)
    for k in range(flat.shape[-1]):z[:,k]=norm.ppf((rankdata(flat[:,k])-.375)/(len(flat)+.25))
    z=z.reshape(s.shape);fold=abs(s-np.median(flat,axis=0));fz=np.empty_like(flat)
    for k in range(flat.shape[-1]):fz[:,k]=norm.ppf((rankdata(fold[:,:,k].ravel())-.375)/(len(flat)+.25))
    rh=np.maximum(rhat(z),rhat(fz.reshape(s.shape)));bulk=ess(z)
    q=np.quantile(flat,[.05,.95],axis=0);tail=np.minimum(ess((s<=q[0]).astype(float)),ess((s<=q[1]).astype(float)))
    return pd.DataFrame({'parameter':names,'rhat':rh,'ess_bulk':bulk,'ess_tail':tail})

def fit_meta(d,x):
    m,s,q=unpack(x);h=d[d.title!='HOPE'].copy();b=h[h.use_baseline].copy();ids={name:i for i,name in enumerate(h.nct.unique())}
    text=' '.join(map(str,np.r_[m,s,q]))+'\n'+f'2 {len(ids)} {len(h)} {len(b)}\n'
    for r in h.itertuples():text+=f'{ids[r.nct]} {r.chg} {r.chg_se} {r.w_em} {int(r.arm=="active")} {.5 if r.arm=="active" else 1}\n'
    for r in b.itertuples():text+=f'{ids[r.nct]} {r.base_mmd} {r.baseline_se} {r.w_em}\n'
    path=OUT/'diagnostics/meta_input.txt';path.write_text(text)
    chains,names=run_chains('meta',path);diag=diagnostics(chains[:,:,:-6],names[:-6]);diag.to_csv(OUT/'diagnostics/meta_mcmc.csv',index=False)
    print('Historical model: max Rhat',diag.rhat.max(),'minimum bulk ESS',diag.ess_bulk.min(),flush=True)
    flat=chains.reshape(-1,chains.shape[-1]);eta=flat[:,-6:];pmean=eta.mean(axis=0);pcov=np.cov(eta.T)+np.eye(6)*1e-8
    np.savetxt(OUT/'tables/predictive_eta_mean.csv',pmean[None,:],delimiter=',')
    np.savetxt(OUT/'tables/predictive_eta_cov.csv',pcov,delimiter=',')
    pd.DataFrame(clinical(eta),columns=NAMES).to_csv(OUT/'draws/historical_predictive_clinical.csv.gz',index=False)
    table(clinical(eta)).to_csv(OUT/'tables/historical_predictive_summary.csv',index=False)
    study=[]
    for j,name in enumerate(h.title.unique()):
        e=flat[:,27+j*6:33+j*6];t=table(clinical(e));t['study']=name;study.append(t)
    pd.concat(study).to_csv(OUT/'tables/historical_study_summary.csv',index=False)
    # Both underlying arm means and observational posterior replications.
    r=np.random.default_rng(4302);ppc=[]
    for kind,frame in [('change',h),('baseline',b)]:
        for idx,row in frame.iterrows():
            j=ids[row.nct];c=clinical(flat[:,27+j*6:33+j*6])
            if kind=='change':
                a=int(row.arm=='active');mu=row.w_em*(c[:,2]+a*c[:,4])+(1-row.w_em)*(c[:,3]+a*c[:,5]);obs=row.chg;se=row.chg_se
            else:mu=row.w_em*c[:,0]+(1-row.w_em)*c[:,1];obs=row.base_mmd;se=row.baseline_se
            rep=mu+r.normal(size=len(mu))*se
            ppc.append({'kind':kind,'row':idx,'study':row.title,'arm':row.arm,'w_em':row.w_em,
             'observed':obs,'se':se,'mu05':np.quantile(mu,.05),'mu50':np.median(mu),'mu95':np.quantile(mu,.95),
             'rep05':np.quantile(rep,.05),'rep50':np.median(rep),'rep95':np.quantile(rep,.95)})
    pd.DataFrame(ppc).to_csv(OUT/'tables/historical_arm_ppc.csv',index=False)
    return chains,names,pmean,pcov

def fit_hope(d,pmean,pcov):
    h=d[d.title=='HOPE'];b=h[h.use_baseline];em=np.log(20);es=np.log(10)/(2*norm.ppf(.95))
    text=' '.join(map(str,pmean))+'\n'+' '.join(map(str,pcov.ravel()))+f'\n{em} {es} {len(h)} {len(b)}\n'
    for row in h.itertuples():text+=f'{row.chg} {row.chg_se} {row.w_em} {row.dose}\n'
    for row in b.itertuples():text+=f'{row.base_mmd} {row.baseline_se} {row.w_em}\n'
    path=OUT/'diagnostics/hope_input.txt';path.write_text(text)
    chains,names=run_chains('hope',path,warm=8000,keep=20000);diag=diagnostics(chains,names);diag.to_csv(OUT/'diagnostics/hope_mcmc.csv',index=False)
    print('HOPE model: max Rhat',diag.rhat.max(),'minimum bulk ESS',diag.ess_bulk.min(),flush=True)
    flat=chains.reshape(-1,7);clinical_draws=clinical(flat[:,:6]);full=np.column_stack([clinical_draws,flat[:,6]])
    pd.DataFrame(full,columns=NAMES+['ed50']).to_csv(OUT/'draws/emax_posterior.csv.gz',index=False)
    table(full,NAMES+['ed50']).to_csv(OUT/'tables/emax_posterior_summary.csv',index=False)
    r=np.random.default_rng(4402);gp=r.multivariate_normal(pmean,pcov,100000);gaussian=clinical(gp)
    pd.DataFrame(gaussian,columns=NAMES).to_csv(OUT/'draws/hope_prior_clinical.csv.gz',index=False)
    table(gaussian).to_csv(OUT/'tables/hope_prior_summary.csv',index=False)
    edprior=r.lognormal(em,es,100000)
    pd.DataFrame([{'logmean':em,'logsd':es,'q05':np.exp(em+es*norm.ppf(.05)),
                  'q50':20,'q95':np.exp(em+es*norm.ppf(.95))}]).to_csv(OUT/'tables/ed50_prior.csv',index=False)
    ppc=[]
    for idx,row in h.iterrows():
        c=clinical_draws;g=row.dose/(row.dose+flat[:,6]);mu=row.w_em*(c[:,2]+g*c[:,4])+(1-row.w_em)*(c[:,3]+g*c[:,5])
        rep=mu+r.normal(size=len(mu))*row.chg_se
        ppc.append({'dose':row.dose,'observed':row.chg,'se':row.chg_se,'mu05':np.quantile(mu,.05),
          'mu50':np.median(mu),'mu95':np.quantile(mu,.95),'rep05':np.quantile(rep,.05),'rep95':np.quantile(rep,.95)})
    pd.DataFrame(ppc).to_csv(OUT/'tables/hope_arm_ppc.csv',index=False)
    return chains,full,gaussian,edprior

def prior_only_validation(x):
    # An independent prior-only chain checks the hierarchical sampler target.
    m,s,q=unpack(x);path=OUT/'diagnostics/prior_only_input.txt'
    path.write_text(' '.join(map(str,np.r_[m,s,q]))+'\n2 3 0 0\n')
    raw=OUT/'diagnostics/prior_only.csv'
    subprocess.run([str(ROOT/'reference_sampler'),'meta',str(path),str(raw),'8891','10000','50000','1'],check=True)
    df=pd.read_csv(raw);v=df.iloc[:,12:27].to_numpy();L=np.zeros((len(v),6,6));L[:,0,0]=1;k=0
    for i in range(1,6):
        remain=np.ones(len(v))
        for j in range(i):
            rho=np.tanh(v[:,k]);k+=1;L[:,i,j]=rho*remain;remain*=np.sqrt(1-rho*rho)
        L[:,i,i]=remain
    omega=L@L.transpose(0,2,1);off=omega[:,np.triu_indices(6,1)[0],np.triu_indices(6,1)[1]]
    rows=[]
    for i in range(6):
        for name,actual,expected in [('mu_mean',df[f'mu{i}'].mean(),m[i]),('mu_sd',df[f'mu{i}'].std(),s[i]),
          ('tau_mean',df[f'tau{i}'].mean(),q[i]*np.sqrt(2/np.pi)),('eta_new_variance',df[f'eta_new{i}'].var(),s[i]**2+q[i]**2)]:
            rows.append({'check':name,'dimension':i+1,'observed':actual,'expected':expected})
    rows.extend([{'check':'omega_pair_mean','dimension':0,'observed':off.mean(),'expected':0},
                 {'check':'omega_pair_variance','dimension':0,'observed':off.var(),'expected':1/9}])
    pd.DataFrame(rows).to_csv(OUT/'diagnostics/prior_only_validation.csv',index=False);raw.unlink()

def probability(success,eligible=None):
    if eligible is not None:success=success[eligible]
    n=len(success);k=int(np.sum(success));p=k/n if n else np.nan;z=norm.ppf(.975)
    if n:
        center=(p+z*z/(2*n))/(1+z*z/n);half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    else:center=half=np.nan
    return dict(n=n,successes=k,probability=p,mcse=np.sqrt(p*(1-p)/n) if n else np.nan,lower=center-half,upper=center+half)

def programs(full,n=50000):
    r=np.random.default_rng(4601);ids=r.integers(len(full),size=n);truth=full[ids];doses=np.array([0,20,100,300])
    candidates=np.quantile(full[:,-1],[.1,.5,.9]);g=doses[:,None]/(doses[:,None]+candidates)
    contrasts=-(g-g.mean(axis=0));contrasts/=np.sqrt((contrasts**2).sum(axis=0))
    nr=np.random.default_rng(4501);null=np.max(nr.normal(size=(200000,4))@contrasts,axis=1);critical=np.quantile(null,.95)
    # Independent null check, not reusing the calibration sample.
    nr=np.random.default_rng(4502);rate=np.mean(np.max(nr.normal(size=(100000,4))@contrasts,axis=1)>critical)
    pd.DataFrame([{'critical_value':critical,'ed50_10':candidates[0],'ed50_50':candidates[1],'ed50_90':candidates[2],
                  'independent_null_rejection':rate,'null_mcse':np.sqrt(rate*(1-rate)/100000)}]).to_csv(OUT/'tables/mcp_calibration.csv',index=False)
    def means(ds,pop):
        return truth[:,2+pop,None]+truth[:,4+pop,None]*ds/(ds+truth[:,-1,None])
    mu=(means(doses,0)+means(doses,1))/2
    y=mu+r.normal(size=mu.shape)*5/np.sqrt(140)
    delta=y-y[:,[0]];sig=np.max(y@contrasts/(5/np.sqrt(140)),axis=1)>critical
    qualifying=(delta<=-2)&(doses>0);go=sig&qualifying.any(axis=1)
    selected=np.min(np.where(qualifying,doses,np.inf),axis=1);selected[~go]=np.nan
    stage2=pd.DataFrame({'program':np.arange(1,n+1),'posterior_draw':ids+1,'go':go,'selected_dose':selected,'significant':sig})
    for i,dose in enumerate(doses):stage2[f'delta_{dose}']=delta[:,i]
    stage2.to_csv(OUT/'draws/phase2_records.csv.gz',index=False)
    direct=[];chosen=[];ph3=[]
    for pop in range(2):
        ds=np.array([0,20,200]);mu=means(ds,pop);y=mu+r.normal(size=mu.shape)*5/np.sqrt(233)
        dh=y[:,1:]-y[:,[0]];td=mu[:,1:]-mu[:,[0]];pv=2*norm.cdf(-abs(dh/(5*np.sqrt(2/233))));direct.append((dh,pv,td))
        for i,dose in enumerate(ds[1:]):ph3.append(pd.DataFrame({'program':np.arange(1,n+1),'strategy':'P1_fixed','population':['EM','CM'][pop],
              'dose':dose,'delta':dh[:,i],'p_two':pv[:,i],'true_delta':td[:,i]}))
        # Generate selected-dose phase 3 observations ONLY in advancing programs.
        mu0=truth[go,2+pop];mua=mu0+truth[go,4+pop]*selected[go]/(selected[go]+truth[go,-1])
        y0=mu0+r.normal(size=go.sum())*5/np.sqrt(174);ya=mua+r.normal(size=go.sum())*5/np.sqrt(174)
        dh2=np.full(n,np.nan);pv2=np.full(n,np.nan);dh2[go]=ya-y0;pv2[go]=2*norm.cdf(-abs((ya-y0)/(5*np.sqrt(2/174))))
        chosen.append((dh2,pv2));ph3.append(pd.DataFrame({'program':np.arange(1,n+1)[go],'strategy':'P2_selected','population':['EM','CM'][pop],
              'dose':selected[go],'delta':dh2[go],'p_two':pv2[go],'true_delta':mua-mu0}))
    pd.concat(ph3).to_csv(OUT/'draws/phase3_records.csv.gz',index=False)
    rows=[];derisk=[];eventframes=[]
    for threshold in [0,-2,-3,-4]:
        sd=[(a[0]<0)&(a[0]<=threshold)&(a[1]<.05) for a in direct]
        ss=[(a[0]<0)&(a[0]<=threshold)&(a[1]<.05) for a in chosen]
        s1={'EM':sd[0].any(axis=1),'CM':sd[1].any(axis=1),'both':sd[0].any(axis=1)&sd[1].any(axis=1),'both_same_dose':(sd[0]&sd[1]).any(axis=1)}
        s2={'EM':ss[0],'CM':ss[1],'both':ss[0]&ss[1],'both_same_dose':ss[0]&ss[1]}
        for strategy,events in [('P1_direct',s1),('P1_screened',{k:v&go for k,v in s1.items()}),('P2_selected',s2)]:
            eventframes.append(pd.DataFrame({'program':np.arange(1,n+1),'strategy':strategy,'threshold':threshold,'go':go,**events}))
            for pop,s in events.items():
                for estimand,eligible in [('from_start',None),('given_go',go)]:
                    rows.append(dict(strategy=strategy,threshold=threshold,population=pop,estimand=estimand,**probability(s,eligible)))
        for pop,s in s1.items():
            g=go.astype(float);s=s.astype(float);pgo=g.mean();pa=s.mean();ps=s[go].mean();influence=g*(s-ps)/pgo-(s-pa)
            derisk.append(dict(threshold=threshold,population=pop,p_go=pgo,direct=pa,given_go=ps,given_no_go=s[~go].mean(),
              gated_from_start=np.mean(g*s),conditional_gain=ps-pa,gain_mcse=influence.std(ddof=1)/np.sqrt(n)))
    pd.DataFrame(rows).to_csv(OUT/'tables/program_probabilities.csv',index=False)
    pd.DataFrame(derisk).to_csv(OUT/'tables/fixed_design_derisking.csv',index=False)
    pd.concat(eventframes).to_csv(OUT/'draws/program_events.csv.gz',index=False)
    pd.DataFrame([probability(go)]).to_csv(OUT/'tables/phase2_go.csv',index=False)
    pd.DataFrame({'dose':[20,100,300],'count':[np.sum(selected==d) for d in [20,100,300]],
                  'fraction_of_go':[np.mean(selected[go]==d) for d in [20,100,300]]}).to_csv(OUT/'tables/phase2_selection.csv',index=False)
    pd.DataFrame({'strategy':['P1_direct','P1_screened','P2_selected'],'expected_n':[1398,560+1398*go.mean(),560+696*go.mean()]}).to_csv(OUT/'tables/expected_participants.csv',index=False)
    truth_df=pd.DataFrame(truth,columns=NAMES+['ed50']);truth_df.insert(0,'posterior_draw',ids+1);truth_df.insert(0,'program',np.arange(1,n+1));truth_df.to_csv(OUT/'draws/program_truth.csv.gz',index=False)
    # Shared event algebra is an essential check on the interpretation of de-risking.
    assert all(np.mean(s&go)<=np.mean(s)+1e-12 for s in s1.values())
    print('Programs complete. Go:',go.mean(),flush=True)

def main():
    x,prior,valid=calibration();d=arms();meta,meta_names,pmean,pcov=fit_meta(d,x)
    prior_only_validation(x)
    hope,full,gaussian,edprior=fit_hope(d,pmean,pcov);programs(full)
    from make_figures import figures
    figures(x,prior,valid,d,meta,meta_names,hope,full,gaussian,edprior)
    inputs=hashlib.sha256((PKG/'inst/extdata/supplied_arms.csv').read_bytes()).hexdigest()
    (OUT/'execution_manifest.json').write_text(json.dumps({'execution':'Independent C++/Python reference; R and CmdStan not run',
       'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,
       'supplied_input_sha256':inputs,'calibration_training':20000,'calibration_validation':100000,
       'meta_chains':4,'meta_warmup_per_chain':15000,'meta_retained_per_chain':30000,
       'hope_chains':4,'hope_warmup_per_chain':8000,'hope_retained_per_chain':20000,
       'programs':50000,'warning':'Quantile targets and variance allocation remain provisional assumptions.'},indent=2))
    print('Reference analysis finished.',flush=True)
if __name__=='__main__':main()

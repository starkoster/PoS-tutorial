"""Publication figures for the transparently labeled numerical reference run."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'results';F=OUT/'figures'
BLUE='#001965';LIGHT='#6F8DCC';ORANGE='#C66A2E';GREY='#8E99AB'
NAMES=['baseline_em','baseline_cm','placebo_em','placebo_cm','emax_em','emax_cm','active_em','active_cm']
LABELS=['EM baseline','CM baseline','EM placebo change','CM placebo change','EM maximum drug effect','CM maximum drug effect','EM active change','CM active change']
plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','axes.spines.top':False,
 'axes.spines.right':False,'axes.titlecolor':BLUE,'axes.labelcolor':'#333333',
 'axes.titleweight':'bold','pdf.fonttype':42,'ps.fonttype':42,'legend.frameon':False})

def save(fig,name):
    fig.savefig(F/f'{name}.pdf',bbox_inches='tight');fig.savefig(F/f'{name}.png',dpi=180,bbox_inches='tight');plt.close(fig)
def density(ax,x,color,label=None,limits=None):
    # Subsample for plotting only; summaries use every draw.
    x=x[np.linspace(0,len(x)-1,min(10000,len(x))).astype(int)]
    limits=limits if limits is not None else np.quantile(x,[.001,.999])
    grid=np.linspace(*limits,300);ax.plot(grid,gaussian_kde(x)(grid),color=color,lw=1.5,label=label)
def marginals(sets,name,title,targets=None):
    fig,axs=plt.subplots(2,3,figsize=(10.5,6),layout='constrained')
    for i,ax in enumerate(axs.flat):
        vals=np.concatenate([x[:,i] for _,x,_ in sets]);lo,hi=np.quantile(vals,[.001,.999])
        if targets is not None:lo=min(lo,min(targets[i]));hi=max(hi,max(targets[i]))
        for label,x,color in sets:density(ax,x[:,i],color,label,limits=(lo,hi))
        if targets is not None:
            for j,v in enumerate(targets[i]):ax.axvline(v,color=ORANGE,ls='--' if j==1 else ':',lw=1)
        ax.set(title=LABELS[i],xlabel='MMDs',ylabel='Density',xlim=(lo,hi))
    axs[0,0].legend(fontsize=9);fig.suptitle(title,color=BLUE,fontsize=14,fontweight='bold');save(fig,name)

def figures(x,prior,valid,d,meta,meta_names,hope,full,gaussian,edprior):
    # Imported here so the execution script can also be run as a module.
    from run_reference import clinical, pred, noise, TARGET
    F.mkdir(parents=True,exist_ok=True)
    marginals([('Calibrated prior',prior,BLUE)],'01_prior_calibration','Clinical targets and calibrated new-study prior',TARGET)
    marginals([('Mean uncertainty',pred(x,valid,'mean'),GREY),
               ('Heterogeneity at fixed mean',pred(x,valid,'heterogeneity'),LIGHT),
               ('Full new-study prior',prior,BLUE)],'03_prior_layers','Separate contributions to clinical uncertainty')
    initial_identity=pred(x,noise(50000,4214,True))
    crossplots([('Independent latent deviations',initial_identity),('LKJ(2) integrated',prior)],'02_prior_em_cm')
    mf=meta.reshape(-1,meta.shape[-1]);hp=clinical(mf[:,-6:])
    marginals([('Initial new-study prior',prior,GREY),('Historical new-study prediction',hp,BLUE)],
              '04_meta_marginals','Historical update on the clinical scale')
    # Six marginal panels include all historical studies and initial/posterior
    # predictive intervals. Study points are latent study parameter posteriors.
    studies=[clinical(mf[:,27+6*j:33+6*j]) for j in range(3)]
    sets=[prior,*studies,hp];names=['Initial new-study prior','CONQUER posterior','DELIVER posterior','LY3451838 posterior','New-study posterior prediction']
    fig,axs=plt.subplots(2,3,figsize=(11,6.5),layout='constrained')
    for i,ax in enumerate(axs.flat):
        for row,xx in enumerate(sets):
            q=np.quantile(xx[:,i],[.05,.5,.95]);ax.plot(q[[0,2]],[row,row],color=GREY if row==0 else BLUE,lw=1.5)
            ax.scatter(q[1],row,color=ORANGE if row in [0,4] else BLUE,s=22,zorder=3)
        ax.set(yticks=range(5),yticklabels=names if i%3==0 else ['']*5,title=LABELS[i],xlabel='MMDs: median and central 90% interval')
        ax.invert_yaxis()
    fig.suptitle('Historical studies, priors and posterior predictions',color=BLUE,fontsize=14,fontweight='bold');save(fig,'05_meta_studies')
    crossplots([('Initial new-study prior',prior),('Historical new-study prediction',hp)],'06_meta_em_cm',studies)
    marginals([('Empirical predictive distribution',hp,BLUE),('Gaussian transfer approximation',gaussian,ORANGE)],
              '07_transfer_check','Check of the Gaussian transfer approximation')
    marginals([('Transferred HOPE prior',gaussian,GREY),('HOPE posterior',full,BLUE)],
              '08_hope_update','Same-drug update with HOPE')
    for adjusted,name in [(True,'09_emax_adjusted'),(False,'10_emax_change')]:
        fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained');ds=np.linspace(0,750,151)
        for pop,ax in enumerate(axs):
            curve=full[:,4+pop,None]*ds/(ds+full[:,-1,None])
            if not adjusted:curve+=full[:,2+pop,None]
            q=np.quantile(curve,[.05,.5,.95],axis=0)
            ax.fill_between(ds,q[0],q[2],color=LIGHT,alpha=.35,label='Central 90% interval');ax.plot(ds,q[1],color=BLUE,label='Posterior median')
            ax.axhline(-2,color=ORANGE,ls=':',lw=1);ax.set(title=['Episodic migraine','Chronic migraine'][pop],xlabel='Dose (mg)',
                ylabel='Active minus placebo (MMDs)' if adjusted else 'Change from baseline (MMDs)')
        axs[0].legend(fontsize=8);save(fig,name)
    fig,ax=plt.subplots(figsize=(6,3.5),layout='constrained')
    grid=(0,100);density(ax,edprior,GREY,'Separate lognormal prior',grid);density(ax,full[:,-1],BLUE,'HOPE posterior',grid)
    ax.set(xlabel='ED50 (mg)',ylabel='Density',title='Compound-specific potency');ax.legend();save(fig,'11_ed50_update')
    historical_input_figure(d)
    arm_ppc_figure()
    program_figures()
    # Trace plots for selected difficult parameters, not decorative diagnostics.
    for arr,names,name in [(meta[:,:,:12],meta_names[:12],'18_meta_traces'),(hope,['eta1','eta2','eta3','eta4','eta5','eta6','ED50'],'19_hope_traces')]:
        fig,axs=plt.subplots(len(names),1,figsize=(10,max(7,len(names)*1.1)),sharex=True,layout='constrained')
        for k,ax in enumerate(np.atleast_1d(axs)):
            for chain in range(4):ax.plot(arr[chain,::40,k],lw=.5,alpha=.7)
            ax.set_ylabel(names[k],rotation=0,labelpad=25,fontsize=8)
        np.atleast_1d(axs)[-1].set_xlabel('Retained iteration / 40 (four chains)');save(fig,name)

def crossplots(sets,name,studies=None):
    fig,axs=plt.subplots(len(sets),3,figsize=(10,3.2*len(sets)),layout='constrained',squeeze=False)
    pairs=[(0,1,'Baseline'),(2,3,'Placebo change'),(6,7,'Active change')]
    for row,(stage,draws) in enumerate(sets):
        for col,(a,b,label) in enumerate(pairs):
            ax=axs[row,col];idx=np.linspace(0,len(draws)-1,12000).astype(int);p=draws[idx]
            ax.hexbin(p[:,a],p[:,b],gridsize=45,mincnt=1,bins='log',cmap='Blues',rasterized=True)
            lo=min(np.quantile(p[:,a],.001),np.quantile(p[:,b],.001));hi=max(np.quantile(p[:,a],.999),np.quantile(p[:,b],.999))
            ax.plot([lo,hi],[lo,hi],color=GREY,ls=':',lw=1)
            if studies is not None and row==len(sets)-1:
                for j,s in enumerate(studies):ax.scatter(np.median(s[:,a]),np.median(s[:,b]),color=ORANGE,marker=['o','s','^'][j],s=35,label=['CONQUER','DELIVER','LY3451838'][j])
            ax.set(xlabel='EM (MMDs)',ylabel='CM (MMDs)',title=f'{stage}\n{label}; r = {np.corrcoef(draws[:,a],draws[:,b])[0,1]:.2f}',xlim=(lo,hi),ylim=(lo,hi))
    if studies is not None:axs[-1,0].legend(fontsize=7,loc='upper left')
    save(fig,name)

def historical_input_figure(d):
    h=d[d.title!='HOPE'];fig,axs=plt.subplots(1,3,figsize=(11,4.5),layout='constrained')
    for ax,kind,title in zip(axs,['baseline','placebo','active'],['Baseline means','Placebo changes','Active changes']):
        sub=h[h.use_baseline] if kind=='baseline' else h[h.arm==kind]
        labels=[]
        for i,row in enumerate(sub.itertuples()):
            mu=row.base_mmd if kind=='baseline' else row.chg;se=row.baseline_se if kind=='baseline' else row.chg_se
            ax.errorbar(mu,i,xerr=1.96*se,fmt='o',color=BLUE,ms=4,capsize=2)
            labels.append(f'{row.title} / '+('EM' if row.CM==0 else 'CM' if row.CM==1 else 'mixed')+('/'+row.arm if kind=='baseline' else ''))
        ax.set(yticks=range(len(sub)),yticklabels=labels,title=title,xlabel='Observed mean and 95% interval (MMDs)');ax.invert_yaxis()
        ax.tick_params(axis='y',labelsize=7)
    save(fig,'12_historical_inputs')

def arm_ppc_figure():
    p=pd.read_csv(OUT/'tables/historical_arm_ppc.csv');fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for ax,kind in zip(axs,['baseline','change']):
        sub=p[p.kind==kind];labels=[]
        for i,row in enumerate(sub.itertuples()):
            ax.plot([row.rep05,row.rep95],[i,i],color=LIGHT,lw=3,alpha=.6)
            ax.plot([row.mu05,row.mu95],[i,i],color=BLUE,lw=1.5)
            ax.scatter(row.observed,i,color=ORANGE,s=25,zorder=4)
            labels.append(f'{row.study} / {row.arm} / '+('EM' if row.w_em==1 else 'CM' if row.w_em==0 else 'mixed'))
        ax.set(yticks=range(len(sub)),yticklabels=labels,title='Baseline' if kind=='baseline' else 'Change from baseline',xlabel='MMDs');ax.invert_yaxis();ax.tick_params(axis='y',labelsize=8)
    fig.suptitle('Arm checks: observation (orange), mean interval (blue), replication interval (light blue)',color=BLUE,fontsize=11)
    save(fig,'13_historical_ppc')

def program_figures():
    p=pd.read_csv(OUT/'tables/program_probabilities.csv');fig,axs=plt.subplots(1,3,figsize=(11,4),layout='constrained')
    strategies=['P1_direct','P1_screened','P2_selected'];labels=['Direct P3','P2 gate + fixed P3','P2 + selected P3']
    for ax,pop in zip(axs,['EM','CM','both']):
        for j,strategy in enumerate(strategies):
            sub=p[(p.population==pop)&(p.strategy==strategy)&(p.estimand=='from_start')].set_index('threshold').loc[[0,-2,-3,-4]]
            y=sub.probability.to_numpy();x=np.arange(4)+(j-1)*.25
            ax.bar(x,y,width=.23,color=[BLUE,LIGHT,ORANGE][j],label=labels[j]);ax.errorbar(x,y,yerr=[y-sub.lower,sub.upper-y],color='#333333',fmt='none',lw=.7,capsize=2)
        ax.set(title={'EM':'EM trial','CM':'CM trial','both':'Both trials'}[pop],xticks=range(4),xticklabels=['Stat.','-2','-3','-4'],ylim=(0,1),xlabel='Observed delta threshold (MMDs)');ax.yaxis.set_major_formatter(PercentFormatter(1))
    axs[0].set_ylabel('Success from program start');axs[1].legend(loc='upper center',bbox_to_anchor=(.5,-.2),ncol=3,fontsize=8);save(fig,'14_program_comparison')
    q=pd.read_csv(OUT/'tables/fixed_design_derisking.csv');fig,axs=plt.subplots(1,3,figsize=(11,4),layout='constrained')
    for ax,pop in zip(axs,['EM','CM','both']):
        sub=q[q.population==pop].set_index('threshold').loc[[0,-2,-3,-4]]
        for col,label,color in [('direct','Direct P3',GREY),('given_go','P3 given Go',BLUE),('gated_from_start','Gate and P3 success',ORANGE)]:
            ax.plot(range(4),sub[col],marker='o',color=color,label=label)
        ax.set(title=pop,xticks=range(4),xticklabels=['Stat.','-2','-3','-4'],ylim=(0,1),xlabel='Observed delta threshold (MMDs)');ax.yaxis.set_major_formatter(PercentFormatter(1))
    axs[0].set_ylabel('Probability');axs[1].legend(loc='upper center',bbox_to_anchor=(.5,-.2),ncol=3,fontsize=8);save(fig,'15_fixed_design_derisking')
    s=pd.read_csv(OUT/'tables/phase2_selection.csv');fig,ax=plt.subplots(figsize=(6,3.5),layout='constrained')
    ax.bar(s.dose.astype(str),s.fraction_of_go,color=[BLUE,LIGHT,ORANGE]);ax.yaxis.set_major_formatter(PercentFormatter(1));ax.set(xlabel='Selected dose (mg)',ylabel='Fraction among Go programs',title='Phase 2 dose selection');save(fig,'16_phase2_selection')
    ph=pd.read_csv(OUT/'draws/phase3_records.csv.gz');fig,axs=plt.subplots(1,3,figsize=(10.5,3.8),layout='constrained')
    subsets=[('P1_fixed',20,'Direct P3: 20 mg'),('P1_fixed',200,'Direct P3: 200 mg'),('P2_selected',None,'Selected-dose P3 after Go')]
    for ax,(strategy,dose,label) in zip(axs,subsets):
        sub=ph[ph.strategy==strategy]
        if dose is not None:sub=sub[sub.dose==dose]
        pairs=sub.pivot(index='program',columns='population',values='delta');ax.hexbin(pairs.EM,pairs.CM,gridsize=45,bins='log',mincnt=1,cmap='Blues',rasterized=True)
        ax.axvline(-2,color=ORANGE,ls=':');ax.axhline(-2,color=ORANGE,ls=':');ax.set(title=label,xlabel='Observed EM delta (MMDs)',ylabel='Observed CM delta (MMDs)')
    save(fig,'17_phase3_joint_effects')

if __name__=='__main__':
    import run_reference as rr
    p=pd.read_csv(OUT/'tables/calibrated_prior.csv');x=np.r_[p.m,np.log(np.sqrt(p.s**2+p.q**2))]
    prior=pd.read_csv(OUT/'draws/prior_clinical.csv.gz').to_numpy();valid=rr.noise(100000,4211)
    meta=np.stack([pd.read_csv(OUT/f'draws/meta_chain_{i+1}.csv.gz').to_numpy() for i in range(4)])
    names=pd.read_csv(OUT/'draws/meta_chain_1.csv.gz',nrows=1).columns.tolist()
    hope=np.stack([pd.read_csv(OUT/f'draws/hope_chain_{i+1}.csv.gz').to_numpy() for i in range(4)])
    full=pd.read_csv(OUT/'draws/emax_posterior.csv.gz').to_numpy();gaussian=pd.read_csv(OUT/'draws/hope_prior_clinical.csv.gz').to_numpy()
    edprior=np.random.default_rng(4403).lognormal(np.log(20),np.log(10)/(2*rr.norm.ppf(.95)),100000)
    figures(x,prior,valid,pd.read_csv(OUT/'tables/analysis_arms.csv'),meta,names,hope,full,gaussian,edprior)

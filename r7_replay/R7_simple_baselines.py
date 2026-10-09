#!/usr/bin/env python3
"""R7 science audit: run ordinary held-out interpolation on exactly the R3
nine steel-contact loops and the R3 10 seeded index splits; no test-response
information enters input generation or fitting. Scores use published R3 metric.
Same-loop point reconstruction ONLY; not prospective phase availability.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator, CubicSpline
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parent
DATA=ROOT/'R3'/'03_Code_Replay'/'inputs'
R3CSV=ROOT/'R3'/'02_Results_and_Methods'/'R3_target_blind_10seed_per_split.csv'
END={"AF10":1522500,"AF19":997500,"AF24":1575000,"AF28":1575000,"AF29":1575000,"AF36":1575000,"AF38":1575000,"AF39":1575000,"AF40":1575000}
MODELS=('Linear (clamped)','PCHIP (clamped)','Cubic spline (clamped)')

def predict(y,u,tr,te,method):
    ix=np.argsort(u[tr]); xtrain=u[tr][ix]; vtrain=y[tr][ix]
    xmin,xmax=xtrain[0],xtrain[-1]
    t=np.clip(u[te],xmin,xmax)  # constant extrapolation; no observed held-out boundary used
    if method=='Linear (clamped)': return np.interp(t,xtrain,vtrain)
    if method=='PCHIP (clamped)': return PchipInterpolator(xtrain,vtrain,extrapolate=False)(t)
    if method=='Cubic spline (clamped)': return CubicSpline(xtrain,vtrain,bc_type='natural',extrapolate=False)(t)
    raise ValueError(method)

def met(x,f,px,pf,te):
    ex=(px-x[te])/np.ptp(x); ef=(pf-f[te])/np.ptp(f)
    return {'nrmse_pair':float(np.sqrt(np.mean((ex*ex+ef*ef)/2))),
            'nrmse_force':float(np.sqrt(np.mean(ef*ef))),
            'nmae_force':float(np.mean(np.abs(ef))),
            'r2_force':float(1-np.sum((pf-f[te])**2)/np.sum((f[te]-np.mean(f[te]))**2))}

def main():
    orig=pd.read_csv(R3CSV);orig=orig[(orig.method=='index')]
    assert orig.shape[0]==360, orig.shape
    assert set(orig.experiment)==set(END),set(orig.experiment)
    assert set(orig.seed)==set(range(10))
    rows=[]; provenance=[]
    for exp,cycle in sorted(END.items()):
        p=DATA/f'{exp}_{cycle}cycle.mat'
        a=np.asarray(loadmat(p,squeeze_me=True)['hyst'],dtype=float)
        assert a.shape==(571,2) and np.isfinite(a).all()
        x,f=a.T; u=np.linspace(0,1,len(x))
        provenance.append({'experiment':exp,'file':p.name,'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'n':len(x)})
        for seed in range(10):
            rng=np.random.default_rng(seed);perm=rng.permutation(len(u));ntr=int(round(.7*len(u)))
            tr,te=perm[:ntr],perm[ntr:]
            assert len(tr)==400 and len(te)==171
            edge=np.count_nonzero((u[te]<u[tr].min()) | (u[te]>u[tr].max()))
            for model in MODELS:
                px=predict(x,u,tr,te,model);pf=predict(f,u,tr,te,model)
                rows.append({'experiment':exp,'method':'index','seed':seed,'model':model,'n_train':len(tr),'n_test':len(te),'n_endpoint_heldout':int(edge),**met(x,f,px,pf,te)})
    df=pd.DataFrame(rows);df.to_csv(ROOT/'R7_simple_baseline_per_split.csv',index=False)
    orig=orig.copy();orig['n_endpoint_heldout']=np.nan
    combined=pd.concat([orig,df],ignore_index=True)
    assert combined.shape[0]==630
    combined.to_csv(ROOT/'R7_all_7_models_per_split.csv',index=False)
    condition=combined.groupby(['experiment','model'],as_index=False)[['nrmse_pair','nrmse_force']].mean()
    condition.to_csv(ROOT/'R7_seven_models_per_condition.csv',index=False)
    agg=condition.groupby('model',as_index=False)[['nrmse_pair','nrmse_force']].mean()
    agg=agg.sort_values('nrmse_pair');agg.to_csv(ROOT/'R7_seven_models_aggregate.csv',index=False)
    # Bootstrap among NINE fixed experimental conditions, paired case resampling
    names=['FCM-TS','SC9-TS','SC16-TS','ANFIS-TS',*MODELS]
    conditions=sorted(END)
    metrics=['nrmse_pair','nrmse_force']; rg=np.random.default_rng(20261008)
    B=50000; draws=rg.integers(0,len(conditions),size=(B,len(conditions)))
    rows=[]
    for metric in metrics:
        mat=np.array([[float(condition.query('experiment==@e and model==@m')[metric].iloc[0]) for m in names] for e in conditions])
        assert mat.shape==(9,7)
        means=mat[draws].mean(axis=1)
        for k,m in enumerate(names):
            a=float(mat[:,k].mean()); lo,hi=np.quantile(means[:,k],[.025,.975]); w=sum(np.argmin(mat,axis=1)==k)
            rows.append({'metric':metric,'model':m,'mean':a,'ci_95_lower':float(lo),'ci_95_upper':float(hi),'condition_winner_count':int(w)})
        # Pair differences of best fuzzy (ANFIS joint, FCM force) against simple options
        fuzzy='ANFIS-TS' if metric=='nrmse_pair' else 'FCM-TS'
        for candidate in MODELS:
            k,j=names.index(fuzzy),names.index(candidate)
            d=means[:,k]-means[:,j];lo,hi=np.quantile(d,[.025,.975])
            rows.append({'metric':metric,'model':f'Paired difference: {fuzzy} - {candidate}',
               'mean':float(mat[:,k].mean()-mat[:,j].mean()),'ci_95_lower':float(lo),'ci_95_upper':float(hi),'condition_winner_count':np.nan})
    boot=pd.DataFrame(rows);boot.to_csv(ROOT/'R7_bootstrap_7_models.csv',index=False)
    pd.DataFrame(provenance).to_csv(ROOT/'R7_input_sha256.csv',index=False)
    assert pd.read_csv(ROOT/'R7_simple_baseline_per_split.csv').shape[0]==270
    meta={'n_conditions':9,'n_seeds':10,'fuzzy_methods_frozen':4,'simple_models':list(MODELS),
      'fixed_seed':20261008,'bootstrap_iterations':B,'loop_samples':571,
      'training_n':400,'test_n':171,'eval_metric_normalization':'full loop span only for scoring, never model fitting',
      'endpoint_policy':'clip test phase to min and max TRAIN phase; constant endpoint value',
      'training_target_usage':'only train displacement and train force',
      'seed_protocol':'numpy.default_rng(seed).permutation(571), seed=0..9, first 400 train',
      'fuzzy_source':'frozen 360 archived R3 rows, not retrained',
      'no_prospective_use_claim':True}
    (ROOT/'R7_protocol.json').write_text(json.dumps(meta,indent=2))
    print('QA: INPUT NINE / ALL MODELS x TEN SEEDS',len(df),len(combined),'PASS')
    print('R7_AGGREGATE_BEGIN');print(agg.to_string(index=False,float_format=lambda n:f'{n:.8f}'))
    print('R7_PAIRED_CI_BEGIN');print(boot[boot.model.str.contains('Paired difference')].to_string(index=False,float_format=lambda n:f'{n:.8f}'))
    print('Boundary heldout test counts min/max',int(df['n_endpoint_heldout'].min()),int(df['n_endpoint_heldout'].max()))
if __name__=='__main__': main()

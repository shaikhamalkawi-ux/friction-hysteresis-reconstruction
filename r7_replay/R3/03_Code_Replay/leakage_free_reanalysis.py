#!/usr/bin/env python3
"""56 H2 R2 methodological sensitivity; original PoliTO numeric hyst arrays.

No target-based normalization or target-derived u is used to *fit* or evaluate
heldout force on these protocols. Native sample order is retained; interpolation
of the full target curve before splitting is forbidden.

The split is an interpolative point holdout within a known loop, NOT a new-cycle
or new-apparatus generalization estimate.
"""
import argparse,hashlib,sys,time,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.optimize import minimize
from importlib.machinery import SourceFileLoader
BASE=SourceFileLoader('baseline_fuzzy',str(Path(__file__).resolve().parent/'upstream_R1_Pilot01_pipeline.py')).load_module()
PREMISES=BASE.PREMISES
# Explicitly locked nine end-of-test files. Mid-cycle files may coexist in inputs
# and MUST NOT expand the confirmatory 3×3 cohort at replay time.
END_OF_TEST_FILES = (
    'AF10_1522500cycle.mat', 'AF19_997500cycle.mat',
    'AF24_1575000cycle.mat', 'AF28_1575000cycle.mat',
    'AF29_1575000cycle.mat', 'AF36_1575000cycle.mat',
    'AF38_1575000cycle.mat', 'AF39_1575000cycle.mat',
    'AF40_1575000cycle.mat',
)

def feature_u(x,mode):
    if mode=='index':
        return np.linspace(0.,1.,len(x))
    if mode=='x_displacement':
        d=np.abs(np.diff(x))
        if np.sum(d)<1e-20: raise ValueError('no displacement variation')
        return np.r_[0.,np.cumsum(d)]/np.sum(d)
    raise ValueError(mode)

def normalize_train(y,train):
    lo=float(np.min(y[train])); span=float(np.ptp(y[train]))
    if span<=0: raise ValueError('zero train dynamic range')
    return (y-lo)/span,lo,span

def fitted_predict(u,x,f,tr,model,mode,maxiter=46):
    # Fit in normalized ORIGINAL ordered samples, with TRAIN-ONLY scalers.
    xn,x0,xsp=normalize_train(x,tr)
    fn,f0,fsp=normalize_train(f,tr)
    if model!='ANFIS-TS':
        c,s=PREMISES[model]
        ax,bx=BASE.fit_ts(u,xn,c,s,tr)
        af,bf=BASE.fit_ts(u,fn,c,s,tr)
        nit=0
    else:
        c0,s0=PREMISES[model]; nr=len(c0)
        p0=np.r_[c0,np.log(s0)]
        bounds=[(0.,1.)]*nr+[(np.log(0.005),np.log(0.2))]*nr
        def unpack(p):
            c=p[:nr]; s=np.exp(p[nr:]); z=np.argsort(c)
            return c[z],s[z]
        def loss(p):
            c,s=unpack(p)
            ax,bx=BASE.fit_ts(u,xn,c,s,tr,ridge=1e-10)
            af,bf=BASE.fit_ts(u,fn,c,s,tr,ridge=1e-10)
            px=BASE.predict_ts(u[tr],c,s,ax,bx)
            pf=BASE.predict_ts(u[tr],c,s,af,bf)
            if mode=='x_displacement':
                return np.mean((pf-fn[tr])**2) # only force is target, x defines u
            return np.mean(0.5*((px-xn[tr])**2+(pf-fn[tr])**2))
        o=minimize(loss,p0,method='L-BFGS-B',bounds=bounds,
                   options={'maxiter':maxiter,'ftol':1e-12,'maxls':30})
        c,s=unpack(o.x)
        ax,bx=BASE.fit_ts(u,xn,c,s,tr,ridge=1e-10)
        af,bf=BASE.fit_ts(u,fn,c,s,tr,ridge=1e-10)
        nit=int(o.nit)
    px=x0+xsp*BASE.predict_ts(u,c,s,ax,bx)
    pf=f0+fsp*BASE.predict_ts(u,c,s,af,bf)
    return px,pf,nit

def metrics(x,f,px,pf,test,mode):
    ex=(px[test]-x[test])/np.ptp(x)
    ef=(pf[test]-f[test])/np.ptp(f)
    combined=float(np.sqrt(np.mean(0.5*(ex**2+ef**2))))
    force=float(np.sqrt(np.mean(ef**2)))
    # force-only metric valid for both protocols. Combined metric only for index.
    return {'nrmse_pair':combined if mode=='index' else np.nan,
            'nrmse_force':force,
            'nmae_force':float(np.mean(np.abs(ef))),
            'r2_force':float(1-np.sum((pf[test]-f[test])**2)/np.sum((f[test]-np.mean(f[test]))**2))}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--seeds',type=int,default=1)
    ap.add_argument('--methods',default='index,x_displacement')
    args=ap.parse_args()
    root=Path(__file__).parent; data=root/'inputs'; modes=args.methods.split(',')
    allrows=[]; hashrows=[]
    tic=time.time()
    for path in (data / name for name in END_OF_TEST_FILES):
        exp=path.name.split('_')[0]
        dat=path.read_bytes(); sha=hashlib.sha256(dat).hexdigest()
        h=loadmat(path,squeeze_me=True)['hyst']
        x=np.asarray(h[:,0],float); f=np.asarray(h[:,1],float)
        assert h.shape==(571,2) and np.isfinite(h).all()
        hashrows.append({'case':exp,'file':path.name,'sha256':sha,'N':len(x),'x_range':float(np.ptp(x)),'force_range':float(np.ptp(f))})
        for mode in modes:
            u=feature_u(x,mode)
            for seed in range(args.seeds):
                rng=np.random.default_rng(seed)
                p=rng.permutation(len(u)); ntr=int(round(0.7*len(u)))
                tr,te=p[:ntr],p[ntr:]
                for model in ('FCM-TS','SC9-TS','SC16-TS','ANFIS-TS'):
                    px,pf,nit=fitted_predict(u,x,f,tr,model,mode)
                    met=metrics(x,f,px,pf,te,mode)
                    allrows.append({'experiment':exp,'method':mode,'seed':seed,
                                    'model':model,'n_train':len(tr),'n_test':len(te),'optimizer_iters':nit,**met})
                    print(f'{exp} {mode:14s} seed={seed} {model:8s} force_NRMSE={met["nrmse_force"]:.6f} pair_NRMSE={met["nrmse_pair"]:.6f}',flush=True)
    df=pd.DataFrame(allrows); df.to_csv(root/'R2_target_blind_holdout_per_seed.csv',index=False)
    agg=df.groupby(['method','model'],as_index=False).agg(N=('nrmse_force','size'),mean_force_nrmse=('nrmse_force','mean'),median_force_nrmse=('nrmse_force','median'),mean_pair_nrmse=('nrmse_pair','mean'),mean_r2_force=('r2_force','mean'))
    agg.to_csv(root/'R2_target_blind_holdout_aggregate.csv',index=False)
    pd.DataFrame(hashrows).to_csv(root/'R2_raw_loop_provenance.csv',index=False)
    print('AGGREGATE\n',agg.to_string(index=False),flush=True)
    print('Runtime_s:',time.time()-tic,flush=True)
if __name__=='__main__':main()

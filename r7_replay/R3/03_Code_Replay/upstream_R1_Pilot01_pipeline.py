#!/usr/bin/env python3
"""
56 H2 R1 Pilot 01 — external friction-hysteresis replication.

Inputs:
    inputs/AF19_900cycle.mat
    inputs/AF19_504000cycle.mat
    inputs/AF19_997500cycle.mat
    inputs/AF10_900cycle.mat
    inputs/AF10_735000cycle.mat
    inputs/AF10_1522500cycle.mat
    inputs/AF36_900cycle.mat
    inputs/AF36_787500cycle.mat
    inputs/AF36_1575000cycle.mat

The external dataset is CC BY 4.0:
Fantetti et al., Mendeley Data, DOI 10.17632/gy587m7gx7.1.

This script keeps the H2 architecture counts fixed and does not treat the
metal-metal friction experiment as pneumatic-cylinder evidence.
"""
from pathlib import Path
import argparse, math
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.optimize import minimize

PREMISES = {
    "FCM-TS": (
        np.array([0.0736,0.2422,0.4139,0.5861,0.7578,0.9264]),
        np.array([0.0944]*6),
    ),
    "SC9-TS": (
        np.array([0.1018,0.2588,0.3422,0.4274,0.5125,0.5977,0.7145,0.8314,0.9399]),
        np.array([0.0354]*9),
    ),
    "SC16-TS": (
        np.array([0.0501,0.1285,0.2137,0.2571,0.3005,0.3673,0.4341,0.5008,
                  0.5693,0.6244,0.6795,0.7346,0.7913,0.8531,0.9132,0.9699]),
        np.array([0.0198]*16),
    ),
    "ANFIS-TS": (
        np.array([0.0849,0.2532,0.4428,0.4823,0.4902,0.6427,0.7536,0.7572,1.0000]),
        np.array([0.0368,0.0195,0.0335,0.0358,0.0311,0.0293,0.0384,0.0368,0.0345]),
    ),
}

CASES = [
    ("AF19","early",900),
    ("AF19","mid",504000),
    ("AF19","late",997500),
    ("AF10","early",900),
    ("AF10","mid",735000),
    ("AF10","late",1522500),
    ("AF36","early",900),
    ("AF36","mid",787500),
    ("AF36","late",1575000),
]

def load_loop(path):
    d=loadmat(path,squeeze_me=True)
    h=np.asarray(d["hyst"],dtype=float)
    return h[:,0],h[:,1]

def path_resample(x,f,n=None):
    if n is None: n=len(x)
    rx=np.ptp(x); rf=np.ptp(f)
    xn=(x-x.min())/rx
    fn=(f-f.min())/rf
    ds=np.sqrt(np.diff(xn)**2+np.diff(fn)**2)
    s=np.r_[0.0,np.cumsum(ds)]
    u0=s/s[-1]
    keep=np.r_[True,np.diff(u0)>1e-12]
    ug=np.linspace(0.0,1.0,n)
    xr=np.interp(ug,u0[keep],x[keep])
    fr=np.interp(ug,u0[keep],f[keep])
    return ug,xr,fr

def weights(u,c,s):
    mu=np.exp(-0.5*((u[:,None]-c[None,:])/s[None,:])**2)
    return mu/np.maximum(mu.sum(axis=1,keepdims=True),1e-15)

def fit_ts(u,y,c,s,idx,ridge=0.0):
    us=u[idx]; ys=y[idx]
    w=weights(us,c,s)
    A=np.column_stack([z for i in range(len(c)) for z in (w[:,i]*us,w[:,i])])
    if ridge>0:
        coef=np.linalg.solve(A.T@A+ridge*np.eye(A.shape[1]),A.T@ys)
    else:
        coef=np.linalg.lstsq(A,ys,rcond=None)[0]
    return coef[0::2],coef[1::2]

def predict_ts(u,c,s,a,b):
    w=weights(u,c,s)
    return (w*(a[None,:]*u[:,None]+b[None,:])).sum(axis=1)

def fixed_fit_predict(u,x,f,train_idx,c,s):
    xn=(x-x.min())/np.ptp(x)
    fn=(f-f.min())/np.ptp(f)
    ax,bx=fit_ts(u,xn,c,s,train_idx)
    af,bf=fit_ts(u,fn,c,s,train_idx)
    xhn=predict_ts(u,c,s,ax,bx)
    fhn=predict_ts(u,c,s,af,bf)
    return x.min()+np.ptp(x)*xhn, f.min()+np.ptp(f)*fhn

def adaptive_anfis(u,x,f,train_idx,maxiter=46):
    xn=(x-x.min())/np.ptp(x)
    fn=(f-f.min())/np.ptp(f)
    c0,s0=PREMISES["ANFIS-TS"]
    nr=len(c0)
    p0=np.r_[c0,np.log(s0)]
    bounds=[(0.0,1.0)]*nr+[(math.log(0.005),math.log(0.2))]*nr

    def unpack(p):
        c=p[:nr].copy()
        s=np.exp(p[nr:])
        o=np.argsort(c)
        return c[o],s[o]

    def objective(p):
        c,s=unpack(p)
        ax,bx=fit_ts(u,xn,c,s,train_idx,ridge=1e-10)
        af,bf=fit_ts(u,fn,c,s,train_idx,ridge=1e-10)
        px=predict_ts(u[train_idx],c,s,ax,bx)
        pf=predict_ts(u[train_idx],c,s,af,bf)
        return np.mean(0.5*((px-xn[train_idx])**2+(pf-fn[train_idx])**2))

    opt=minimize(objective,p0,method="L-BFGS-B",bounds=bounds,
                 options={"maxiter":maxiter,"ftol":1e-12,"maxls":30})
    c,s=unpack(opt.x)
    ax,bx=fit_ts(u,xn,c,s,train_idx,ridge=1e-10)
    af,bf=fit_ts(u,fn,c,s,train_idx,ridge=1e-10)
    px=predict_ts(u,c,s,ax,bx)
    pf=predict_ts(u,c,s,af,bf)
    xh=x.min()+np.ptp(x)*px
    fh=f.min()+np.ptp(f)*pf
    return xh,fh,c,s,opt.nit

def point_metrics(x,f,xh,fh,idx):
    rx=np.ptp(x); rf=np.ptp(f)
    ex=(x[idx]-xh[idx])/rx
    ef=(f[idx]-fh[idx])/rf
    ep=np.sqrt(0.5*(ex**2+ef**2))
    nrmse=float(np.sqrt(np.mean(ep**2)))
    nmae=float(np.mean(0.5*(np.abs(ex)+np.abs(ef))))
    nmax=float(ep.max())
    ss=np.sum((f[idx]-f[idx].mean())**2)
    r2=float(1-np.sum((f[idx]-fh[idx])**2)/ss) if ss>0 else np.nan
    return nrmse,nmae,nmax,r2

def loop_area(x,f):
    return abs(np.trapezoid(f,x))

def physical_metrics(x,f,xh,fh):
    nrmse=point_metrics(x,f,xh,fh,np.arange(len(x)))[0]
    a0=loop_area(x,f); a1=loop_area(xh,fh)
    area_err=100*abs(a1-a0)/max(a0,1e-15)
    closure=100*np.sqrt(0.5*(((xh[-1]-xh[0])/np.ptp(x))**2+
                            ((fh[-1]-fh[0])/np.ptp(f))**2))
    return nrmse,area_err,closure

def reversal_split(x,radius=12):
    test=set()
    for center in (int(np.argmin(x)),int(np.argmax(x))):
        test.update(range(max(0,center-radius),min(len(x),center+radius+1)))
    te=np.array(sorted(test),dtype=int)
    mask=np.ones(len(x),dtype=bool); mask[te]=False
    tr=np.where(mask)[0]
    return tr,te

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--data-dir",type=Path,default=Path("inputs"))
    ap.add_argument("--out-dir",type=Path,default=Path("rerun"))
    args=ap.parse_args()
    args.out_dir.mkdir(parents=True,exist_ok=True)

    rows=[]; revrows=[]
    for exp,stage,cycle in CASES:
        path=args.data_dir/f"{exp}_{cycle}cycle.mat"
        x0,f0=load_loop(path)
        u,x,f=path_resample(x0,f0,len(x0))
        n=len(u); ntr=int(round(0.7*n))

        for model in ("FCM-TS","SC9-TS","SC16-TS","ANFIS-TS"):
            # Full descriptive fit
            if model=="ANFIS-TS":
                xhf,fhf,_,_,full_nit=adaptive_anfis(u,x,f,np.arange(n))
            else:
                c,s=PREMISES[model]
                xhf,fhf=fixed_fit_predict(u,x,f,np.arange(n),c,s)
                full_nit=0
            full_nrmse,area_err,closure=physical_metrics(x,f,xhf,fhf)

            cv=[]
            for seed in range(10):
                rng=np.random.default_rng(seed)
                p=rng.permutation(n); tr=p[:ntr]; te=p[ntr:]
                if model=="ANFIS-TS":
                    xh,fh,_,_,nit=adaptive_anfis(u,x,f,tr)
                else:
                    c,s=PREMISES[model]
                    xh,fh=fixed_fit_predict(u,x,f,tr,c,s); nit=0
                cv.append((*point_metrics(x,f,xh,fh,te),nit))
            cv=np.asarray(cv,float)

            rows.append({
                "experiment":exp,"stage":stage,"cycle":cycle,
                "model":"ANFIS-TS (adaptive)" if model=="ANFIS-TS" else model,
                "full_nrmse":full_nrmse,
                "area_err_pct":area_err,
                "closure_err_pct":closure,
                "random_test_nrmse_mean":cv[:,0].mean(),
                "random_test_nrmse_sd":cv[:,0].std(ddof=1),
                "random_test_nmae_mean":cv[:,1].mean(),
                "random_test_nmax_mean":cv[:,2].mean(),
                "random_test_r2F_mean":cv[:,3].mean(),
                "opt_iter_mean":cv[:,4].mean(),
                "full_opt_iter":full_nit,
            })

            tr,te=reversal_split(x,12)
            if model=="ANFIS-TS":
                xh,fh,_,_,nit=adaptive_anfis(u,x,f,tr)
            else:
                c,s=PREMISES[model]
                xh,fh=fixed_fit_predict(u,x,f,tr,c,s); nit=0
            r=point_metrics(x,f,xh,fh,te)
            revrows.append({
                "experiment":exp,"stage":stage,"cycle":cycle,"model":model,
                "n_test":len(te),"rev_nrmse":r[0],"rev_nmae":r[1],
                "rev_nmax":r[2],"rev_r2F":r[3],"opt_iter":nit
            })

    pd.DataFrame(rows).to_csv(args.out_dir/"pilot_results.csv",index=False)
    pd.DataFrame(revrows).to_csv(args.out_dir/"reversal_holdout.csv",index=False)

if __name__=="__main__":
    main()

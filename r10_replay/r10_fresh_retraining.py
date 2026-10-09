#!/usr/bin/env python3
"""H2 R10: new, explicitly documented training of four TS configurations.
No archived source-paper center or consequent value is used in parameter fitting.
This is a retrospective description of a figure-digitized 600-point curve.
"""
from pathlib import Path
import argparse, hashlib, json, time
import numpy as np
import pandas as pd
from scipy.optimize import minimize

SEED_ORDER = range(10)
MODELS = ('FCM-TS', 'SC9-TS', 'SC16-TS', 'ANFIS-TS')
SOURCE = Path(__file__).parent/'source_R7'/'Presliding_Figure_Code_Package'/'data'/'Hysteresis_600pts.csv'
OUT = Path(__file__).parent/'output'
OPT_ITER=46

def fcm_1d(u, k=6, m=2, maxiter=500):
    # Original H2 employs FCM with m=2; this routine uses only u of fit subset.
    vals=np.asarray(u,float)
    c=np.quantile(vals, np.linspace(0.07, 0.93,k))
    for it in range(maxiter):
        dist=np.abs(vals[:,None]-c[None,:]); dist=np.maximum(dist,1e-12)
        memberships=dist**(-2/(m-1))
        memberships /= memberships.sum(axis=1,keepdims=True)
        mm=memberships**m
        cnew=(mm*vals[:,None]).sum(axis=0)/mm.sum(axis=0)
        if np.max(np.abs(cnew-c))<1e-12: c=cnew;break
        c=cnew
    ix=np.argsort(c); c=c[ix]; mm=mm[:,ix]
    sigma=np.sqrt((mm*(vals[:,None]-c[None,:])**2).sum(axis=0)/mm.sum(axis=0))
    if np.any(sigma<=1e-12):raise ValueError('degenerate FCM spread')
    return c,sigma,int(it+1)

def subtractive_1d(u,k):
    # Classical Chiu-like density/depletion with fixed architecture K.
    # r_a=1/K is fixed by the stated K model granularity; depletion 1.25*r_a
    # follows the source paper. Gaussian TS sigma=r_a/sqrt(8), reproducing
    # the exponent exp(-4 (u-c)^2/r_a^2) on normalized input.
    vals=np.sort(np.unique(np.asarray(u,float)))
    ra=1.0/k; rb=1.25*ra
    dd=vals[:,None]-vals[None,:]
    P=np.exp(-4*dd**2/ra**2).sum(axis=1)
    c=[]
    for j in range(k):
        p=P.copy()
        if c:
            for z in c: p[np.argmin(np.abs(vals-z))]=-np.inf
        idx=int(np.argmax(p)); power=float(P[idx])
        if not np.isfinite(power) or power<=0:raise ValueError(f'SC{k} insufficient positive potential at center {j+1}')
        center=float(vals[idx]); c.append(center)
        P=P-power*np.exp(-4*(vals-center)**2/rb**2)
    c=np.sort(np.array(c)); sig=np.full(k,ra/np.sqrt(8))
    return c,sig

def basis(u,c,s):
    u=np.asarray(u,float); c=np.asarray(c,float); s=np.asarray(s,float)
    z=(u[:,None]-c[None,:])/s[None,:]
    logmf=-0.5*z**2
    logmf -= logmf.max(axis=1,keepdims=True)
    mf=np.exp(logmf); w=mf/mf.sum(axis=1,keepdims=True)
    A=np.stack((w*u[:,None],w),axis=2).reshape(len(u),-1)
    return A,w

def get_linear(u,c,s,yy):
    A,w=basis(u,c,s)
    sol,_,rank,_=np.linalg.lstsq(A,yy,rcond=None)
    return sol,A,w,int(rank)

def adaptive(u,yy,c0,s0,budget=OPT_ITER):
    k=len(c0); init=np.r_[c0,np.log(s0)]
    bounds=[(0,1)]*k+[(np.log(.005),np.log(.2))]*k
    def fun_and_grad(p):
        c=p[:k]; s=np.exp(p[k:])
        beta,A,w,rank=get_linear(u,c,s,yy)
        pred=A@beta; res=pred-yy
        obj=float(np.mean(res*res))
        if not np.isfinite(obj):return 1e8,np.ones(len(p))*1e6
        b=beta.reshape(k,2,2) # each k antecedent: u and intercept; channels x,F
        local=u[:,None,None]*b[:,0,:][None,:,:]+b[:,1,:][None,:,:]
        dw=w[:,:,None]*(local-pred[:,None,:])
        grads_c=(dw*((u[:,None,None]-c[None,:,None])/(s[None,:,None]**2)))
        grads_logs=dw*((u[:,None,None]-c[None,:,None])**2/(s[None,:,None]**2))
        grad=np.r_[2*np.mean(np.sum(res[:,None,:]*grads_c,axis=2),axis=0),
                   2*np.mean(np.sum(res[:,None,:]*grads_logs,axis=2),axis=0)] / yy.shape[1]
        return obj,grad
    optim=minimize(fun_and_grad,init,method='L-BFGS-B',jac=True,bounds=bounds,
                   options={'maxiter':int(budget),'ftol':1e-12,'maxls':30,'gtol':1e-8})
    return optim.x[:k],np.exp(optim.x[k:]),optim

def fitted(u,x,f,train,model,budget=OPT_ITER):
    tr=np.array(train,dtype=int)
    if np.ptp(x[tr])==0 or np.ptp(f[tr])==0:raise ValueError('zero training span')
    low=np.array([x[tr].min(),f[tr].min()]);scale=np.array([np.ptp(x[tr]),np.ptp(f[tr])]);
    yy=(np.column_stack((x[tr],f[tr]))-low)/scale
    if model=='FCM-TS':
        c,s,nit=fcm_1d(u[tr]);info={'opt_iters':0,'opt_converged':True,'fcm_iters':nit}
    else:
        k=9 if model in ('SC9-TS','ANFIS-TS') else 16
        c,s=subtractive_1d(u[tr],k)
        info={'opt_iters':0,'opt_converged':True,'fcm_iters':0}
    if model=='ANFIS-TS':
        c,s,optim=adaptive(u[tr],yy,c,s,budget)
        info={'opt_iters':int(optim.nit),'opt_converged':bool(optim.success),
              'termination':str(optim.message),'fcm_iters':0}
    beta,A,w,rank=get_linear(u[tr],c,s,yy)
    cond=float(np.linalg.cond(A))
    B,_=basis(u,c,s)
    pred=(B@beta)*scale+low
    assert np.isfinite(pred).all() and rank==A.shape[1],f'Rank deficiency for {model}, rank {rank}/{A.shape[1]}'
    info.update({'design_rank':rank,'n_columns':A.shape[1],'condition':cond})
    return pred,c,s,beta,info

def scores(x,f,xy,indices):
    p=np.asarray(xy,float);t=np.asarray(indices,dtype=int)
    ex=(p[t,0]-x[t])/np.ptp(x);ef=(p[t,1]-f[t])/np.ptp(f)
    e2=(ex*ex+ef*ef)/2
    return {'nrmse_joint':float(np.sqrt(np.mean(e2))),
            'nrmse_force':float(np.sqrt(np.mean(ef**2))),
            'nmae_joint':float(np.mean((np.abs(ex)+np.abs(ef))/2)),
            'nmae_force':float(np.mean(np.abs(ef))),
            'nmax_joint':float(np.max(np.sqrt(e2))),
            'r2_force':float(1-np.dot(p[t,1]-f[t],p[t,1]-f[t])/np.sum((f[t]-f[t].mean())**2))}

def geometry(x,f,xy):
    xa,fa=xy[:,0],xy[:,1]
    a_ref=float(np.trapezoid(f,x));a_pred=float(np.trapezoid(fa,xa))
    endpoint=xy[-1]-xy[0]
    observed=np.array([x[-1]-x[0],f[-1]-f[0]])
    metric=lambda d:float(100*np.sqrt(0.5*np.sum((d/np.array([np.ptp(x),np.ptp(f)]))**2)))
    return {'signed_area_measured_Nmm':a_ref,'signed_area_pred_Nmm':a_pred,
            'loop_area_error_pct':float(100*abs(a_pred-a_ref)/abs(a_ref)),
            'closure_error_pct':metric(endpoint-observed),
            'predicted_endpoint_chord_pct':metric(endpoint),
            'observed_endpoint_chord_pct':metric(observed)}

def run(do_splits=True):
    OUT.mkdir(exist_ok=True)
    raw=SOURCE.read_bytes();df=pd.read_csv(SOURCE)
    u=df.u.to_numpy(); x=df.position_mm.to_numpy();f=df.force_N.to_numpy()
    assert len(df)==600 and np.all(np.isfinite(df.to_numpy()))
    assert np.allclose(u,np.linspace(0,1,600),rtol=0,atol=6e-14)
    model_rows=[];split_rows=[];params=[];pred_rows=[]
    for model in MODELS:
        tick=time.time()
        preds,c,s,beta,info=fitted(u,x,f,np.arange(600),model)
        point=scores(x,f,preds,np.arange(600));geo=geometry(x,f,preds)
        model_rows.append({'model':model,'fit':'fresh u-only premise construction + fresh consequent fit',
                           **point,**geo,**info})
        for j in range(len(c)):
            params.append({'model':model,'rule':j+1,'center':c[j],'sigma':s[j],
                           'coefficient_u_position_norm':beta.reshape(len(c),2,2)[j,0,0],
                           'coefficient_intercept_position_norm':beta.reshape(len(c),2,2)[j,1,0],
                           'coefficient_u_force_norm':beta.reshape(len(c),2,2)[j,0,1],
                           'coefficient_intercept_force_norm':beta.reshape(len(c),2,2)[j,1,1]})
        for i in range(len(u)):
            pred_rows.append({'model':model,'index':i,'u':u[i],'x_observed_mm':x[i],
                              'force_observed_N':f[i],'x_fitted_mm':preds[i,0],
                              'force_fitted_N':preds[i,1]})
        print('FULLFIT',model,'joint %.8f'%point['nrmse_joint'],'area %.5f%%'%geo['loop_area_error_pct'],
              'closure %.5f%%'%geo['closure_error_pct'],'cond %.3g'%info['condition'],
              'opt',info['opt_iters'],'time %.2fs'%(time.time()-tick),flush=True)
    if do_splits:
        for seed in SEED_ORDER:
            p=np.random.default_rng(seed).permutation(len(u));nt=round(0.7*len(u));tr=p[:nt];te=p[nt:]
            for model in MODELS:
                tick=time.time();xy,c,s,beta,info=fitted(u,x,f,tr,model)
                for kind,ix in [('train',tr),('test',te)]:
                    split_rows.append({'model':model,'seed':seed,'subset':kind,'n':len(ix),
                                       **scores(x,f,xy,ix),**info})
                print('SPLIT',seed,model,'test %.8f'%split_rows[-1]['nrmse_joint'],
                      'opt',info['opt_iters'],info['opt_converged'],'%.2fs'%(time.time()-tick),flush=True)
    pd.DataFrame(model_rows).to_csv(OUT/'R10_new_primary_fullfit.csv',index=False)
    pd.DataFrame(params).to_csv(OUT/'R10_fresh_model_parameters.csv',index=False)
    pd.DataFrame(pred_rows).to_csv(OUT/'R10_new_primary_predicted_600.csv',index=False)
    if split_rows:
        pd.DataFrame(split_rows).to_csv(OUT/'R10_new_primary_point_splits.csv',index=False)
        summ=(pd.DataFrame(split_rows).query("subset=='test'").groupby('model',as_index=False)
              .agg(n_splits=('seed','nunique'),nrmse_joint_mean=('nrmse_joint','mean'),
                   nrmse_joint_sd=('nrmse_joint','std'),nrmse_force_mean=('nrmse_force','mean'),
                   nmae_joint_mean=('nmae_joint','mean'),nmax_joint_mean=('nmax_joint','mean'),
                   condition_median=('condition','median'), optimizer_iterations_mean=('opt_iters','mean')))
        summ.to_csv(OUT/'R10_new_primary_point_splits_summary.csv',index=False)
    protocol={
      'status':'new retrospective figure-derived 600-point recalculation, NOT recovery of historical optimizer',
      'source':'Hysteresis_600pts.csv (figure digitized Nouri via Recchia)',
      'source_sha256':hashlib.sha256(raw).hexdigest(),
      'predictor':'u recorded normalized sample-path index on figure-digitized/resampled 600-point curve; original digitization and resampling used plotted force; not prospective response-blind input for source figure',
      'fit':'FCM six centers m=2 u-only; SC9/SC16 Chiu-like density/depletion ra=1/K, rb=1.25ra and sigma=ra/sqrt(8) u-only; ANFIS from freshly computed SC9 init, L-BFGS-B joint normalized loss maxiter=46',
      'rationale':'Four architecture rule counts and FCM exponent inherited from H2; SC radii/widths tied to count with mathematical kernel matching (no response tuning); 46 maxiter inherited from archived manuscript',
      'least_squares':'numpy SVD lstsq no ridge; reject rank deficiency rather than retune',
      'validation':'fresh descriptive fullfit plus ten random 70/30 point holdouts seeds 0..9; training-only scales for prediction; joint train fit equal normalized channels (same source metric); no independent predictive validity asserted',
      'metrics':'joint NRMSE = sqrt(mean((ex^2+ef^2)/2)); NMAE = mean((abs(ex)+abs(ef))/2); force NRMSE; area = |100*(Ahat-A)/A|, A signed trapezoidal integral of F against x (N-mm); closure=100*RMS normalized difference in predicted and observed end-to-start chords, observed chord zero; max uses normalized pointwise e',
      'original_source_table7':'unresolved archived source values are excluded from conclusions, preserved separately verbatim',
      'versions':{'numpy':np.__version__,'pandas':pd.__version__},
      'fit_parameters':{'FCM_m':2,'gaussian_width':'weighted FCM standard deviation or SC kernel matching','ANFIS_maxiter':46,'L_BFGS_B_bounds_centers':[0,1],'L_BFGS_B_bounds_sigma':[0.005,0.2]}
    }
    (OUT/'R10_primary_protocol.json').write_text(json.dumps(protocol,indent=2,ensure_ascii=False)+'\n')
    print('FINISHED',len(model_rows),'full models',len(split_rows),'split rows',flush=True)

if __name__=='__main__':
    pa=argparse.ArgumentParser();pa.add_argument('--no-splits',action='store_true');a=pa.parse_args();run(not a.no_splits)

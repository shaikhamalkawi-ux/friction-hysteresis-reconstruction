#!/usr/bin/env python3
"""H2 R11 supplementary analysis. Run from an extracted R11 work/replay archive.

Extended nine-condition, 10-seed within-loop point completion with train-only
reversal-neighborhood definition and completed-curve signed-loop-area proxy;
time-order examination of five later cycles after one early-cycle fit in each
of THREE physical conditions. No held-out force is used in model fitting.

This does NOT assess independent pneumatic-seal validation or true deployment
where phase/cycle length may be unavailable. The selected 9 loops and older
R7/R10 numerical results remain frozen.
"""
from __future__ import annotations
import argparse, hashlib, json, math, sys
from pathlib import Path
from importlib.machinery import SourceFileLoader
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parent.parent
R7 = ROOT / 'r7_replay'
SCI = ROOT / 'science'
OUT = SCI / 'results'
OUT.mkdir(parents=True,exist_ok=True)
R = SourceFileLoader('R11_r3',str(R7/'R3/03_Code_Replay/leakage_free_reanalysis.py')).load_module()
B = SourceFileLoader('R11_baseline',str(R7/'R7_simple_baselines.py')).load_module()
MODELS=('Linear (clamped)','PCHIP (clamped)','Cubic spline (clamped)','FCM-TS','SC9-TS','SC16-TS','ANFIS-TS')


def sha(p:Path)->str: return hashlib.sha256(p.read_bytes()).hexdigest()

def polygon_signed_area(x:np.ndarray,y:np.ndarray)->float:
    """Closed polygon 1/2 sum(x_i*y_{i+1}-x_{i+1}*y_i), signed; x,y calibrated native axes."""
    return float(0.5*np.sum(x*np.roll(y,-1)-np.roll(x,-1)*y))

def paths_model_fit(u,x,f,tr,model):
    if model in B.MODELS:
        te=np.arange(len(u));px=B.predict(x,u,tr,te,model);pf=B.predict(f,u,tr,te,model)
        return np.asarray(px),np.asarray(pf),0
    return R.fitted_predict(u,x,f,tr,model,'index',maxiter=46)

def run_within(nseeds=10):
    detail=[];input_records=[]
    archived=pd.read_csv(R7/'R7_all_7_models_per_split.csv')
    archived=archived[archived.method=='index'].copy()
    assert archived.shape[0]==630
    for exp,cyc in sorted(B.END.items()):
        path=B.DATA/f'{exp}_{cyc}cycle.mat'
        a=np.asarray(loadmat(path,squeeze_me=True)['hyst'],float);assert a.shape==(571,2)
        input_records.append({'task':'within','condition':exp,'snapshot':cyc,'filename':path.name,'sha256':sha(path),'rows':len(a)})
        x,f=a.T; u=np.linspace(0,1,len(x));A=polygon_signed_area(x,f);rx=np.ptp(x);rf=np.ptp(f)
        assert abs(A)>1e-8,(exp,A)
        for seed in range(nseeds):
            rnd=np.random.default_rng(seed).permutation(len(u));tr,te=rnd[:400],rnd[400:]
            assert len(te)==171
            # Define reversal regions using TRAIN displacement only (not heldout force/displacement).
            extremes=[int(tr[np.argmin(x[tr])]),int(tr[np.argmax(x[tr])])]
            mask=np.zeros(len(u),bool)
            for k in extremes: mask[max(0,k-20):min(len(u),k+21)]=True
            test_near=te[mask[te]]
            assert len(test_near)>=8, (exp,seed,test_near)
            for model in MODELS:
                px,pf,niter=paths_model_fit(u,x,f,tr,model)
                m=B.met(x,f,px[te],pf[te],te)
                e_near=(pf[test_near]-f[test_near])/rf
                nrmse_near=float(np.sqrt(np.mean(e_near**2)))
                completed_x=x.copy();completed_f=f.copy()
                completed_x[te]=px[te];completed_f[te]=pf[te]
                Acomp=polygon_signed_area(completed_x,completed_f)
                rel_area=100.*abs(Acomp-A)/abs(A)
                rev_coverage=float(len(test_near)/len(te))
                detail.append({'condition':exp,'seed':seed,'model':model,'n_training':400,'n_heldout':171,
                    'n_reversal_test':len(test_near),'reversal_fraction_heldout':rev_coverage,
                    'reversal_tip_train_min':extremes[0],'reversal_tip_train_max':extremes[1],
                    'nrmse_pair':m['nrmse_pair'],'nrmse_force':m['nrmse_force'],
                    'reversal_force_nrmse':nrmse_near,
                    'true_signed_area_native_force_displacement':A,
                    'completed_signed_area_native_force_displacement':Acomp,
                    'completed_area_error_pct':rel_area,'opt_iters':niter})
        print('WITHIN COMPLETE',exp,flush=True)
    d=pd.DataFrame(detail)
    assert len(d)==9*nseeds*7,len(d)
    p=OUT/'R11_within_reversal_area_per_split.csv'; d.to_csv(p,index=False,float_format='%.14g')
    # Compare unchanged registered R7 seven-model point errors against R11 recomputation.
    if nseeds==10:
        m=d.merge(archived,on=['condition' if 'condition' in archived.columns else 'experiment','seed','model']) if False else None
        archive=archived.rename(columns={'experiment':'condition'})
        merged=d.merge(archive,on=['condition','seed','model'],suffixes=('_r11','_r7'),validate='one_to_one')
        assert len(merged)==630
        deltas={col:float(np.max(np.abs(merged[col+'_r11']-merged[col+'_r7']))) for col in ['nrmse_pair','nrmse_force']}
        assert deltas['nrmse_pair']<2e-8 and deltas['nrmse_force']<2e-8,deltas
        (OUT/'R11_R7_score_parity.json').write_text(json.dumps({'repeated_archived_scores':len(merged),'maximum_absolute_error':deltas,'status':'PASS'},indent=2)+'\n')
    per_case=d.groupby(['condition','model'],as_index=False)[['nrmse_pair','nrmse_force','reversal_force_nrmse','completed_area_error_pct']].mean()
    per_case.to_csv(OUT/'R11_within_reversal_area_per_condition.csv',index=False,float_format='%.12g')
    agg=per_case.groupby('model',as_index=False)[['nrmse_pair','nrmse_force','reversal_force_nrmse','completed_area_error_pct']].agg(['mean','median'])
    agg.columns=['model']+[f'{key}_{stat}' for key,stat in agg.columns.tolist()[1:]] if False else ['model']+[str(k)+'_'+str(v) for k,v in agg.columns.tolist()[1:]]
    agg.to_csv(OUT/'R11_within_reversal_area_aggregate.csv',index=False,float_format='%.12g')
    # Paired nine-condition bootstrap under a fixed RNG. Do not misinterpret as 90 independent loops.
    condition_ids=sorted(per_case.condition.unique());rng=np.random.default_rng(20261009)
    idx=rng.integers(0,9,size=(50000,9))
    comparisons=[]
    for metric in ['reversal_force_nrmse','completed_area_error_pct','nrmse_force']:
        mat=np.array([[per_case[(per_case.condition==c)&(per_case.model==model)][metric].iloc[0] for model in MODELS] for c in condition_ids])
        assert mat.shape==(9,7)
        for a_name,b_name in [('ANFIS-TS','Linear (clamped)'),('SC16-TS','PCHIP (clamped)'),('FCM-TS','Linear (clamped)')]:
            ia,ib=MODELS.index(a_name),MODELS.index(b_name)
            differences=mat[:,ia]-mat[:,ib];means=differences[idx].mean(axis=1)
            lo,hi=np.quantile(means,[.025,.975]);comparisons.append({'metric':metric,'model_A':a_name,'model_B':b_name,'difference_A_minus_B':differences.mean(),'bootstrap_95_low':lo,'bootstrap_95_high':hi,'n_conditions':9})
    pd.DataFrame(comparisons).to_csv(OUT/'R11_paired_within_condition_bootstrap.csv',index=False,float_format='%.12g')
    pd.DataFrame(input_records).to_csv(OUT/'R11_within_source_SHA256.csv',index=False)
    print('WITHIN AGGREGATE\n',agg.to_string(index=False,float_format=lambda v:f'{v:.6g}'))
    return d,agg

TEMPORAL_CYCLES={
 'AF10':[900,304500,609000,892500,1207500,1522500],
 'AF19':[900,199500,399000,609000,787500,997500],
 'AF36':[900,315000,630000,945000,1260000,1575000],
}

def run_temporal():
    rows=[];sha_records=[];earlyfit_records=[]
    for exp,vals in TEMPORAL_CYCLES.items():
        D=SCI/'temporal_inputs'/exp
        values=[]
        for v in vals:
            path=D/f'{v}cycle.mat';assert path.exists(),path
            a=np.asarray(loadmat(path,squeeze_me=True)['hyst'],float);assert a.shape==(571,2) and np.isfinite(a).all()
            sha_records.append({'task':'temporal','condition':exp,'snapshot':v,'filename':f'{exp}/{path.name}','sha256':sha(path),'rows':571})
            values.append(a)
        early=values[0];u=np.linspace(0,1,571);x0,f0=early.T
        cache={}
        for model in (*MODELS[3:],'Unfitted early-cycle template'):
            if model=='Unfitted early-cycle template':xp,fp=x0.copy(),f0.copy();nit=0
            else: xp,fp,nit=paths_model_fit(u,x0,f0,np.arange(571),model)
            cache[model]=(xp,fp,nit)
            earlyfit_records.append({'condition':exp,'model':model,'optimizer_iters':nit,'in_sample_force_nrmse':float(np.sqrt(np.mean(((fp-f0)/np.ptp(f0))**2)))})
        for k in range(1,6):
            x,f=values[k].T
            true_area=polygon_signed_area(x,f)
            for model,(xp,fp,nit) in cache.items():
                ef=(fp-f)/np.ptp(f);ex=(xp-x)/np.ptp(x)
                force_rmse=float(np.sqrt(np.mean(ef**2)));pair_rmse=float(np.sqrt(np.mean((ef**2+ex**2)/2)))
                area_pred=polygon_signed_area(xp,fp)
                rows.append({'condition':exp,'test_cycle':vals[k],'fraction_of_final':vals[k]/vals[-1],
                    'train_cycle':900,'model':model,'n_test':571,'n_training':571,
                    'force_nrmse':force_rmse,'joint_nrmse':pair_rmse,
                    'signed_area_observed':true_area,'signed_area_predicted':area_pred,
                    'loop_area_error_pct':100*abs(area_pred-true_area)/abs(true_area),
                    'inference_phase':'fixed index i/570; measured cycle lengths 571, not inferred from outcomes',
                    'mean_force_shift':float(np.mean(f)-np.mean(f0))})
        print('TEMPORAL COMPLETE',exp,flush=True)
    frame=pd.DataFrame(rows)
    assert len(frame)==3*5*5
    frame.to_csv(OUT/'R11_temporal_extended_per_case.csv',index=False,float_format='%.14g')
    s=frame.groupby(['condition','model'],as_index=False)[['force_nrmse','joint_nrmse','loop_area_error_pct']].mean()
    s.to_csv(OUT/'R11_temporal_extended_condition_means.csv',index=False,float_format='%.12g')
    overall=s.groupby('model',as_index=False)[['force_nrmse','joint_nrmse','loop_area_error_pct']].mean()
    overall.to_csv(OUT/'R11_temporal_extended_aggregate.csv',index=False,float_format='%.12g')
    (OUT/'R11_temporal_fit_diagnostics.csv').write_text(pd.DataFrame(earlyfit_records).to_csv(index=False))
    pd.DataFrame(sha_records).to_csv(OUT/'R11_temporal_source_SHA256.csv',index=False)
    # Existing early-late mean should have parity with R3 archive up to numerical tolerance.
    last=frame.groupby('condition')['test_cycle'].max().to_dict()
    final=frame[frame.apply(lambda row:row.test_cycle==last[row.condition],axis=1)]
    previous=pd.read_csv(R7/'R3/02_Results_and_Methods/R3_cross_cycle_extended_summary.csv')
    previous=previous[previous.scenario=='early_to_late'].set_index('model')
    translate={'Unfitted early-cycle template':'unfitted_phase_template'}
    parity={}
    for model in final.model.unique():
        b=translate.get(model,model)
        newmean=final.loc[final.model==model,'force_nrmse'].mean()
        oldmean=float(previous.loc[b,'mean_force_nrmse'])
        parity[model]={'new':float(newmean),'archived':oldmean,'abs_diff':abs(newmean-oldmean)}
    (OUT/'R11_temporal_R3_early_late_parity.json').write_text(json.dumps(parity,indent=2)+'\n')
    assert max(v['abs_diff'] for v in parity.values())<5e-8,parity
    print('TEMPORAL AGGREGATE\n',overall.to_string(index=False,float_format=lambda v:f'{v:.7g}'))
    return frame,overall


def summarize_blocked():
    df=pd.read_csv(R7/'R7_blocked_gap_per_case.csv')
    assert len(df)==27*7
    out=df.groupby('model').agg(n=('nrmse_force','size'),mean_force=('nrmse_force','mean'),median_force=('nrmse_force','median'),
                                q90_force=('nrmse_force',lambda a:a.quantile(.9)),n_force_error_gt1=('nrmse_force',lambda a:int((a>1).sum())),
                                mean_joint=('nrmse_pair','mean'),median_joint=('nrmse_pair','median')).reset_index()
    out.to_csv(OUT/'R11_blocked_gap_robust_summary.csv',index=False,float_format='%.12g')
    print('BLOCKED GAP AGGREGATE\n',out.to_string(index=False,float_format=lambda v:f'{v:.6g}'))
    return out


def write_protocol():
    p={'protocol_version':'H2_R11_EVIDENCE_EXTENSION','within_loop':{'n_distinct_experimental_conditions':9,'n_seeds_each':10,'n_train':400,'n_test':171,'n_samples_per_loop':571,
         'reversal_neighborhood_halfwidth_indices':20,'tip_detection':'train-only observed x global minimum and global maximum index',
         'area_metric':'closed signed polygon, completed curve replaces heldout x and force with fitted predictions, train observations retained',
         'percentage_denominator':'absolute signed area of observed full loop, evaluation only',
         'report':'mean per condition first, aggregate across nine conditions, bootstrap at condition level, B=50000, seed=20261009'},
       'time_ordered':{'condition_names':list(TEMPORAL_CYCLES),'cycles':TEMPORAL_CYCLES,'n_distinct_conditions':3,'n_later_cycles_each':5,
         'training':'one early cycle per condition, all 571 points; no later fit/recalibration',
         'baseline':'unfitted phase-aligned early measurement used at all future target cycles','predictor':'known/common index i/570; no prediction of future phase alignment',
         'status':'exploratory repeated timepoints within three experimental conditions, not 15 independent experimental subjects'},
       'blocked_gap':'reanalysis of all 27 published R7 within-loop 171-contiguous-point tests, no winner-only selection',
       'no_new_pneumatic_seal_data':True,'old_results_unchanged':True,
       'not_physical_dissipation_certificate':True,
       'requires_feature_availability':'predefined matched sample index and complete length of 571 points'},
    (OUT/'R11_protocol_and_scope.json').write_text(json.dumps(p,indent=2)+'\n')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--tasks',default='within,temporal,blocked');ap.add_argument('--seeds',type=int,default=10)
    args=ap.parse_args();write_protocol()
    for t in args.tasks.split(','):
        if t=='within': run_within(args.seeds)
        elif t=='temporal':run_temporal()
        elif t=='blocked':summarize_blocked()
        else:raise ValueError(t)
    print('R11_REPLAY_PASS',args.tasks,flush=True)
if __name__=='__main__':main()

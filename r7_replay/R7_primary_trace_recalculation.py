#!/usr/bin/env python3
"""New primary-figure-digitization calibration, not a reconstruction of historical optimizers.
All results explicitly belong to one recovered 600-point plotted trajectory.
"""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from importlib.machinery import SourceFileLoader
ROOT=Path(__file__).resolve().parent
source=ROOT/'Presliding_Figure_Code_Package'/'data'/'Hysteresis_600pts.csv'
replay=SourceFileLoader('archived_figures',str(ROOT/'Presliding_Figure_Code_Package'/'python'/'presliding_figures.py')).load_module()
R=SourceFileLoader('r3',str(ROOT/'R3/03_Code_Replay/leakage_free_reanalysis.py')).load_module()
a=pd.read_csv(source)
u=a.u.to_numpy();x=a.position_mm.to_numpy();f=a.force_N.to_numpy();tr=np.arange(len(u))
assert len(u)==600 and np.allclose(u,np.linspace(0,1,600))
models=['FCM-TS','SC9-TS','SC16-TS','ANFIS-TS']
old_claim={'FCM-TS':0.02835,'SC9-TS':0.01782,'SC16-TS':0.00360,'ANFIS-TS':0.00397}
rows=[];line=[]
for model in models:
    xold,fold=replay.model_prediction(model,u,x,f)
    xnew,fnew,nit=R.fitted_predict(u,x,f,tr,model,'index',maxiter=46)
    for status,xp,fp in [('rounded-parameter replay',xold,fold),('new full-data fit',xnew,fnew)]:
        ex=(xp-x)/np.ptp(x);ef=(fp-f)/np.ptp(f)
        nr=float(np.sqrt(np.mean((ex*ex+ef*ef)/2)))
        nforce=float(np.sqrt(np.mean(ef**2)))
        true_area=float(np.trapezoid(f,x)); pred_area=float(np.trapezoid(fp,xp))
        area_err=float(100*abs(pred_area-true_area)/max(abs(true_area),1e-15))
        d0=np.array([x[-1]-x[0],f[-1]-f[0]])
        d1=np.array([xp[-1]-xp[0],fp[-1]-fp[0]])
        closure_vec=float(np.sqrt(.5*np.sum(((d1-d0)/np.array([np.ptp(x),np.ptp(f)]))**2)))
        rows.append({'model':model,'status':status,'nrmse_pair':nr,'nrmse_force':nforce,
                     'signed_area_pct_difference':area_err,'closure_vector_norm_discrepancy':closure_vec,'optimizer_iters':int(nit) if status.startswith('new') else np.nan,
                     'historical_table7_claimed_nrmse':old_claim[model],'delta_from_historical_NRMSE':nr-old_claim[model]})
    line.append({'model':model,'historical':old_claim[model],'archived_round_replay':rows[-2]['nrmse_pair'],'new_fitted':rows[-1]['nrmse_pair']})
pd.DataFrame(rows).to_csv(ROOT/'R7_primary_600point_fullfit_recalculated.csv',index=False)
pd.DataFrame(line).to_csv(ROOT/'R7_Table7_original_vs_verified.csv',index=False)
(ROOT/'R7_primary_600point_method.json').write_text(json.dumps({'source_csv':str(source.relative_to(ROOT)),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
 'limits':'600 points digitized from published plotted experiment, not raw sensors; new fit not historical recovered model',
 'new_fixed_premise_parameters':'from R3 known original source rounded premise centers/widths; consequent train all 600',
 'new_adaptive_method':'L-BFGS-B maxiter=46, archived ANFIS centers used as initial values',
 'metric':'RMS of equally weighted squared displacement and force normalized by their full 600-point channel ranges',
 'area':'signed trapezoidal integral of force over displacement with percentage difference; not independent energy measurement',
 'closure':'norm of difference between predicted and observed net endpoint displacement-force vector, normalized by full spans; NOT original historical closure metric'},indent=2))
print('PASS: 600 points, original four rounded-parameter replays, four independently fitted descriptions')
print(pd.DataFrame(rows).to_string(index=False,float_format=lambda v:f'{v:.7f}'))

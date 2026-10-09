#!/usr/bin/env python3
"""Fixed-premise TS Tikhonov sensitivity, on the 27 prespecified blocked-gap tests.
Sensitivity grid fixed before reading evaluation outcomes in script. Not used for
selecting or tuning one preferred penalty; includes all results.
"""
from pathlib import Path
from importlib.machinery import SourceFileLoader
import pandas as pd, numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parent
R=SourceFileLoader('r3',str(ROOT/'R3/03_Code_Replay/leakage_free_reanalysis.py')).load_module()
B=SourceFileLoader('b',str(ROOT/'R7_simple_baselines.py')).load_module()
rows=[];alphas=[0.0,1e-10,1e-8,1e-6,1e-4]
for exp,cycle in sorted(B.END.items()):
 a=loadmat(B.DATA/f'{exp}_{cycle}cycle.mat',squeeze_me=True)['hyst']; x,f=a.T;u=np.linspace(0,1,571)
 for st in [60,200,340]:
  te=np.arange(st,st+171); tr=np.r_[np.arange(0,st),np.arange(st+171,571)]
  xn,x0,xr=R.normalize_train(x,tr);fn,f0,fr=R.normalize_train(f,tr)
  for model in ['FCM-TS','SC9-TS','SC16-TS']:
   c,s=R.PREMISES[model]
   for alpha in alphas:
    ax,bx=R.BASE.fit_ts(u,xn,c,s,tr,ridge=alpha)
    af,bf=R.BASE.fit_ts(u,fn,c,s,tr,ridge=alpha)
    px=x0+xr*R.BASE.predict_ts(u[te],c,s,ax,bx)
    pf=f0+fr*R.BASE.predict_ts(u[te],c,s,af,bf)
    met=B.met(x,f,px,pf,te)
    rows.append({'experiment':exp,'start':st,'model':model,'ridge':alpha,**met})
res=pd.DataFrame(rows);res.to_csv(ROOT/'R7_fixed_TS_ridge_sensitivity_all.csv',index=False)
agg=res.groupby(['model','ridge'],as_index=False).agg(n=('nrmse_pair','size'),mean_pair=('nrmse_pair','mean'),median_pair=('nrmse_pair','median'),max_pair=('nrmse_pair','max'))
agg.to_csv(ROOT/'R7_fixed_TS_ridge_sensitivity_summary.csv',index=False)
print('PASS',len(res),'27 blocked gaps x 3 models x 5 frozen penalty settings')
print(agg.to_string(index=False,float_format=lambda x:f'{x:.6g}'))

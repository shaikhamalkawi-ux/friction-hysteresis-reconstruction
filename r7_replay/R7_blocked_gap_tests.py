#!/usr/bin/env python3
"""Prespecified internal 171-point contiguous gaps in 571-point PoliTO loops;
400 remaining points for training. Three gaps per nine conditions, fixed starts
60, 200, 340. Ends remain observed to test interpolation rather than extrapolation.
No heldout force or position enters fit and no reparameterization by target.
This evaluates geometry completion with known sample index, NOT cycle transfer.
"""
from pathlib import Path
import pandas as pd,numpy as np
from scipy.io import loadmat
from importlib.machinery import SourceFileLoader
ROOT=Path(__file__).resolve().parent
R=SourceFileLoader('r3',str(ROOT/'R3/03_Code_Replay/leakage_free_reanalysis.py')).load_module()
B=SourceFileLoader('base',str(ROOT/'R7_simple_baselines.py')).load_module()
all=[]
for exp,cycle in sorted(B.END.items()):
  a=np.asarray(loadmat(B.DATA/f'{exp}_{cycle}cycle.mat',squeeze_me=True)['hyst'],float);x,f=a.T
  assert a.shape==(571,2)
  u=np.linspace(0,1,len(x))
  for start in [60,200,340]:
    test=np.arange(start,start+171);train=np.r_[np.arange(0,start),np.arange(start+171,571)]
    assert len(train)==400 and len(test)==171 and train.min()==0 and train.max()==570
    for model in ['FCM-TS','SC9-TS','SC16-TS','ANFIS-TS',*B.MODELS]:
      if model in B.MODELS:
        px=B.predict(x,u,train,test,model);pf=B.predict(f,u,train,test,model);it=0
      else:
        xp,fp,it=R.fitted_predict(u,x,f,train,model,'index',maxiter=46)
        px,pf=xp[test],fp[test]
      y=B.met(x,f,px,pf,test)
      all.append({'experiment':exp,'start':start,'model':model,'n_train':400,'n_test':171,'opt_iters':it,**y})
    print('FINISHED',exp,start,flush=True)
result=pd.DataFrame(all); result.to_csv(ROOT/'R7_blocked_gap_per_case.csv',index=False)
summary=result.groupby('model')[['nrmse_pair','nrmse_force']].mean().sort_values('nrmse_pair')
summary.to_csv(ROOT/'R7_blocked_gap_summary.csv')
print('PASS:',len(result),'=9 conditions x 3 internal windows x 7 methods')
print(summary.to_string(float_format=lambda v:f'{v:.7f}'))

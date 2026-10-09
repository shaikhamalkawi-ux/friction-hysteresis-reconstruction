#!/usr/bin/env python3
"""Predeclared longer-budget sensitivity for nonconverged baseline; not primary retuning."""
from pathlib import Path
import numpy as np,pandas as pd
from r10_fresh_retraining import SOURCE, fitted, scores,geometry
D=pd.read_csv(SOURCE);u=D.u.to_numpy();x=D.position_mm.to_numpy();f=D.force_N.to_numpy()
results=[]
for seed in [-1]+list(range(10)):
    if seed<0: tr=np.arange(600);te=tr
    else:
        p=np.random.default_rng(seed).permutation(600);tr=p[:420];te=p[420:]
    for budget in (46,150):
        pred,c,s,coef,info=fitted(u,x,f,tr,'ANFIS-TS',budget)
        r={'seed':seed,'budget':budget,'n_test':len(te),**scores(x,f,pred,te),**info}
        if seed<0:r.update(geometry(x,f,pred))
        results.append(r)
        print('ANFIS budget sensitivity',seed,budget,'NRMSE joint',r['nrmse_joint'],
              'iterations',info['opt_iters'],'converged',info['opt_converged'],flush=True)
OUT=Path(__file__).parent/'output';pd.DataFrame(results).to_csv(OUT/'R10_ANFIS_46_vs_150_budget.csv',index=False)
print('DONE longer budget, primary remains 46 regardless of outcome')

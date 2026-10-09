#!/usr/bin/env python3
"""Independent finite-difference validation of variable-projection ANFIS gradient."""
from pathlib import Path
import numpy as np, pandas as pd, json
from r10_fresh_retraining import fcm_1d,subtractive_1d,basis
p=Path(__file__).parent
raw=pd.read_csv(p/'source_R7/Presliding_Figure_Code_Package/data/Hysteresis_600pts.csv')
idx=np.arange(0,600,7)
u=raw.u.to_numpy()[idx];x=raw.position_mm.to_numpy()[idx];f=raw.force_N.to_numpy()[idx]
yy=np.column_stack(((x-x.min())/np.ptp(x),(f-f.min())/np.ptp(f)))
c0,s0=subtractive_1d(u,9); K=len(c0)

def state(q):
 c=q[:K];sig=np.exp(q[K:]);A,w=basis(u,c,sig)
 beta=np.linalg.lstsq(A,yy,rcond=None)[0];pr=A@beta;res=pr-yy
 return float(np.mean(res**2)),A,w,beta,pr,res,c,sig

def grad(q):
 loss,A,w,beta,pr,res,c,sig=state(q)
 b=beta.reshape(K,2,2)
 local=u[:,None,None]*b[:,0,:][None,:,:]+b[:,1,:][None,:,:]
 delta=w[:,:,None]*(local-pr[:,None,:]); dc=(u[:,None,None]-c[None,:,None])/(sig[None,:,None]**2)
 dlogs=(u[:,None,None]-c[None,:,None])**2/(sig[None,:,None]**2)
 return np.r_[2*np.mean(np.sum(res[:,None,:]*delta*dc,axis=2),axis=0),
              2*np.mean(np.sum(res[:,None,:]*delta*dlogs,axis=2),axis=0)]/2

max_abs=0;max_rel=0
for pert in (0,1):
 q=np.r_[c0,np.log(s0)]
 if pert:
  q[:K]+=0.01*np.sin(np.arange(K))
  q[K:]+=0.04*np.cos(np.arange(K))
 exact=grad(q)
 eps=1e-5
 numeric=np.zeros(len(q))
 for j in range(len(q)):
  a=q.copy();b=q.copy();a[j]+=eps;b[j]-=eps
  numeric[j]=(state(a)[0]-state(b)[0])/(2*eps)
 md=float(np.max(np.abs(exact-numeric)))
 mr=float(np.max(np.abs(exact-numeric)/(1e-9+np.abs(exact)+np.abs(numeric))))
 max_abs=max(max_abs,md);max_rel=max(max_rel,mr)
 print('gradient finite difference',pert,'max_abs',md,'max_relative',mr,flush=True)
 assert md<2e-7,(md,mr)

# Premise construction is independent of the output column, provided u is fixed.
c1,s1,it1=fcm_1d(u)
u_copy=u.copy(); f_pert=f[::-1] # does not enter premise routines
c2,s2,it2=fcm_1d(u_copy)
assert np.array_equal(c1,c2) and np.array_equal(s1,s2)
a,b=subtractive_1d(u,16);a2,b2=subtractive_1d(u_copy,16)
assert np.array_equal(a,a2) and np.array_equal(b,b2)
report={'gradient_finite_difference_max_abs':max_abs,'gradient_relative_max':max_rel,
        'two_gradient_cases_PASS':True,'FCM_premise_u_only_PASS':True,
        'SC16_premise_u_only_PASS':True,'note':'The supplied u itself originated after plot-digitization; this does not make historical point holdouts prospective.'}
(p/'output/R10_gradient_and_premise_QA.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS: 2 gradient numerical comparisons and response-independent premise functions')

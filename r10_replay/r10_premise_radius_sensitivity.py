#!/usr/bin/env python3
"""Source-anchored radius sensitivity, keeping premises newly estimated from index data."""
from pathlib import Path
import numpy as np,pandas as pd
from r10_fresh_retraining import SOURCE,OUT,basis,scores,geometry
D=pd.read_csv(SOURCE);u=D.u.to_numpy();x=D.position_mm.to_numpy();f=D.force_N.to_numpy()
low=np.array([x.min(),f.min()]);scale=np.array([np.ptp(x),np.ptp(f)]);yy=(np.column_stack((x,f))-low)/scale
rows=[]
for model,k,sigma_original in [('SC9-TS',9,0.0354),('SC16-TS',16,0.0198)]:
 for kind,ra in [('count_derived',1./k),('archived_premise_width_equivalent',float(sigma_original*np.sqrt(8)))]:
  vals=np.sort(u)
  d=vals[:,None]-vals[None,:]; P=np.exp(-4*d*d/ra**2).sum(axis=1)
  centers=[];rb=1.25*ra
  for j in range(k):
   idx=int(np.argmax(P))
   if P[idx]<=0:raise ValueError('bad potential')
   c=float(vals[idx]);centers.append(c)
   P=P-P[idx]*np.exp(-4*(vals-c)**2/rb**2)
  cc=np.sort(np.asarray(centers));ss=np.full(k,ra/np.sqrt(8))
  A,_=basis(u,cc,ss);bb,_,rank,_=np.linalg.lstsq(A,yy,rcond=None)
  pred=A@bb*scale+low
  row={'model':model,'radius_specification':kind,'r_a':ra,'r_b':rb,'sigma':ss[0],
       **scores(x,f,pred,np.arange(len(x))),**geometry(x,f,pred),'rank':rank,'kappa':np.linalg.cond(A)}
  rows.append(row)
  print(model,kind,'radius',ra,'joint',row['nrmse_joint'],'area',row['loop_area_error_pct'],'closure',row['closure_error_pct'],flush=True)
OUT.mkdir(exist_ok=True)
pd.DataFrame(rows).to_csv(OUT/'R10_SC_radius_source_anchor_sensitivity.csv',index=False)

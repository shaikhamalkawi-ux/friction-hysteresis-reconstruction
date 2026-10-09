#!/usr/bin/env python3
from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat
plt.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':8,'savefig.bbox':'tight'})
D=Path(__file__).resolve().parent
O=D/'figures';O.mkdir(exist_ok=True)
R=D/'results'
# Figure 1 - actual measured loop heterogeneity; 1 micrometre stroke source, native displacement units.
fig,axs=plt.subplots(1,3,figsize=(10.7,3.1),layout='constrained')
for ax,c in zip(axs,['AF19','AF10','AF36']):
  root=D.parent/'r7_replay'/'R3/03_Code_Replay/inputs'
  cy={'AF19':997500,'AF10':1522500,'AF36':1575000}[c]
  a=np.asarray(loadmat(root/f'{c}_{cy}cycle.mat',squeeze_me=True)['hyst'],float)
  ax.plot(a[:,0],a[:,1],linewidth=1.1)
  ax.set(title=c,xlabel='Displacement (source units)',ylabel='Friction force (N)')
  ax.grid(alpha=.20)
fig.savefig(O/'R11_Measured_Contact_Loops.pdf');plt.close(fig)
# Figure 2 - force near two training-derived reversal regions, normalized across selected conditions.
agg=pd.read_csv(R/'R11_within_reversal_area_aggregate.csv')
order=['Linear (clamped)','PCHIP (clamped)','Cubic spline (clamped)','ANFIS-TS','FCM-TS','SC16-TS','SC9-TS']
a=agg.set_index('model').loc[order]
fig,ax=plt.subplots(figsize=(8.5,3.3),layout='constrained')
x=np.arange(len(a));ax.bar(x,a['reversal_force_nrmse_mean'])
ax.set_xticks(x,[v.replace(' (clamped)','').replace('Cubic spline','Cubic') for v in a.index],rotation=18,ha='right')
ax.set_ylabel('Reversal-neighborhood force NRMSE');ax.set_xlabel('Specified model');ax.grid(axis='y',alpha=.2)
fig.savefig(O/'R11_Reversal_Comparison.pdf');plt.close(fig)
# Figure 3 - chronological transfer: difference vs unfitted early-cycle template, conditional on three cases.
pdta=pd.read_csv(R/'R11_temporal_extended_per_case.csv')
base=pdta[pdta.model=='Unfitted early-cycle template'][['condition','test_cycle','force_nrmse']].rename(columns={'force_nrmse':'template'})
m=pdta.merge(base,on=['condition','test_cycle'])
m=m[m.model!='Unfitted early-cycle template'].copy();m['excess']=m.force_nrmse-m.template
# Compare four model mean differences, not four artificially independent time series.
fig,ax=plt.subplots(figsize=(8.7,3.3),layout='constrained')
for model,subset in m.groupby('model'):
  points=subset.assign(time_rank=subset.groupby('condition')['test_cycle'].rank(method='first')).groupby('time_rank').excess.mean()
  ax.plot(points.index.astype(int),1e4*points.values,marker='o',linewidth=1.1,label=model)
ax.axhline(0,color='black',linestyle='--',linewidth=.8)
ax.set(xlabel='Later-cycle observation (1 to 5 after early training)',ylabel='Excess force NRMSE vs template (×10⁻⁴)')
ax.set_xticks(range(1,6));ax.grid(alpha=.22);ax.legend(ncol=2,frameon=False)
fig.savefig(O/'R11_CrossCycle_Excess.pdf');plt.close(fig)
# Figure 4 work completion area error shown separately from force error.
fig,ax=plt.subplots(figsize=(8.5,3.1),layout='constrained')
ax.bar(np.arange(len(a)),a['completed_area_error_pct_mean'])
ax.set_xticks(np.arange(len(a)),[v.replace(' (clamped)','').replace('Cubic spline','Cubic') for v in a.index],rotation=18,ha='right')
ax.set(ylabel='Completed-loop signed-area error (%)',xlabel='Specified model')
ax.grid(axis='y',alpha=.2)
fig.savefig(O/'R11_CompletedLoop_Area.pdf');plt.close(fig)
print('R11 FIGURES:',[p.name for p in sorted(O.glob('*.pdf'))])

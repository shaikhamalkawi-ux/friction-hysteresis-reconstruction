#!/usr/bin/env python3
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['ps.fonttype']=42
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent;O=ROOT/'figures';O.mkdir(exist_ok=True)
rows=pd.read_csv(ROOT/'R7_bootstrap_7_models.csv')
models=['Linear (clamped)','PCHIP (clamped)','Cubic spline (clamped)','ANFIS-TS','SC16-TS','FCM-TS','SC9-TS']
short=['Linear','PCHIP','Natural cubic','ANFIS–TS','SC16–TS','FCM–TS','SC9–TS']
for metric,name,title in [('nrmse_pair','R7_Fig1_heldout_joint','(a) Joint displacement–force error'),('nrmse_force','R7_Fig2_heldout_force','(b) Force-only error')]:
 s=rows[(rows.metric==metric)&(rows.model.isin(models))].set_index('model').loc[models]
 fig,ax=plt.subplots(figsize=(7.6,3.4));y=np.arange(len(models));m=s['mean'].to_numpy();l=s.ci_95_lower.to_numpy();h=s.ci_95_upper.to_numpy();
 ax.errorbar(m,y,xerr=np.array([m-l,h-m]),fmt='o',capsize=3.4,markersize=5.5,color='#23364b',elinewidth=1.15)
 ax.set_yticks(y,short);ax.invert_yaxis();ax.set_xlabel('Held-out NRMSE (lower is better)');ax.set_title(title,loc='left',fontsize=11)
 ax.set_xlim(left=-0.00015);ax.grid(axis='x',alpha=.17);ax.tick_params(labelsize=9)
 fig.tight_layout();fig.savefig(O/(name+'.pdf'),bbox_inches='tight');fig.savefig(O/(name+'.png'),dpi=190,bbox_inches='tight');plt.close(fig)
# plot per-condition comparison: linear versus best fuzzy JOINT
x=pd.read_csv(ROOT/'R7_seven_models_per_condition.csv')
t=x.pivot(index='experiment',columns='model',values='nrmse_pair')
fig,ax=plt.subplots(figsize=(5.3,3.4))
xx=t['Linear (clamped)'].to_numpy();yy=t['ANFIS-TS'].to_numpy();ax.scatter(xx,yy,color='#2b5878',s=46)
for i,name in enumerate(t.index):ax.annotate(name,(xx[i],yy[i]),xytext=(3,3),textcoords='offset points',fontsize=7.5)
lim=max(np.max(xx),np.max(yy))*1.07;ax.plot([0,lim],[0,lim],'--',color='#777777',lw=0.8)
ax.set_xlim(0,lim);ax.set_ylim(0,lim);ax.set_aspect('equal',adjustable='box');ax.set_xlabel('Linear interpolation joint NRMSE');ax.set_ylabel('ANFIS–TS joint NRMSE')
ax.grid(alpha=.14);fig.tight_layout();fig.savefig(O/'R7_Fig3_paired_conditions.pdf',bbox_inches='tight');fig.savefig(O/'R7_Fig3_paired_conditions.png',dpi=190,bbox_inches='tight');plt.close(fig)
# retrospective cycle transfer summary
z=pd.read_csv(ROOT/'R3/02_Results_and_Methods/R3_cross_cycle_extended_summary.csv')
a=z[(z.scenario=='early_to_late')].copy();order=['unfitted_phase_template','FCM-TS','SC9-TS','SC16-TS','ANFIS-TS']
a=a.set_index('model').loc[order]
fig,ax=plt.subplots(figsize=(5.9,3.1));ax.barh(range(len(a)),a.mean_force_nrmse.values,color='#4c728a');ax.set_yticks(range(len(a)),['Unfitted template','FCM–TS','SC9–TS','SC16–TS','ANFIS–TS']);ax.invert_yaxis();ax.set_xlim(0,.105);ax.set_xlabel('Mean force NRMSE across three conditions');ax.set_title('Early-cycle fitting → late-cycle evaluation',fontsize=10.5)
fig.tight_layout();fig.savefig(O/'R7_Fig4_cycle_transfer.pdf',bbox_inches='tight');fig.savefig(O/'R7_Fig4_cycle_transfer.png',dpi=190,bbox_inches='tight');plt.close(fig)
# user source-derived illustrative digitized loop, NOT raw measurements
q=pd.read_csv(ROOT/'Presliding_Figure_Code_Package/data/Hysteresis_600pts.csv')
fig,ax=plt.subplots(figsize=(5.6,3.6));ax.plot(q.position_mm,q.force_N,linewidth=1.15,color='#2d5869');ax.set_xlabel('Digitized displacement (mm)');ax.set_ylabel('Digitized friction force (N)');ax.grid(alpha=.2)
fig.tight_layout();fig.savefig(O/'R7_Fig0_600point_trace.pdf',bbox_inches='tight');fig.savefig(O/'R7_Fig0_600point_trace.png',dpi=190,bbox_inches='tight');plt.close(fig)
print('FIGURE QA:',len(list(O.glob('*.pdf'))),'vector PDFs,',len(list(O.glob('*.png'))),'PNG previews')

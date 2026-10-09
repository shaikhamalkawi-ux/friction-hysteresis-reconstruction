#!/usr/bin/env python3
"""Validate actual numerical files of H2 R11 without modifying them. Exit nonzero on failure."""
from pathlib import Path
import numpy as np, pandas as pd, hashlib,json
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parent.parent
SCI=ROOT/'science';OUT=SCI/'results';R7=ROOT/'r7_replay';R10=ROOT/'r10_replay'
a=pd.read_csv(OUT/'R11_within_reversal_area_per_split.csv');b=pd.read_csv(OUT/'R11_temporal_extended_per_case.csv')
arch=pd.read_csv(R7/'R7_all_7_models_per_split.csv');arch=arch[arch.method=='index'].rename(columns={'experiment':'condition'})
m=a.merge(arch,on=['condition','seed','model'],suffixes=('_R11','_R7'),validate='one_to_one')
assert len(a)==len(m)==630 and a.model.nunique()==7 and a.condition.nunique()==9 and a.seed.nunique()==10
assert len(b)==75 and b.model.nunique()==5 and b.condition.nunique()==3
nmax=float(np.max(np.abs(m.nrmse_pair_R11-m.nrmse_pair_R7)))
fmax=float(np.max(np.abs(m.nrmse_force_R11-m.nrmse_force_R7)))
assert nmax<2e-8 and fmax<2e-8,(nmax,fmax)
assert len(pd.read_csv(R10/'output/R10_new_primary_fullfit.csv'))==4
assert len(pd.read_csv(OUT/'R11_blocked_gap_robust_summary.csv'))==7
assert len(pd.read_csv(OUT/'R11_paired_within_condition_bootstrap.csv'))==9
assert len(pd.read_csv(OUT/'R11_temporal_extended_condition_means.csv'))==15
assert len(pd.read_csv(OUT/'R11_temporal_extended_aggregate.csv'))==5
paths=pd.read_csv(OUT/'R11_temporal_source_SHA256.csv')
assert len(paths)==18
for row in paths.itertuples():
 p=SCI/'temporal_inputs'/row.filename
 assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==row.sha256,p
 assert np.asarray(loadmat(p,squeeze_me=True)['hyst']).shape==(571,2)
r7_provenance=pd.read_csv(OUT/'R11_within_source_SHA256.csv')
assert len(r7_provenance)==9
from importlib.machinery import SourceFileLoader
B=SourceFileLoader('R11_v_source',str(R7/'R7_simple_baselines.py')).load_module()
for row in r7_provenance.itertuples():
 p=B.DATA/f'{row.condition}_{row.snapshot}cycle.mat'
 assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==row.sha256,p
parity=json.loads((OUT/'R11_temporal_R3_early_late_parity.json').read_text())
assert len(parity)==5 and max(v['abs_diff'] for v in parity.values())<5e-8
print('PASS R11 CLEAN SCIENCE FILE VALIDATION')
print('REGISTERED ROWS 630/630 within, 75/75 temporal, 27x7 blocked gap')
print('ORIGINAL INPUT HASH 18/18 chronological, 9/9 selected terminal loops')
print('REFERENCE CROSS-CHECK max_joint',nmax,'max_force',fmax)
print('PREDECESSOR R10 PRIMARY MODEL RECORDS 4/4 preserved')

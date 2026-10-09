#!/usr/bin/env python3
"""Data-free checks for derived published model result tables (no raw data required)."""
from pathlib import Path
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[1]
within=pd.read_csv(R/'science/results/R11_within_reversal_area_per_split.csv')
assert len(within)==630 and within.condition.nunique()==9 and within.model.nunique()==7
assert within.groupby(['condition','seed']).size().eq(7).all()
assert within.groupby('condition').seed.nunique().eq(10).all()
assert within[['nrmse_pair','nrmse_force','reversal_force_nrmse','completed_area_error_pct']].notna().all().all()
chron=pd.read_csv(R/'science/results/R11_temporal_extended_per_case.csv')
assert len(chron)==75 and chron.condition.nunique()==3 and chron.model.nunique()==5
assert chron.groupby(['condition','test_cycle']).size().eq(5).all()
assert chron.groupby('condition').test_cycle.nunique().eq(5).all()
for file,n in [('R11_within_source_SHA256.csv',9),('R11_temporal_source_SHA256.csv',18)]:
    x=pd.read_csv(R/'science/results'/file)
    assert len(x)==n and x.sha256.str.fullmatch(r'[0-9a-f]{64}').all()
summary=pd.read_csv(R/'r7_replay/R7_seven_models_aggregate.csv')
assert len(summary)==7
assert summary.nrmse_pair.min()>0
assert summary.sort_values('nrmse_pair').model.iloc[0]=='Linear (clamped)'
fresh=pd.read_csv(R/'r10_replay/output/R10_new_primary_fullfit.csv')
assert len(fresh)==4 and set(fresh.model)=={'FCM-TS','SC9-TS','SC16-TS','ANFIS-TS'}
assert np.isfinite(fresh.nrmse_joint).all()
print('PASS: 630 within-loop records; 75 temporal records; 9+18 registered input hashes; 7 model aggregate; 4 primary fresh fits.')
print('NOTE: This is a data-free consistency check, NOT a numerical rerun of training.')

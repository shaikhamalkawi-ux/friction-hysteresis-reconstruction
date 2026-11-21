#!/usr/bin/env python3
"""56HA public derived-score verification: no laboratory measurement files needed.

Compares frozen expected Imperial27 numerical score tables with archived prior G5
source-replayed score tables, and recomputes published 27-condition selector
statistics from those frozen scores. This is NOT a fresh raw-source rerun.
"""
from __future__ import annotations
import csv, hashlib, json, math, statistics, sys
from pathlib import Path
from collections import defaultdict
ROOT=Path(__file__).resolve().parent
PRIMARY='P1_Imperial27_SevenModel_FourMasks.csv'
TEN='P1_Imperial27_PostOutcome_TenSeed_PointMasks.csv'
FILES=[PRIMARY,TEN]
KEY=('IC','mask','model')
NUMERIC=('E_joint','E_force','E_x')
CATEGORICAL=('status','optimizer_success','optimizer_iter','optimizer_evaluations',
             'optimizer_termination','n_observed','n_test','cycle_id','cohort')

def load(path):
    with path.open(newline='',encoding='utf-8') as f:
        rows=list(csv.DictReader(f))
    by_key={tuple(r[k] for k in KEY):r for r in rows}
    if len(rows)!=len(by_key):raise ValueError(f'Duplicate keys in {path.name}')
    return rows,by_key

def sha_check():
    saved=json.loads((ROOT/'RELEASE_SHA256.json').read_text(encoding='utf-8'))
    if len(saved)!=7:raise ValueError('Expected seven immutable released input files')
    for rel,h in saved.items():
        p=ROOT/rel
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=h:
            raise ValueError(f'Input SHA-256 mismatch or missing file: {rel}')
    return len(saved)

def compare():
    max_abs={x:0.0 for x in NUMERIC};sizes=[]
    for name,n in [(PRIMARY,756),(TEN,1890)]:
        a,ea=load(ROOT/'expected'/name)
        b,eb=load(ROOT/'replayed'/name)
        if len(a)!=n or len(b)!=n or ea.keys()!=eb.keys():
            raise ValueError(f'Key or row mismatch in {name}')
        for key in ea:
            x,y=ea[key],eb[key]
            for col in CATEGORICAL:
                if x[col]!=y[col]:raise ValueError(f'Categorical mismatch: {name} {key} {col}')
            for col in NUMERIC:
                xv,yv=float(x[col]),float(y[col])
                if not math.isfinite(xv) or not math.isfinite(yv):
                    raise ValueError(f'Nonfinite score in {name} {key} {col}')
                delta=abs(xv-yv)
                max_abs[col]=max(max_abs[col],delta)
                if not math.isclose(xv,yv,abs_tol=1e-6,rel_tol=1e-12):raise ValueError(f'Numerical score mismatch: {name} {key} {col}: {delta}')
        sizes.append(n)
    return {'rows':sizes,'max_abs_by_score':max_abs}

def selection():
    rows,_=load(ROOT/'expected'/PRIMARY)
    ten,_=load(ROOT/'expected'/TEN)
    ics=sorted({int(r['IC']) for r in rows})
    if len(ics)!=27:raise ValueError('Not exactly 27 source conditions')
    point, gaps,ten_points=defaultdict(list),defaultdict(list),defaultdict(list)
    for r in rows:
        key=(int(r['IC']),r['model']); v=float(r['E_joint'])
        if r['task']=='point':point[key].append(v)
        elif r['task']=='gap':gaps[key].append(v)
        else: raise ValueError('Unknown mask task')
    for r in ten:
        ten_points[(int(r['IC']),r['model'])].append(float(r['E_joint']))
    changed=0;changed_ten=0;regrets=[]; by_ic={}
    for ic in ics:
        point_methods={model:values[0] for (k,model),values in point.items() if k==ic and len(values)==1}
        gap_methods={model:statistics.mean(values) for (k,model),values in gaps.items() if k==ic and len(values)==3}
        ten_methods={model:statistics.mean(values) for (k,model),values in ten_points.items() if k==ic and len(values)==10}
        if any(len(z)!=7 for z in (point_methods,gap_methods,ten_methods)):
            raise ValueError(f'Invalid number of models or masks for IC{ic}')
        pw=min(point_methods,key=point_methods.get);gw=min(gap_methods,key=gap_methods.get)
        tw=min(ten_methods,key=ten_methods.get)
        changed+=pw!=gw;changed_ten+=tw!=gw
        regret=gap_methods[pw]-gap_methods[gw];regrets.append(regret);by_ic[ic]=(pw,gw,regret)
    if changed!=20 or changed_ten!=19:raise ValueError('Frozen selection anchors changed')
    if abs(statistics.median(regrets)-0.03469866134152921)>1e-10:
        raise ValueError('Frozen median regret changed')
    if abs(max(regrets)-0.6062264174856786)>1e-10:raise ValueError('Frozen maximum regret changed')
    crows,_=load_case_rows(ROOT/'evidence/Imperial27_casewise_geometry_and_regret.csv')
    for r in crows:
        ic=int(r['IC']); v=by_ic[ic]
        if r['point_winner']!=v[0] or r['gap_oracle_hindsight']!=v[1]:
            raise ValueError(f'Case winner changed for IC{ic}')
        if abs(float(r['gap_regret'])-v[2])>1e-10:
            raise ValueError(f'Case regret changed for IC{ic}')
    return {'source_conditions':len(ics),'selection_changed':changed,
            'ten_seed_selection_changed':changed_ten,
            'median_gap_regret':statistics.median(regrets),'maximum_gap_regret':max(regrets)}

def load_case_rows(path):
    with path.open(newline='',encoding='utf-8') as f:
        rows=list(csv.DictReader(f))
    by_id={int(r['IC']):r for r in rows}
    if len(rows)!=27 or len(by_id)!=27:raise ValueError('Casewise ledger is not 27 unique IC conditions')
    return rows,by_id

def main():
    if sys.flags.optimize:raise ValueError('Do not run with -O: fixed checks must remain enabled')
    hashes=sha_check()
    scores=compare(); decisions=selection()
    with (ROOT/'provenance/provider_source_sha256.csv').open(newline='',encoding='utf-8') as f:
        source=list(csv.DictReader(f))
    if len(source)!=27 or len({x['sha256'] for x in source})!=27:
        raise ValueError('Source-provenance register incomplete')
    out={'status':'DERIVED_LEDGER_PASS','sha256_input_files':hashes,'score_table_comparison':scores,
         'decisions':decisions,'raw_source_replay_performed':False,
         'scope':'Checks two previously retained source-derived score ledgers; does not access provider measurements'}
    print(json.dumps(out,indent=2));return 0

if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,KeyError,ValueError,ZeroDivisionError) as exc:
        print('DERIVED_LEDGER_FAIL: '+str(exc),file=sys.stderr)
        raise SystemExit(2)

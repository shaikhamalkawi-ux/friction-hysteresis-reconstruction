#!/usr/bin/env python3
"""R20 secondary descriptive density sweep for ARTS sources; 9 source files, 6 non-R19 cycles.
Common outer held-out 30% per source-cycle/seed. Nested fitting budgets select
phase-spaced candidate ranks, not responses. No machine-physic claims.
"""
from pathlib import Path
import csv,json,zipfile,hashlib,os
import numpy as np
from scipy.interpolate import PchipInterpolator,CubicSpline
ROOT=Path(__file__).resolve().parents[1]
_source = os.environ.get('H2_BRFD_SOURCE_ZIP')
if not _source:
    raise SystemExit('Missing H2_BRFD_SOURCE_ZIP. See README_FIRST.md; obtain the nine original ARTS recordings from the upstream repository for local use only.')
SRC=Path(_source).expanduser()
if not SRC.is_file():
    raise SystemExit(f'Original-source ZIP not found: {SRC}. See README_FIRST.md.')
INNER='H2_Research08_ChronologyAudit_20261013/dependencies/R07/dependencies/R05/dependencies/R03/code/'
OUT=ROOT/'analysis';CYCLES=(0,1,2,4,5,6);SEEDS=(0,1,2);BUDGETS=(20,40,80,150,300,'all70pct')
METHODS=['Linear','PCHIP','Natural cubic']

def predict(u,f,tr,te,method):
 q=u[tr];val=f[tr];sort=np.argsort(q);q=q[sort];val=val[sort];target=np.clip(u[te],q[0],q[-1])
 if method=='Linear':return np.interp(target,q,val)
 if method=='PCHIP':return PchipInterpolator(q,val,extrapolate=False)(target)
 return CubicSpline(q,val,bc_type='natural',extrapolate=False)(target)

def csvwrite(p,rows):
 with p.open('w',newline='') as fp:
  wr=csv.DictWriter(fp,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
rows=[];sources=[]
with zipfile.ZipFile(SRC) as z:
 m=json.loads(z.read(INNER+'locked_source_manifest.json'))
 assert len(m['records'])==9
 for e in m['records']:
  raw=z.read(INNER+'raw/'+e['source_path'])
  sha=hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()
  assert sha==e['git_blob_sha1']
  a=np.asarray(json.loads(raw)['data'],float);N=e['samples_per_cycle'];assert a.shape==(7*N+1,6)
  for ci in CYCLES:
   f=a[ci*N:(ci+1)*N,2];span=np.ptp(f);assert span>0;u=np.arange(N)/(N-1)
   for seed in SEEDS:
    perm=np.random.default_rng(seed).permutation(N);cand=perm[:int(np.floor(.7*N))];te=perm[len(cand):]
    rank=np.argsort(cand)
    sorted_candidates=cand[rank]
    for budget in BUDGETS:
     if budget=='all70pct': tr=cand
     else:
      # Consecutive candidate rank selection is frozen independent of responses.
      chosen=np.linspace(0,len(sorted_candidates)-1,budget,dtype=int);assert len(set(chosen))==budget
      tr=sorted_candidates[chosen]
     for method in METHODS:
      pf=predict(u,f,tr,te,method)
      err=float(np.sqrt(np.mean(((pf-f[te])/span)**2)))
      rows.append({'recording':e['case'],'source_group':e['case'].split('_')[0],
      'cycle':ci,'seed':seed,'budget':budget,'n_train':len(tr),'n_test':len(te),
      'method':method,'force_nrmse':err,'source_sha1':sha})
  print(e['case'],'done',flush=True)
assert len(rows)==9*6*3*6*3
csvwrite(OUT/'cycle_out_density_scores.csv',rows)
# Parity all candidate (the R20 fixed 70% random-mask experiment)
base=list(csv.DictReader((OUT/'cycle_out_all_scores.csv').open()))
maxdiff=0
for row in rows:
 if row['budget']!='all70pct':continue
 b=[r for r in base if r['recording']==row['recording'] and int(r['cycle_index'])==row['cycle'] and r['task']=='random' and r['mask']==str(row['seed']) and r['model']==row['method']]
 assert len(b)==1
 maxdiff=max(maxdiff,abs(float(b[0]['force_nrmse'])-row['force_nrmse']))
assert maxdiff<1e-12,maxdiff
summary=[]
for budget in BUDGETS:
 for method in METHODS:
  byfile=[]
  for e in m['records']:
   rr=[r for r in rows if r['budget']==budget and r['method']==method and r['recording']==e['case']]
   assert len(rr)==18
   byfile.append(np.mean([r['force_nrmse'] for r in rr]))
  summary.append({'budget':budget,'method':method,'files':9,'cycles_per_file':6,'seeds_per_cycle':3,
  'mean_file_first_force_nrmse':float(np.mean(byfile)),
  'median_file_first_force_nrmse':float(np.median(byfile)),
  'worst_file_first_force_nrmse':float(np.max(byfile))})
csvwrite(OUT/'cycle_out_density_summary.csv',summary)
(OUT/'density_QA.json').write_text(json.dumps({'source_files_verified':9,'total_score_rows':len(rows),'parity_all70pct_max_abs_error':maxdiff,
 'experimental_units':'9 dependent source recordings, 6 cycles each',
 'protocol':'same outer withheld indices across budgets, no withheld-response tuning; phase-spaced candidate indices',
 'classification':'posthoc sensitivity; no independent physical validation','no_original_json_redistributed':True},indent=2)+'\n')
print('PARITY',maxdiff,'SCORE',len(rows))
for b in BUDGETS:
 print(b,[(r['method'],round(r['mean_file_first_force_nrmse'],6)) for r in summary if r['budget']==b])

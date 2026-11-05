#!/usr/bin/env python3
"""H2 R20: source-verified cycle-out robustness test for classical loop completion.

R19 reported only cycle index 3. This locks evaluation to cycle indices
{0,1,2,4,5,6} of the same previously known 9 ARTS recordings, i.e. not
previously unstudied original files and not independent physical experiments.
No comparison models/thresholds are selected using these scores.
Neither raw JSON files nor materialized source samples are written to outputs.
"""
from pathlib import Path
import numpy as np, csv, json, zipfile, hashlib, itertools, sys, os
from scipy.interpolate import PchipInterpolator, CubicSpline

WORK=Path(__file__).resolve().parents[1]
_source = os.environ.get('H2_BRFD_SOURCE_ZIP')
if not _source:
    raise SystemExit('Missing H2_BRFD_SOURCE_ZIP. See README_FIRST.md; obtain the nine original ARTS recordings from the upstream repository for local use only.')
SRC=Path(_source).expanduser()
if not SRC.is_file():
    raise SystemExit(f'Original-source ZIP not found: {SRC}. See README_FIRST.md.')
INNER='H2_Research08_ChronologyAudit_20261013/dependencies/R07/dependencies/R05/dependencies/R03/code/'
OUT=WORK/'analysis'
SEEDS=(0,1,2)
CYCLES=(0,1,2,4,5,6) # never evaluated in R19 loop-completion pilot, though historical H2 saw recordings
METHODS=('Linear','PCHIP','Natural cubic')
LENGTHS=(.1,.2,.3,.4)
LOCATIONS=(.15,.5,.85)

def prediction(u,vals,train,method,score):
    q=u[train];v=vals[train]; inds=np.argsort(q); q=q[inds];v=v[inds]
    ti=np.clip(u[score],q[0],q[-1]);
    if method=='Linear': return np.interp(ti,q,v)
    if method=='PCHIP': return PchipInterpolator(q,v,extrapolate=False)(ti)
    if method=='Natural cubic': return CubicSpline(q,v,bc_type='natural',extrapolate=False)(ti)
    raise ValueError(method)

def writedata(path,rows):
    if rows:
       with path.open('w',newline='') as f:
          wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)

rows=[];provenance=[]
with zipfile.ZipFile(SRC) as z:
    assert z.testzip() is None
    man=json.loads(z.read(INNER+'locked_source_manifest.json'))
    assert len(man['records'])==9 and man['repository_commit']=='0e09aa9d950760277e31c08cd93fd01f39a5af09'
    for entry in man['records']:
        raw=z.read(INNER+'raw/'+entry['source_path'])
        blob=hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()
        assert blob==entry['git_blob_sha1']
        obj=json.loads(raw)
        a=np.asarray(obj['data'],dtype=float)
        N=int(entry['samples_per_cycle'])
        assert a.shape==(7*N+1,6) and np.isfinite(a).all()
        assert obj['displacement units']=='in' and obj['friction force units']=='kip'
        assert np.all(np.diff(a[:,0])>0)
        source_group='dataset2' if entry['case'].startswith('dataset2') else 'dataset3'
        provenance.append({'source_recording':entry['case'],'source_group':source_group,'sha256':hashlib.sha256(raw).hexdigest(),'git_blob_sha1':blob,'samples_per_cycle':N,'tested_cycles':6,'training_source_cycle_R19':3,'original_commit':man['repository_commit']})
        for cycle in CYCLES:
            loop=a[cycle*N:(cycle+1)*N]
            u=np.arange(N,dtype=float)/(N-1)
            x=loop[:,1];f=loop[:,2]
            fs=np.ptp(f);xs=np.ptp(x)
            assert fs>0 and xs>0
            masks=[]
            for seed in SEEDS:
                perm=np.random.default_rng(seed).permutation(N)
                ntr=int(np.floor(.7*N))
                tr=np.sort(perm[:ntr]);te=np.sort(perm[ntr:])
                masks.append(('random',str(seed),None,None,tr,te))
            for frac,loc in itertools.product(LENGTHS,LOCATIONS):
                L=int(round(frac*N));start=int(round(loc*(N-L)))
                te=np.arange(start,start+L)
                tr=np.r_[np.arange(start),np.arange(start+L,N)]
                masks.append(('gap',f'{frac:.1f}_{loc:.2f}',frac,loc,tr,te))
            for task,mask_id,fraction,location,tr,te in masks:
                assert len(np.intersect1d(tr,te))==0
                for meth in METHODS:
                    pf=prediction(u,f,tr,meth,te);px=prediction(u,x,tr,meth,te)
                    ef=np.sqrt(np.mean(((pf-f[te])/fs)**2))
                    ex=np.sqrt(np.mean(((px-x[te])/xs)**2))
                    joint=np.sqrt(.5*(ef*ef+ex*ex))
                    rows.append({'recording':entry['case'],'source_group':source_group,'cycle_index':cycle,
                     'task':task,'mask':mask_id,'gap_fraction':fraction,'gap_location':location,
                     'train_points':len(tr),'test_points':len(te),'model':meth,'force_nrmse':float(ef),
                     'displacement_nrmse':float(ex),'joint_nrmse':float(joint),'source_git_sha1':blob})
        print(entry['case'],'done',flush=True)
assert len(rows)==9*6*(3+12)*3, len(rows)
writedata(OUT/'cycle_out_all_scores.csv',rows)
writedata(OUT/'verified_source_register.csv',provenance)
# Exact scoring parity: same curve and mask as existing R19 pilot by manually recomputing its one central gap.
with zipfile.ZipFile(SRC) as z:
 e=man['records'][0];a=np.asarray(json.loads(z.read(INNER+'raw/'+e['source_path']))['data'],float)
 N=e['samples_per_cycle'];f=a[3*N:4*N,2];u=np.arange(N)/(N-1)
 frac=.3;loc=.5;L=int(round(frac*N));start=int(round(loc*(N-L)))
 te=np.arange(start,start+L);tr=np.r_[np.arange(start),np.arange(start+L,N)]
 got=np.sqrt(np.mean(((prediction(u,f,tr,'PCHIP',te)-f[te])/np.ptp(f))**2))
 with (WORK/'reference_R19/P1_BRFD_Multigap_Classical_PerScenario.csv').open() as fp:
  check=[r for r in csv.DictReader(fp) if r['recording']==e['case'] and r['model']=='PCHIP' and r['gap_fraction']=='0.3' and r['gap_location']=='0.5']
 assert len(check)==1
 parity=abs(got-float(check[0]['force_nrmse']))
 assert parity<1e-13,parity
# file/condition first, cycle/mask summarized within recording
summary=[]
for task in ['random','gap']:
 for mask in (['0','1','2'] if task=='random' else [f'{frac:.1f}_{loc:.2f}' for frac,loc in itertools.product(LENGTHS,LOCATIONS)]):
  by=[r for r in rows if r['task']==task and r['mask']==mask]
  recording_avg={}
  for file in [p['source_recording'] for p in provenance]:
   for meth in METHODS:
    src=[r for r in by if r['recording']==file and r['model']==meth]
    assert len(src)==6
    recording_avg[file,meth]=float(np.mean([r['force_nrmse'] for r in src]))
  for meth in METHODS:
   file_scores=np.array([recording_avg[p['source_recording'],meth] for p in provenance])
   wins=sum(recording_avg[p['source_recording'],meth]==min(recording_avg[p['source_recording'],m] for m in METHODS) for p in provenance)
   summary.append({'task':task,'mask':mask,'method':meth,'recordings':9,'cycles_per_recording':6,
   'mean_recording_first_force_nrmse':float(file_scores.mean()),
   'median_recording_first_force_nrmse':float(np.median(file_scores)),
   'largest_recording_first_force_nrmse':float(file_scores.max()),
   'recording_wins':wins})
writedata(OUT/'cycle_out_condition_first.csv',summary)
# Within-mask winning order; masks are repeated within the same source cohort
mask_mean=[]
for task in ['random','gap']:
 masks=sorted(set(r['mask'] for r in summary if r['task']==task))
 for mask in masks:
  rrs=[r for r in summary if r['task']==task and r['mask']==mask]
  mean_winner=min(rrs,key=lambda r:r['mean_recording_first_force_nrmse'])
  median_winner=min(rrs,key=lambda r:r['median_recording_first_force_nrmse'])
  mask_mean.append({'task':task,'mask':mask,'mean_winner':mean_winner['method'],'median_winner':median_winner['method'],
   'cubic_mean':next(r['mean_recording_first_force_nrmse'] for r in rrs if r['method']=='Natural cubic'),
   'linear_mean':next(r['mean_recording_first_force_nrmse'] for r in rrs if r['method']=='Linear'),
   'pchip_mean':next(r['mean_recording_first_force_nrmse'] for r in rrs if r['method']=='PCHIP')})
writedata(OUT/'task_mask_winners.csv',mask_mean)
# Aggregate across 12 gaps via recording-first, to avoid pseudoreplication
source_case={}
for p in provenance:
 file=p['source_recording'];group=p['source_group']
 for task in ['random','gap']:
  for meth in METHODS:
   selected=[r['force_nrmse'] for r in rows if r['recording']==file and r['task']==task and r['model']==meth]
   source_case[file,task,meth]=float(np.mean(selected))
agg=[]
for task in ['random','gap']:
 for meth in METHODS:
  arr=np.array([source_case[p['source_recording'],task,meth] for p in provenance]);
  agg.append({'task':task,'model':meth,'recordings':9,'score_per_recording':'mean over 6 cycles and 3 masks (random) or 12 masks (gap)',
   'mean_force_nrmse':float(arr.mean()),'median_file_mean_force_nrmse':float(np.median(arr)),
   'mean_dataset2':float(np.mean([source_case[p['source_recording'],task,meth] for p in provenance if p['source_group']=='dataset2'])),
   'mean_dataset3':float(np.mean([source_case[p['source_recording'],task,meth] for p in provenance if p['source_group']=='dataset3']))})
writedata(OUT/'overall_recording_first.csv',agg)
qa={'release':'R20 source-locked cycle-out stress; derived-only','upstream_commit':man['repository_commit'],
    'source_checks':9,'source_files':9,'excluded_prior_pilot_cycle':3,'tested_cycles_per_file':6,
    'case_cycle_mask_model_records':len(rows),'random_masks':list(SEEDS),'gap_fractions':list(LENGTHS),'gap_locations':list(LOCATIONS),
    'parity_check_R19_PCHIP_mid30_abs_difference':parity,'no_raw_release':True,
    'retrospective_exploratory':True,'not_independent_physical_replicates':True,
    'main_conclusion':{'mask_winners':mask_mean,'recording_first_overall':agg}}
(OUT/'cycle_out_QA.json').write_text(json.dumps(qa,indent=2)+'\n')
print(json.dumps({'rows':len(rows),'winners':mask_mean,'overall':agg,'parity':parity},indent=2))

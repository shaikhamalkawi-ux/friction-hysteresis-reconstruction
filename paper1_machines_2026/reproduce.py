#!/usr/bin/env python3
"""Portable ARTS source-to-score replay or a clearly labelled derived-only audit.

No original measurements are included. Supply a locally acquired ARTS ZIP whose
nine JSON files match the pinned Git blobs. Original scientific scripts are
copied without modification into an isolated execution directory.
"""
from __future__ import annotations
import argparse, hashlib, json, os, shutil, subprocess, sys, zipfile
from pathlib import Path
import pandas as pd
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'source_pipeline'))
from compare import compare_frames

def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mode',choices=['raw','derived'],required=True)
    ap.add_argument('--source-zip',type=Path)
    ap.add_argument('--out',type=Path,required=True,help='New or empty output directory')
    args=ap.parse_args()
    if sys.flags.optimize: raise ValueError('Run without -O: frozen scientific assertions must remain enabled')
    if args.mode=='raw' and (args.source_zip is None or not args.source_zip.is_file()):
        raise FileNotFoundError('Raw mode requires the original local ARTS source ZIP; no saved scores are substituted')
    source_checks=[]
    if args.mode=='raw':
        # Trust the release manifest, never an identity claim supplied inside the raw ZIP.
        frozen_manifest=json.loads((HERE/'source_pipeline/frozen_code/P1_BRFD_PublicSource_9file_Manifest.json').read_text())
        prefix='H2_Research08_ChronologyAudit_20261013/dependencies/R07/dependencies/R05/dependencies/R03/code/'
        with zipfile.ZipFile(args.source_zip) as archive:
            if archive.testzip() is not None:raise ValueError('Source ZIP CRC failure')
            supplied_manifest=json.loads(archive.read(prefix+'locked_source_manifest.json'))
            if supplied_manifest!=frozen_manifest:raise ValueError('Source manifest differs from the pinned release manifest')
            for entry in frozen_manifest['records']:
                raw=archive.read(prefix+'raw/'+entry['source_path'])
                blob=hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()
                if blob!=entry['git_blob_sha1']:raise ValueError('Pinned Git-blob mismatch: '+entry['source_path'])
                source_checks.append({'recording':entry['case'],'source_path':entry['source_path'],
                                      'git_blob_sha1':blob,'sha256':hashlib.sha256(raw).hexdigest()})
    out=args.out.resolve()
    if out.exists() and any(out.iterdir()): raise FileExistsError('Output directory is nonempty; baseline overwrite is prohibited')
    out.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    qa={'mode':args.mode,'paper':'Machines','scientific_scope':'ARTS cycle-out, density and cross-geometry audits only','raw_replay':False,'source_identity':source_checks}
    if args.mode=='raw':
        w=out/'source_replay'; (w/'analysis').mkdir(parents=True)
        shutil.copytree(HERE/'source_pipeline/frozen_code',w/'code')
        (w/'reference_R19').mkdir()
        shutil.copy2(HERE/'source_pipeline/expected/P1_BRFD_Multigap_Classical_PerScenario.csv',w/'reference_R19')
        env['H2_BRFD_SOURCE_ZIP']=str(args.source_zip.resolve())
        for script in ('held_cycle_robustness.py','held_cycle_density.py'):
            result=subprocess.run([sys.executable,str(w/'code'/script)],env=env,capture_output=True,text=True)
            (out/(script+'.log')).write_text(result.stdout+'\n'+result.stderr)
            if result.returncode: raise RuntimeError(f'{script} failed; see its log. No PASS is issued.')
        qa['tables']={}
        for name in ('cycle_out_all_scores.csv','cycle_out_density_scores.csv','cycle_out_condition_first.csv','cycle_out_density_summary.csv','overall_recording_first.csv'):
            qa['tables'][name]=compare_frames(pd.read_csv(HERE/'source_pipeline/expected'/name),pd.read_csv(w/'analysis'/name),rtol=1e-12,atol=1e-12)
        qa.update(raw_replay=True,source_recordings=9,cycle_mask_model_rows=2430,density_rows=2916)
    audit=out/'cross_geometry_audit';audit.mkdir()
    shutil.copy2(HERE/'MR1_Derived_Audit/run_strengthening.py',audit)
    (audit/'inputs').mkdir()
    inp=(out/'source_replay/analysis/cycle_out_all_scores.csv') if args.mode=='raw' else (HERE/'MR1_Derived_Audit/inputs/P1_cycle_out_all_scores.csv')
    shutil.copy2(inp,audit/'inputs/P1_cycle_out_all_scores.csv')
    result=subprocess.run([sys.executable,str(audit/'run_strengthening.py'),'--paper','P1'],env=env,capture_output=True,text=True)
    (out/'cross_geometry.log').write_text(result.stdout+'\n'+result.stderr)
    if result.returncode: raise RuntimeError('Cross-geometry audit failed')
    qa['audit_tables']={}
    for name in ('P1_gap_transfer_penalty.csv','P1_leave_one_recording_out.csv','P1_recording_first_errors.csv'):
        qa['audit_tables'][name]=compare_frames(pd.read_csv(HERE/'MR1_Derived_Audit/outputs'/name),pd.read_csv(audit/'outputs'/name),rtol=1e-12,atol=1e-12)
    qa['status']='SOURCE_TO_ARTS_RESULT_PASS' if args.mode=='raw' else 'DERIVED_ONLY_PASS'
    qa['not_claimed']='Complete replay of all PoliTO/Imperial fuzzy-model, chronological or auxiliary analyses'
    (out/'REPLAY_QA.json').write_text(json.dumps(qa,indent=2)+'\n')
    print(json.dumps(qa,indent=2));return 0
if __name__=='__main__':
    try: raise SystemExit(main())
    except (ValueError,FileNotFoundError,FileExistsError,RuntimeError,zipfile.BadZipFile,KeyError) as exc:
        print('FAIL: '+str(exc),file=sys.stderr);raise SystemExit(2)

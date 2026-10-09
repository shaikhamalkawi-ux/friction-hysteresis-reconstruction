#!/usr/bin/env python3
"""Stage user-supplied original measured files by registered SHA-256. No network access."""
from __future__ import annotations
from pathlib import Path
import argparse,csv,hashlib,shutil
ROOT=Path(__file__).resolve().parents[1]
EXPECTED_DIGITIZED_SHA='6765a4d5e9486d8f6866110d1a7a5f54a9a993f87fef456d8350c71a530b607e'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()
def register(path):
    with path.open(newline='',encoding='utf8') as f:return list(csv.DictReader(f))
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,help='Local folder containing authorized Fantetti V1 .mat files')
    p.add_argument('--digitized-source',type=Path,help='Locally authorized H2 600-point digitized CSV')
    a=p.parse_args()
    if not a.source and not a.digitized_source:p.error('Specify --source and/or --digitized-source')
    if a.source:
        if not a.source.is_dir():p.error(f'Not a folder: {a.source}')
        manifests=[('temporal',ROOT/'science/results/R11_temporal_source_SHA256.csv'),('within',ROOT/'science/results/R11_within_source_SHA256.csv')]
        records=[]
        for role,path in manifests:
            if not path.is_file():raise FileNotFoundError(path)
            records +=[(role,row) for row in register(path)]
        needed={row['sha256'] for role,row in records}
        found={};scanned=0
        for file in sorted(a.source.rglob('*.mat')):
            scanned+=1
            digest=sha(file)
            if digest in needed:
                found.setdefault(digest,file)
        missing=[(role,r) for role,r in records if r['sha256'] not in found]
        if missing:
            raise SystemExit('Source-mismatch: '+str(len(missing))+' registered input references not found (scanned '+str(scanned)+' MAT files). First missing '+missing[0][1]['filename'])
        copied=0
        for role,row in records:
            if role=='temporal':
                target=ROOT/'science/temporal_inputs'/row['filename']
            else:
                target=ROOT/'r7_replay/R3/03_Code_Replay/inputs'/row['filename']
            target.parent.mkdir(parents=True,exist_ok=True)
            if target.exists() and sha(target)==row['sha256']:continue
            shutil.copyfile(found[row['sha256']],target)
            assert sha(target)==row['sha256']
            copied+=1
        print('STAGED',copied,'original files; SHA-matched source MAT files',len(found),'of',len(needed),'unique hashes')
    if a.digitized_source:
        if not a.digitized_source.is_file():p.error('Digitized source does not exist')
        if sha(a.digitized_source)!=EXPECTED_DIGITIZED_SHA:raise SystemExit('Digitized source SHA mismatch')
        targets=[ROOT/'r7_replay/Presliding_Figure_Code_Package/data/Hysteresis_600pts.csv',ROOT/'r10_replay/source_R7/Presliding_Figure_Code_Package/data/Hysteresis_600pts.csv']
        for t in targets:
            t.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(a.digitized_source,t)
        print('STAGED authorized figure-digitized source with verified SHA')
if __name__=='__main__':main()

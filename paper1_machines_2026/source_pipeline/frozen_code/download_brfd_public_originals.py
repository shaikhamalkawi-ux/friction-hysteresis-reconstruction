#!/usr/bin/env python3
"""Fetch nine licensed ORIGINAL ARTS JSON from the pinned upstream source.

Create a LOCAL-ONLY ZIP with source material for rerunning H2 R19.
Do NOT upload or redistribute this resulting raw-source ZIP as part of H2's
journal-facing research archive. All nine bytes are checked to pinned Git blobs.
Run where github.com/raw.githubusercontent.com are reachable (e.g. Colab).
"""
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import quote
import argparse,hashlib,json,zipfile,time

HERE=Path(__file__).resolve().parent
MANIFEST=HERE/'P1_BRFD_PublicSource_9file_Manifest.json'
INNER_ROOT='H2_Research08_ChronologyAudit_20261013/dependencies/R07/dependencies/R05/dependencies/R03/code/'

def fetch_one(base_path,commit,timeout):
    uri='https://raw.githubusercontent.com/ARTS-Laboratory/Dataset-Friction-Damper-with-Backlash/'+commit+'/'+quote(base_path,safe='/')
    req=Request(uri,headers={'User-Agent':'H2-R19-Scientific-SourceReplay/1.0'})
    with urlopen(req,timeout=timeout) as stream:
        data=stream.read(12_000_000)
        if len(data)>=12_000_000:raise RuntimeError('Unexpectedly large source')
    return data

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True,help='LOCAL-ONLY ZIP, never upload with H2 paper')
    ap.add_argument('--timeout',type=int,default=60)
    args=ap.parse_args();m=json.loads(MANIFEST.read_text())
    assert len(m['records'])==9 and m['repository_commit']=='0e09aa9d950760277e31c08cd93fd01f39a5af09'
    out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    if out.exists():raise FileExistsError('Refusing overwrite '+str(out))
    summaries=[]
    with zipfile.ZipFile(out,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        z.writestr(INNER_ROOT+'locked_source_manifest.json',MANIFEST.read_bytes())
        for record in m['records']:
            raw=fetch_one(record['source_path'],m['repository_commit'],args.timeout)
            actual=hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()
            if actual!=record['git_blob_sha1']:
                raise ValueError('Git source-byte mismatch on '+record['source_path'])
            array=json.loads(raw)
            assert len(array['data'])==record['records']
            z.writestr(INNER_ROOT+'raw/'+record['source_path'],raw)
            summaries.append({'case':record['case'],'git_blob_sha1':actual,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()})
    assert zipfile.ZipFile(out).testzip() is None
    print(json.dumps({'status':'9/9 pinned source bytes verified','local_private_source_zip':str(out),'zip_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'records':summaries},indent=2))
    print('Use: H2_BRFD_SOURCE_ZIP='+str(out)+' python code/held_cycle_robustness.py')
    print('Do not include original-source ZIP in public supplementary package.')
if __name__=='__main__':main()

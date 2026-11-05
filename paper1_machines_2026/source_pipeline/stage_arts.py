#!/usr/bin/env python3
"""Stage nine byte-verified ARTS files from an upstream checkout/ZIP for local replay.

This is an input-layout adapter, not a transformation of source measurements.
The output contains third-party measurements and must stay outside release ZIPs.
"""
from __future__ import annotations
import argparse, hashlib, json, os, sys, tempfile, zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
INNER = 'H2_Research08_ChronologyAudit_20261013/dependencies/R07/dependencies/R05/dependencies/R03/code/'
MANIFEST = HERE / 'frozen_code/P1_BRFD_PublicSource_9file_Manifest.json'

def stage(upstream: Path, out: Path) -> dict:
    if out.exists():
        raise FileExistsError('Refusing to overwrite existing source archive')
    manifest = json.loads(MANIFEST.read_text())
    entries = manifest['records']
    if len(entries) != 9:
        raise ValueError('Unexpected source manifest')
    contents = []
    z = None
    try:
        if upstream.is_dir():
            paths = [p for p in upstream.rglob('*.json') if p.is_file()]
            names = [(p.as_posix(), p) for p in paths]
        elif upstream.is_file() and zipfile.is_zipfile(upstream):
            z = zipfile.ZipFile(upstream)
            if z.testzip() is not None:
                raise ValueError('Corrupt upstream ZIP')
            names = [(n, n) for n in z.namelist() if not n.endswith('/')]
        else:
            raise FileNotFoundError('Supply an extracted pinned upstream repository or its ZIP')
        for entry in entries:
            source_path = entry['source_path']
            matches = [key for name, key in names if name == source_path or name.endswith('/' + source_path)]
            if len(matches) != 1:
                raise ValueError(f'Expected one unambiguous source file: {source_path}; found {len(matches)}')
            raw = z.read(matches[0]) if z else matches[0].read_bytes()
            blob = hashlib.sha1(f'blob {len(raw)}\0'.encode() + raw).hexdigest()
            if blob != entry['git_blob_sha1']:
                raise ValueError('Source Git-blob mismatch: ' + source_path)
            if len(json.loads(raw)['data']) != entry['records']:
                raise ValueError('Source record-count mismatch: ' + source_path)
            contents.append((source_path, raw, blob))
    finally:
        if z:
            z.close()
    out.parent.mkdir(parents=True, exist_ok=True)
    # All identities pass before creating the destination; atomic replace never overwrites.
    fd, temporary = tempfile.mkstemp(prefix='.arts-stage-', suffix='.zip', dir=out.parent)
    os.close(fd)
    try:
        with zipfile.ZipFile(temporary, 'w', zipfile.ZIP_DEFLATED) as archive:
            archive.writestr(INNER + 'locked_source_manifest.json', MANIFEST.read_bytes())
            for path, raw, _ in contents:
                archive.writestr(INNER + 'raw/' + path, raw)
        if out.exists():
            raise FileExistsError('Destination appeared during staging; refusing overwrite')
        os.replace(temporary, out)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return {'status': '9_SOURCE_IDENTITIES_PASS', 'commit': manifest['repository_commit'],
            'source_files': [{'source_path': path, 'git_blob_sha1': blob,
                              'sha256': hashlib.sha256(raw).hexdigest()} for path, raw, blob in contents],
            'warning': 'Local original measurements only: exclude this output ZIP from submission/public packages'}

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(stage(args.upstream, args.out), indent=2))
    return 0

if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, FileNotFoundError, FileExistsError, zipfile.BadZipFile) as exc:
        print('FAIL: ' + str(exc), file=sys.stderr)
        raise SystemExit(2)

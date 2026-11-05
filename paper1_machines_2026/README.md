# Machines MR2 — source-to-result computational companion

**New to this code?** Start with [START_HERE.md](START_HERE.md) for a clean environment, tests and source acquisition.

## Scope and boundary
This package accompanies *From Missing Points to Later Cycles: Limits of Friction-Hysteresis Reconstruction*.
MR2 preserves the MR1 scientific results. It exposes and tests the original ARTS scoring
scripts rather than introducing a new estimator. The portable entry point covers the
nine-recording ARTS cycle-out/density tests and the derived cross-geometry audit.
It does **not** certify complete source replay of every PoliTO, Imperial or sequential analysis.
See `RESULT_TO_COMMAND.csv` for the exact coverage.

## Environment
Use Python 3.11 or later in an isolated environment, and install the requirements:

```sh
python -m pip install -r requirements.txt
python verify_release.py
python tests/test_guards.py
```

`environment_tested.json` records the actual release-test environment. The lower version
bounds in `requirements.txt` are not a guarantee of numerical identity in every environment.
All reruns compare against immutable expected outputs and fail when the stated tolerance is exceeded.
Run Python without `-O` so assertions in the unchanged scientific scripts stay enabled.

## Obtain and stage the original ARTS recordings
Provider: ARTS-Laboratory/Dataset-Friction-Damper-with-Backlash.
Pinned commit: `0e09aa9d950760277e31c08cd93fd01f39a5af09`.
The nine exact paths and Git blob identities are in
`source_pipeline/frozen_code/P1_BRFD_PublicSource_9file_Manifest.json`.
No original recordings are supplied in this package.

On a machine with access to the upstream repository, the preserved acquisition helper can fetch
and stage those original bytes:

```sh
python source_pipeline/frozen_code/download_brfd_public_originals.py --out ../local-inputs/ARTS_originals.zip
```

Alternatively, use a local checkout or ZIP of the pinned upstream repository:

```sh
python source_pipeline/stage_arts.py --upstream ../ARTS-upstream-checkout --out ../local-inputs/ARTS_originals.zip
```

The staging adapter finds only the nine manifest-listed paths, requires unique matches,
validates their Git identities, and writes the legacy input layout expected by the unchanged
scripts. It does not resample, clean or modify the source measurements. Keep this output
outside all public/submission packages. Acquisition depends on upstream/network availability;
MR2 verified locally archived originals, not a new successful network download.

## Raw-source execution

```sh
python reproduce.py --mode raw --source-zip ../local-inputs/ARTS_originals.zip --out ../runs/P1_raw
```

A successful run reports `SOURCE_TO_ARTS_RESULT_PASS` in `REPLAY_QA.json`.
The release test reproduced 2,430 cycle/mask/model score rows and 2,916 density rows,
then regenerated the cross-geometry audit from the newly computed scores.
All checked numerical values matched the retained CSV values exactly in the tested environment.
The release manifest is authoritative; a manifest supplied in the raw ZIP cannot override it.
Missing/mismatched sources and nonempty output folders are errors, never fallback requests.

## Derived-only execution

```sh
python reproduce.py --mode derived --out ../runs/P1_derived
```

This mode explicitly reports `DERIVED_ONLY_PASS`: it reruns the cross-geometry audit from
saved scores and is **not** raw-source reproduction. The 108 omission summaries reuse the
same records; neither software replay nor omission makes them independent experiments.

## Historical components and rights
`MR1_Derived_Audit/` preserves the source-derived MR1 calculations. The unchanged earlier compiled companion ZIP remains in the author-controlled Google Drive release-candidate package, rather than being duplicated in GitHub.
The historical reconstruction repository and DOI/provider links are in `OFFICIAL_SOURCES.md`.
No new public DOI, complete manuscript-wide replay or online-prediction validation is claimed.
Original H2 scripts retain the MIT notice; third-party raw data licences are separate.
The original ARTS recordings are CC BY-SA 4.0 and are not redistributed here.

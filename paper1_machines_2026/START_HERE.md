# Start here — From Missing Points to Later Cycles: Limits of Friction-Hysteresis Reconstruction

This public repository contains the **author-generated research code**, source-derived numerical results and explicit reproduction tests for the linked pre-submission paper. It does **not** contain the original experimental force/displacement recordings, nor does it promise complete replay of analyses outside the documented coverage.

## Quick verification (without original laboratory measurements)

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
python verify_release.py
python tests/test_guards.py
python reproduce.py --mode derived --out ../runs/derived_check
```

The final command outputs `DERIVED_ONLY_PASS` if the included derived calculations reproduce. **It is not a test using fresh original experimental data.** The CI workflow executes these commands whenever code is pushed.

## Optional raw-source reproduction (external original files required)

Follow `OFFICIAL_SOURCES.md`, `DATA_AND_RIGHTS.md`, and the main `README.md` for exact upstream links and source-hash admission. Place the acquired source files **outside the GitHub directory** and run:

```bash
python reproduce.py --mode raw --source-zip ../local-inputs/ARTS_originals.zip --out ../runs/p1_raw
```

Incorrect or unavailable source files cause an explicit failure, never an automatic substitution of saved scores. The exact verified source-to-result coverage is described in `RESULT_TO_COMMAND.csv`. See also `PUBLIC_RELEASE_SCOPE.md`.

## Scientific integrity and rights

- **Manuscript**: *From Missing Points to Later Cycles: Limits of Friction-Hysteresis Reconstruction* (under preparation/submission; no journal DOI yet).
- `CITATION.cff` is provisional for the forthcoming software release; no Zenodo DOI has been assigned.
- Original scientific algorithms, evaluation protocols and frozen results are unchanged from the approved release candidate.
- Code is offered under MIT where authored by the H2 research team; third-party experimental data and source-derived tables retain their distinct conditions. See `DATA_AND_RIGHTS.md`.
- Previous manuscript and journal versions remain in the author's project archive. Do not interpret these checks as hardware calibration or independent physical validation.

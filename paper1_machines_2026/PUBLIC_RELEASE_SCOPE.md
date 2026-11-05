# H2 Machines: Task-Conditioned Friction-Hysteresis Reconstruction — Code and Derived Results

## Release status
Public development code snapshot, November 2026. **Not yet a GitHub Release or Zenodo record/DOI.**

## Linked manuscript
From Missing Points to Later Cycles: Limits of Friction-Hysteresis Reconstruction (pre-submission research manuscript; no journal DOI assigned).

## Evidence that can be reproduced
ARTS nine original recordings: 2,430 point/gap score rows, 2,916 density score rows, and derived cross-geometry audits. Historical PoliTO/Imperial results are retained but not fully raw-reproducible from this portable entry point.

Run `python verify_release.py` for the original unchanged manifest, then `python tests/test_guards.py`; a derived-only check can be run as documented in README.md. Raw modes require the original source archives from the provider and will stop if required files are absent.

## Rights and provenance
The authors' original Python code is under MIT. Derived analytical ledgers have separate rights depending on their source. Read `DATA_AND_RIGHTS.md` and `RIGHTS_AND_CONTENT_MAP.md`. Third-party raw measurements are **not distributed**. The source providers remain: https://doi.org/10.17632/gy587m7gx7.1 and https://github.com/ARTS-Laboratory/Dataset-Friction-Damper-with-Backlash.

## Citation metadata
`CITATION.cff` proposes the seven manuscript authors in publication order. It must be approved specifically for software credit and may require ORCIDs before publication. No Zenodo DOI is reserved or asserted.

## Other boundaries
This is an executable software/data provenance package, not a new scientific validation, independent sensor calibration, source-independent method guarantee or representation that the paper was accepted or published.

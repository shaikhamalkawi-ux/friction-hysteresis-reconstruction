# Friction hysteresis reconstruction — research reproducibility

This code accompanies an **unpublished manuscript currently under author review** on within-loop friction-hysteresis reconstruction and cross-cycle transfer. The research compares four fixed Takagi–Sugeno fuzzy-model configurations, three elementary interpolation methods, and an unfitted measured-cycle template. This repository does **not** claim to introduce a new hysteresis state model or to validate pneumatic-seal prediction on independently measured pneumatic cycles.

## Research scope

- Nine selected experimental stainless-steel contact conditions (571 observations per measured loop) are used for within-loop point completion: 10 deterministic 70/30 splits per condition, seven model configurations, 630 model/split records.
- Three distinct experimental contact conditions (AF10, AF19, AF36), six observed cycles each, are used for early-to-later-cycle analysis: five later cycles per condition, four fuzzy configurations plus the unfitted early-cycle template (75 condition/cycle/model records). Repeated cycles within each condition are **dependent**, not independent datasets.
- A 600-point *figure-digitized* pneumatic historical example is retained as a separate illustration. The original recording and historical optimizer states were unavailable. New fitted numerical results are not a reproduction of the historical reported parameters.

All numerical comparisons are specific to the fixed experimental selection, reconstruction protocol, and input definition. A normalized position along a **known, phase-aligned recorded loop** is not an online physical memory state.

## Data provenance and access

The contact measurements originate from Fantetti et al. (2024), [Data in Brief](https://doi.org/10.1016/j.dib.2024.110374) and the associated [Mendeley Data V1 archive](https://doi.org/10.17632/gy587m7gx7.1). Original MAT measurements and third-party illustrations are **not included** in this code-only candidate. Obtain the data from the original provider and follow its license and attribution requirements. The 600-point trace was manually digitized from a curve discussed by Recchia and attributed to Nouri; it is not the original sensor record and is also not redistributed here.

See [docs/DATA_ACCESS.md](docs/DATA_ACCESS.md) for the exact registered input names and SHA-256 provenance certificates, and [docs/REPRODUCE.md](docs/REPRODUCE.md) for code execution instructions.

## Quick verification without external data

```bash
python -m pip install -r requirements.txt
python tools/verify_public_tables.py
```

This checks structure, provenance manifest shape, and summary agreement in the included derived-result tables. It **does not** rerun model fitting without the separately acquired original measurements.

## Full scientific reproduction (requires original inputs)

```bash
python tools/stage_external_inputs.py --source /path/to/extracted_Mendeley_archive
python science/R11_reproduce_extended_evidence.py
python science/R11_verify_frozen_evidence.py
cd r10_replay && bash reproduce_R10_science.sh
```

The R10 historical example additionally requires a separately authorized local copy of the 600-point digitized trajectory: see `docs/DATA_ACCESS.md`. The R7 frozen cross-checks require the registered original input files, and the R3 trained fuzzy rows are retained as provenance records. See `docs/REPRODUCE.md` for details.

## Publication and licensing status

This repository contains a **public pre-release research-code snapshot uploaded at the request of the corresponding author**. Coauthor approval for the final software release and third-party rights review remain pending. No version-of-record software release, software license, Zenodo DOI, or final manuscript citation is asserted. Until the authors approve a software license, the absence of a license means no affirmative permission to reuse the software is granted by this repository. Contact the corresponding author for permission.

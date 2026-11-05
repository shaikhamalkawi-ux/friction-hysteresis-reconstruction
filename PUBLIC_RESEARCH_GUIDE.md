# Public reader guide — H2 Machines (2026)

**This GitHub repository is public.** Anyone may browse, clone, download, inspect, and run the available code without being a manuscript author or signing into GitHub. Maintainers control direct edits; others may fork and propose pull requests.

## Which code belongs to the current paper?

The **2026 paper companion** is in [paper1_machines_2026/](paper1_machines_2026/). Its manuscript is *From Missing Points to Later Cycles: Limits of Friction-Hysteresis Reconstruction* (pre-submission; no article DOI assigned).

The older `r7_replay/`, `r10_replay/`, `science/`, and root-level documentation remain available for historical reproducibility. They are **not interchangeable with the newer paper companion**. Start with [START_HERE.md](paper1_machines_2026/START_HERE.md) and [RESULT_TO_COMMAND.csv](paper1_machines_2026/RESULT_TO_COMMAND.csv).

## Download and run a complete public-derived check

```sh
git clone https://github.com/shaikhamalkawi-ux/friction-hysteresis-reconstruction.git
cd friction-hysteresis-reconstruction/paper1_machines_2026
python -m venv .venv
# Activate: source .venv/bin/activate (macOS/Linux)
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python public_quickcheck.py
```

The last command runs the immutable SHA-256 verifier, negative/admission guard tests, and saved-score audit in a temporary output directory. It must finish with `PUBLIC_DERIVED_CHECK_PASS`. This **does not certify new raw-signal reproduction**. The repository also exposes all individual commands in [README.md](paper1_machines_2026/README.md).

## What evidence is actually covered?

| Scope | Documented verification | Boundary |
| --- | --- | --- |
| ARTS recordings | Nine source recordings and 54 *dependent* measured cycles; **2,430** point/gap/mask-model rows and **2,916** density rows reconstructed from admitted raw source in the MR2 evidence | Does not transfer to every other original measurement program |
| Point masking vs contiguous gaps | Natural cubic won three dense point-mask settings, but **zero of 12** contiguous-gap settings; PCHIP won eight and Linear four, according to the frozen descriptive audit | Task-conditioned finding; not a universal model ranking |
| Derived cross-geometry audit | **12** gap scenarios and **108** leave-one-recording-out repeated summaries | These are repeated analyses of the same nine recordings, not 108 independent experiments |
| PoliTO, Imperial, and auxiliary chronological analyses | Historical materials are retained and bounded | Not all covered by this portable raw replay |

Read [the saved numerical summary](paper1_machines_2026/MR1_Derived_Audit/outputs/P1_summary.json) and [raw verification report](paper1_machines_2026/verification/portable_P1_raw.json) for exact fields and tolerances.

## Recompute from original recordings

Original ARTS records are hosted by [ARTS-Laboratory](https://github.com/ARTS-Laboratory/Dataset-Friction-Damper-with-Backlash). The [pinned source register](paper1_machines_2026/source_pipeline/frozen_code/P1_BRFD_PublicSource_9file_Manifest.json) identifies the nine required upstream files. Follow the raw acquisition/staging commands in the paper companion; keep source recordings outside the cloned repository. The original Fantetti material is at https://doi.org/10.17632/gy587m7gx7.1.

## Public availability, rights, citation and archiving

* All readers can **download ZIP** through the GitHub Code menu or `git clone`. Read permission does not grant direct push access.
* MIT applies to the authors' original software as described in [code/data rights](paper1_machines_2026/DATA_AND_RIGHTS.md). ARTS measurements are CC BY-SA 4.0 and other upstream data have separate conditions; **no original measurement JSON/MAT/DAT/TDMS files are distributed**.
* [Software citation metadata](paper1_machines_2026/CITATION.cff) is provisional pending author/software-credit confirmation. Cite original dataset providers separately.
* **No GitHub Release or Zenodo software DOI is claimed here.** Zenodo is a recommended optional archive once release metadata and rights are approved, not a prerequisite for viewing or running this repository.

Scientific claims, frozen numerical outputs, earlier repository history and original data rights are unchanged by this reader guide or its convenience script.

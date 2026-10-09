# Reproduction guide

**Python:** 3.11 or later. Run commands from the repository root in a clean virtual environment.

1. `python -m pip install -r requirements.txt`
2. `python tools/verify_public_tables.py` — data-free validation of derived-result structure and consistency.
3. Download the original Fantetti et al. Mendeley V1 dataset via the cited DOI and extract it **locally**.
4. Run `python tools/stage_external_inputs.py --source /path/to/extracted_data` to populate ignored source input directories with SHA-matched source files. When this archive version does not match the registered hashes, stop rather than silently substitute another version.
5. Run `python science/R11_reproduce_extended_evidence.py` and `python science/R11_verify_frozen_evidence.py` to recompute R11 held-out, reversal, area, temporal and blocked-gap diagnostics. The computation may require several minutes. The R3/R7 rows needed for parity are included as *derived-output* tables, not retrained on this step.
6. For the separately authorized figure-digitized pneumatic input, stage `--digitized-source /path/to/Hysteresis_600pts.csv`, then `cd r10_replay && bash reproduce_R10_science.sh` to recompute the declared *new* 600-point fitted model metrics.

The original R10 fresh-fitting protocol uses a 46-iteration ANFIS cap; the longer-budget run is a sensitivity. Do not treat failures to achieve convergence as convergence.

### Scientific limitations

- No independent measurement of pneumatic-seal pre-sliding loops is provided by the metallic-contact benchmark.
- Index-aligned repeated cycles do not establish online forecasting when phase alignment is unknown.
- Random point hold-out is distinct from contiguous-gap and complete-cycle hold-out.
- Original third-party figures, source raw MAT files, and internal authorship/permissions documents are excluded from this repository candidate.
- Historical unrecoverable numerical records are retained as history but are not the basis for the current measured-data conclusions.

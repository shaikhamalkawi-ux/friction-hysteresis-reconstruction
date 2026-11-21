# 56HA / Machines Paper A — V13 File5-EP2 public derived-score companion

**Paper:** *Choosing Reconstruction Methods for Missing Intervals in Friction-Hysteresis Tests: Accuracy and Numerical Stability* (author review; journal submission **HOLD**).

This **additive package** updates the original historical GitHub research repository without deleting or modifying its prior code, data, or releases. The scientific baseline is 56HA **V13**; EP2 was an editorial/figure-label correction only.

## What is actually in this folder?

* `expected/`: 756 frozen original 27-condition × 4-mask × 7-procedure scores, plus 1,890 pre-existing post-outcome ten-seed scores.
* `replayed/`: the archived G5 source-to-score replay outputs of exactly the same tasks, retained so the numerical identity test is independently inspectable.
* `evidence/Imperial27_casewise_geometry_and_regret.csv`: 27 condition-level descriptive selector disagreements and retrospective gap penalties.
* `provenance/provider_source_sha256.csv`: identifiers and hashes for 27 provider MATLAB files, **not** the provider measurements themselves.
* `verify_derived.py` and `RELEASE_SHA256.json`: a standard-library-only, fail-closed test for frozen input byte identity, all 2,646 scores within declared 1e-6 absolute/1e-12 relative tolerance and exact optimizer metadata, the 20/27 and 19/27 method changes, and 0.03470 median gap penalty.

The two score sets are **retained previous outputs**. A successful test does not constitute a new independent raw-source run, a new physical experiment or complete reproduction of all reported PoliTO, Imperial, ARTS and EUCENTRE analyses. The historical ARTS raw/derived portable workflow remains under [the separate MR2 software companion](../README.md).

## Run (no third-party data required)

```bash
python paper1_machines_2026/56ha_v13_ep2/verify_derived.py
```

Expected final status: `DERIVED_LEDGER_PASS`. Python 3.10+ standard library is sufficient. The verifier recomputes the published retrospective selection counts from saved scores, rather than merely printing the reported numbers. `RELEASE_SHA256.json` pins all seven input files. All score rows and failed optimizer records remain present; a source omission or a value edit is a failure.

## Raw-source replay and limits

To recompute model scores from actual measurements, obtain the original Fantetti Mendeley Data V1 archive directly from [the provider](https://doi.org/10.17632/gy587m7gx7.1), then use the separately archived original G5 adapter subject to author release approval and to the exact source hashes documented here. The **G5 raw source adapter is not represented as part of this public folder**. This publication of retained derived numerical score tables neither supplies raw `.mat` files nor creates an independent provider dataset. EUCENTRE raw/processed and other extension materials remain subject to source/author access controls.

## Interpretation and rights

Selection differs across scattered-point and contiguous-gap geometries in **20/27** registered Imperial conditions; the ten-seed sensitivity gives **19/27**. Median retrospective gap excess is approximately **0.03470** joint NRMSE. The hindsight gap benchmark uses outcomes hidden from any real deployment and **must not** be described as an operational selector.

Underlying metal-contact measurements: Fantetti et al. (2024), [Mendeley Data V1, CC BY 4.0](https://doi.org/10.17632/gy587m7gx7.1). Derived values here are attributed to that source and the research team; uploading them does **not** relicense the original source data or any other team's content under the root software MIT license. No original raw MAT, third-party ARTS JSON, EUCENTRE raw file, machine photograph or rights-uncleared figure is included.

No GitHub Release, persistent Zenodo DOI, complete original raw replay or journal acceptance is asserted. Coauthor/repository release approvals and additional figure/provenance matters remain under author-review hold. 

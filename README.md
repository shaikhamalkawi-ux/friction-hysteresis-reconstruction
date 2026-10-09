# Friction Hysteresis Reconstruction

**Repository status: under preparation; not a finalized software or data release.**

This repository has been created for the reproducibility materials associated with an unpublished study comparing fuzzy reconstruction models and simple interpolation for measured friction-hysteresis loops.

## Scientific scope

The planned analysis distinguishes:
- **Within-loop point completion** on nine selected metallic-contact experimental conditions.
- **Reversal-neighborhood and loop-area diagnostics** for those same measured loops.
- **Time-ordered cycle-transfer assessments** on three distinct experimental conditions.
- A separate **600-point figure-digitized pneumatic example** used for historical illustration, not as raw sensor data.

Four Takagi–Sugeno fuzzy configurations are compared with linear, PCHIP and cubic-spline interpolation under matched reconstruction splits. The published claim is **not** that one family is universally superior, nor that metal-contact results independently validate pneumatic-seal control.

## Current contents

This public repository currently contains only project documentation. The computational source and derived-result tables have been prepared and tested separately, but are **not yet publicly released**. The complete internal research archive, unlicensed figures, raw MAT measurements, and historical figure-digitized point series must not be uploaded here by default.

## Data provenance

- Fantetti et al. (2024), *Data in Brief*: https://doi.org/10.1016/j.dib.2024.110374
- Associated experimental Mendeley Data, V1: https://doi.org/10.17632/gy587m7gx7.1
- Nouri (2004), original pneumatic friction experiment: https://doi.org/10.1016/S0019-0578(07)60031-7

See [data provenance](docs/DATA_SOURCES.md) and [release requirements](docs/RELEASE_STATUS.md).

## Release and citation

This repository is a staging location, **not** a frozen, citable software release. No Zenodo DOI, repository tag, article acceptance, or software license is claimed. A verified code release and `CITATION.cff` will be added after coauthor and redistribution approval.

**Correspondence:** gmalkawi@hct.ac.ae.

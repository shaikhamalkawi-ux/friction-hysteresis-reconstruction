# H2 R10: one-command figure-derived recalculation

## Frozen evidence and claim boundary

R10 uses the 600-point digitized Nouri/Recchia curve originally supplied by the H2 author; it does **not** reconstruct the Nouri sensor acquisition. Do not claim independent test histories, actual physical force prediction, or reproduction of the archived Table 7 optimization state.

The primary source CSV is retained inside `source_R7/Presliding_Figure_Code_Package/data/Hysteresis_600pts.csv` and has SHA-256 `6765a4d5e9486d8f6866110d1a7a5f54a9a993f87fef456d8350c71a530b607e`.

## Reproduce with Python

Requires Python 3.11+, numpy, pandas, scipy and matplotlib.

From the root of this unpacked archive:

```bash
bash reproduce_R10_science.sh
```

This regenerates `output/R10_new_primary_fullfit.csv`, model parameters and predictions, ten train/test split scores, independent numerical QA, ANFIS 46/150 sensitivity, SC width/radius sensitivity and gradient/premise QA.

The source is figure-digitized, with a retrospective normalized path-coordinate index. Ten 70/30 point partitions within the same digitized trajectory **are not independent physically acquired tests**.

The historical source Table 7 row values are preserved in `source_R7/R7_Table7_original_vs_verified.csv` for audit comparison, but new Table 7 estimates come only from the freshly fitted models.

**Important:** In addition to reproducible headline scores, model rankings are explicitly conditional on the fixed rule counts and newly disclosed premise-construction protocol. ANFIS at the inherited cap of 46 iterations does not satisfy a declared convergence stop; the longer budget is a sensitivity, not a replacement selected by observing test outcomes.

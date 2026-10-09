# Data access and provenance

## Recorded metal–metal contact loops

Source: Fantetti, Botto, Schwingshackl and Zucca (2024), [Data in Brief](https://doi.org/10.1016/j.dib.2024.110374); [Mendeley Data, version 1](https://doi.org/10.17632/gy587m7gx7.1). Confirm current dataset license and source acknowledgements at the host before redistribution. No original measurement MAT files are distributed in this code-only repository candidate.

Required terminal snapshots (contact condition: cycle):

| Condition | Snapshot |
|---|---:|
| AF10 | 1522500 |
| AF19 | 997500 |
| AF24, AF28, AF29, AF36, AF38, AF39, AF40 | 1575000 |

For cycle transfer, use the registered cycle list in `science/results/R11_protocol_and_scope.json` (six per AF10/AF19/AF36). The source SHA-256 references are `science/results/R11_within_source_SHA256.csv` and `science/results/R11_temporal_source_SHA256.csv`. These tables identify **the exact source files**, not public redistribution rights.

The staging tool searches a locally downloaded/extracted original dataset **by source-file hash** and copies only matching locally supplied inputs into the expected `.mat` paths. It does not fetch data or reinterpret source files.

## Historical pneumatic example

The 600-point displacement/force curve was *digitized from a published figure* (Nouri's pneumatic experiment as reproduced by Recchia), not directly recorded on a laboratory test bench by the H2 authors. The locally supplied derived CSV's registered SHA-256 is:

`6765a4d5e9486d8f6866110d1a7a5f54a9a993f87fef456d8350c71a530b607e`

The 600-point CSV is **not distributed** here. If the data owner authorizes access and reuse, specify its path to `tools/stage_external_inputs.py --digitized-source /path/to/Hysteresis_600pts.csv`. This option verifies the exact SHA-256 before staging it. Existing historical Table 7 figures must not be conflated with the documented R10 new-fit calculations.

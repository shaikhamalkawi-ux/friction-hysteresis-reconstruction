H2 R7 SEVEN-MODEL EXPERIMENT: AUTHOR REVIEW

Environment: Python 3.11+, numpy, pandas, scipy, matplotlib, scipy.io.
From this evidence folder:
  python R7_simple_baselines.py
  python R7_make_figures.py
  python R7_blocked_gap_tests.py
  python R7_blocked_ridge_sensitivity.py
  python R7_primary_trace_recalculation.py

R7_simple_baselines.py reruns 270 simple-method held-out results and combines
existing frozen 360 model/split results from R3/02_Results_and_Methods.
R7 blocked test scripts can be slow and generate very large numerical errors
for unregularized, ill-conditioned local TS regressions; all outcomes are
retained as diagnostics, not discarded. Original 600-point source values
are figure digitizations, not native sensor logs.

The nine external metal-contact loops are not pneumatic seals; no physical
transport/online-control claim is supported by the experiments.

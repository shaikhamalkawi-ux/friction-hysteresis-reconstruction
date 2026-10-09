#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
python r10_fresh_retraining.py
python r10_independent_verification.py
python r10_anfis_budget_sensitivity.py
python r10_premise_radius_sensitivity.py
python r10_gradient_and_premise_tests.py
printf 'PASS R10 SCIENCE REPLAY (all source checks and independent numeric checks)\n'

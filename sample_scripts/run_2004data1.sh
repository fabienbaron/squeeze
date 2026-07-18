#!/usr/bin/env bash
#
# 2004-data1: monochromatic image reconstruction on a 64x64 grid.
#
# A single-wavelength dataset (V2 + closure phases). No parametric model is used, so this is a
# plain image reconstruction and a good first test that the build works. The default engine is
# simulated annealing; the result is written to recon_2004data1.fits.
#
# Run extra options by appending them, e.g.:
#   ./sample_scripts/run_2004data1.sh -chains 8 -i randomthr     # 8 annealing chains
#   ./sample_scripts/run_2004data1.sh -acsa                      # adaptive-cooling annealing
#
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

set -x
./bin/squeeze ./sample_data/2004-data1.oifits \
    -w 64 -s 0.2 \
    -n 200 \
    -o recon_2004data1 \
    "$@"

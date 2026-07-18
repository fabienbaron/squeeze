#!/usr/bin/env bash
#
# V1295 Aql (HD 190073): grey SPARCO image reconstruction on a 64x64 grid.
#
# Chromatic dataset over 1.5-1.7 um. The circumstellar environment is reconstructed as a single
# grey image while the star and background are handled by the SPARCO parametric model (Kluska et
# al. 2012), so all wavelength dependence lives in the model rather than the image. Settings match
# the OITOOLS.jl demo example_image_reconstruction_sparco_grey_v1295.jl:
#
#   pixel scale 0.125 mas, 64x64
#   star:        point source, flux fraction 0.48 at lambda_0, spectrum (l/l0)^-4
#   background:  over-resolved (V=0), flux fraction 0.0,       spectrum (l/l0)^-4
#   environment: grey image,          flux fraction 1-0.48-0.0, spectrum (l/l0)^-0.1
#   reference wavelength lambda_0 = 1.6e-6 m
#   regularizer: total variation, weight 1e2
#   data used:   V2 + closure phase (T3amp and the differential OI_VIS tables are disabled)
#
# The six -P values are the SPARCO parameters in their initial state, in order:
#   f_star  UD  env_indx  lambda_0  f_bg  bg_indx
# The whole point of SPARCO is to retrieve these jointly with the image, so -S gives a non-zero
# step size to every parameter that is fitted and 0 to those held fixed:
#   f_star 0.01   UD 0.02   env_indx 0.05   lambda_0 0 (fixed)   f_bg 0 (fixed)   bg_indx 0 (fixed)
# lambda_0 is fixed because it is only the normalisation reference. The background is switched off
# for this dataset (no detectable background), so f_bg is pinned at 0; with no background flux its
# spectral index bg_indx has no effect on the likelihood, so it is fixed too rather than left to
# random-walk. To fit a background instead, give f_bg and bg_indx non-zero steps, e.g. 0.01 and 0.1.
#
# Result is written to recon_v1295aql.fits. Append extra options as needed, e.g.:
#   ./sample_scripts/run_v1295aql.sh -acsa
#
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

set -x
./bin/squeeze ./sample_data/2019_v1295Aql.oifits \
    -w 64 -s 0.125 \
    -novis -not3amp \
    -tv 1e2 \
    -modeltype sparco \
    -P 0.48 0.0  -0.1 1.6e-6 0.0 -4.0 \
    -S 0.01 0.02 0.05 0.0    0.0 0.0 \
    -n 200 \
    -o recon_v1295aql \
    "$@"

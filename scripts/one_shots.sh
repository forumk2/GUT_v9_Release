#!/usr/bin/env bash
# Single-spectrum runs.  No slope coefficients, no threshold shifts, no split
# ratio on the command line: everything comes from the spectrum file.
set -euo pipefail
export PYTHONIOENCODING=utf-8
cd "$(dirname "$0")"

# A) Minimal LR content (fermions + bidoublet + Delta_R), no GUT thresholds.
python rge_scan_v5.py --spectrum spectra/lr_minimal.json --outdir out/lr_minimal

# B) Same LR content with the remaining 10_H + 126_H components as GUT-scale
#    thresholds at an assumed common mass (see the file's description).
python rge_scan_v5.py --spectrum spectra/so10_10_126_thresholds.json --outdir out/so10_10_126_thresholds

# C) Manifestly left-right symmetric content with g_L = g_R imposed at M_I.
python rge_scan_v5.py --spectrum spectra/lr_dparity_symmetric.json --outdir out/lr_dparity_symmetric

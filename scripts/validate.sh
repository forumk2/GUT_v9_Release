#!/usr/bin/env bash
# Reproducibility check: the beta-derivation unit tests must pass and the
# three example spectra must reproduce the reference numbers in REPRODUCE.md.
# This checks numerical consistency of the code, not physical correctness.
set -euo pipefail
export PYTHONIOENCODING=utf-8
cd "$(dirname "$0")"

TIMESTAMP="$(date +%Y%m%d_%H%M%S || echo now)"
OUTDIR="results/validate/${TIMESTAMP}"
mkdir -p "${OUTDIR}"
echo "=== VALIDATE: spectrum-derived scan reproducibility ==="
echo "Results will be saved to: ${OUTDIR}"
echo

echo "[1/3] lr_betas unit tests (SM one- and two-loop coefficients, LR one-loop, threshold sign)"
python test_lr_betas.py | tee "${OUTDIR}/tests.log"
echo

echo "[2/3] Example spectra"
for spec in lr_minimal so10_10_126_thresholds lr_dparity_symmetric; do
  echo "--- ${spec} ---"
  python rge_scan_v5.py --spectrum "spectra/${spec}.json" --outdir "${OUTDIR}/${spec}" | tee "${OUTDIR}/${spec}.log"
  echo
done

echo "[3/3] Summary"
{
  echo "# Validation summary (${TIMESTAMP}) -- numerical reproducibility only"
  for spec in lr_minimal so10_10_126_thresholds lr_dparity_symmetric; do
    echo
    echo "${spec}:"
    grep -E "Best M_I|Best M_GUT|Residual mismatch|alpha_GUT|tau_p|split 1/alpha" "${OUTDIR}/${spec}.log"
  done
} > "${OUTDIR}/summary.txt"
cat "${OUTDIR}/summary.txt"
echo
echo "Done. Evidence saved under: ${OUTDIR}"

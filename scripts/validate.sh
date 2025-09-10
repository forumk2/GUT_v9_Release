#!/usr/bin/env bash
set -euo pipefail

# Ensure UTF-8 for Windows shells
export PYTHONIOENCODING=utf-8

ROOT="$(pwd)"
TIMESTAMP="$(date +%Y%m%d_%H%M%S || echo now)"
OUTDIR="${ROOT}/results/validate/${TIMESTAMP}"
mkdir -p "${OUTDIR}"

echo "=== VALIDATE: v9 winners repro + top10/scaling (if available) ==="
echo "Results will be saved to: ${OUTDIR}"
echo

SPECTRUM="spectrum_balanced_v8.json"

# -------- Winner A (0.06% / ~5.5e35) --------
echo "[1/4] Running Winner A (b3=-5.8, b2=-1.6, bbl=6.5, split=0.9, da=+,+,-)"
python rge_scan_v5.py --two-loop --multiplets "${SPECTRUM}" \
  --b3LR -5.8 --b2L -1.6 --b2R -1.6 --bBL 6.5 --split-r 0.9 \
  --dalpha1 0.001 --dalpha2 0.001 --dalpha3 -0.001 | tee "${OUTDIR}/winner_A.log"

# stash artifacts before they get overwritten
if [[ -f out/scan.csv ]]; then cp out/scan.csv "${OUTDIR}/winner_A_scan.csv"; fi
if [[ -f out/best_fit.json ]]; then cp out/best_fit.json "${OUTDIR}/winner_A_best_fit.json"; fi

# -------- Winner B (0.05% / ~1.7e36) --------
echo
echo "[2/4] Running Winner B (b3=-5.8, b2=-1.7, bbl=6.7, split=1.0, da=+,+,-)"
python rge_scan_v5.py --two-loop --multiplets "${SPECTRUM}" \
  --b3LR -5.8 --b2L -1.7 --b2R -1.7 --bBL 6.7 --split-r 1.0 \
  --dalpha1 0.001 --dalpha2 0.001 --dalpha3 -0.001 | tee "${OUTDIR}/winner_B.log"

if [[ -f out/scan.csv ]]; then cp out/scan.csv "${OUTDIR}/winner_B_scan.csv"; fi
if [[ -f out/best_fit.json ]]; then cp out/best_fit.json "${OUTDIR}/winner_B_best_fit.json"; fi

# -------- Parse key numbers into a simple summary --------
echo
echo "[3/4] Summarizing"
{
  echo "# v9 Validation Summary (${TIMESTAMP})"
  echo
  echo "Winner A (b3=-5.8, b2=-1.6, bbl=6.5, split=0.9, da=+,+,-):"
  awk '/Best MI/{mi=$0} /Best MGUT/{mg=$0} /Residual mismatch/{re=$0} /alpha_GUT/{a=$0} /tau_p/{t=$0} END{print mi RS mg RS re RS a RS t}' "${OUTDIR}/winner_A.log"
  echo
  echo "Winner B (b3=-5.8, b2=-1.7, bbl=6.7, split=1.0, da=+,+,-):"
  awk '/Best MI/{mi=$0} /Best MGUT/{mg=$0} /Residual mismatch/{re=$0} /alpha_GUT/{a=$0} /tau_p/{t=$0} END{print mi RS mg RS re RS a RS t}' "${OUTDIR}/winner_B.log"
} > "${OUTDIR}/summary.txt"

cat "${OUTDIR}/summary.txt"

# -------- If the big sweep CSV exists, run top10 + scaling --------
CSV="out/v9/scan_v9_summary.csv"
if [[ -f "${CSV}" && -f "scripts/top10.py" ]]; then
  echo
  echo "[4/4] Top-10 + scaling check from ${CSV}"
  python scripts/top10.py "${CSV}" -N 10 --by tau \
    --scaling-check --ref-tag b3_-5.8_b2_-1.6_bbl_6.5_split_0.9_da_0.001_0.001_-0.001 \
    | tee "${OUTDIR}/top10.txt"
else
  echo
  echo "[4/4] Skipping top-10/scaling (missing ${CSV} or scripts/top10.py)."
  echo "Run your full sweep (./scripts/try_balanced_v9.sh) and rerun this validator to enable it."
fi

echo
echo "Done. Evidence saved under: ${OUTDIR}"

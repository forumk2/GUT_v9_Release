#!/usr/bin/env bash
# Run the scanner once per spectrum file and tabulate the results.
#
# There is nothing to sweep *within* a spectrum any more: the LR betas, the
# GUT-scale thresholds and the U(1) split are all fixed by the multiplet
# content.  The only thing that can vary between rows is the content itself.
set -euo pipefail
export PYTHONIOENCODING=utf-8

cd "$(dirname "$0")"
SPECTRA=("$@")
if [[ ${#SPECTRA[@]} -eq 0 ]]; then SPECTRA=(spectra/*.json); fi

OUTDIR="out"
SUMMARY="$OUTDIR/spectra_summary.csv"
mkdir -p "$OUTDIR"
echo "tag,log10MGUT,log10MI,residual_pct,alphaGUT,tau" > "$SUMMARY"

for spec in "${SPECTRA[@]}"; do
  tag="$(basename "$spec" .json)"
  echo "=== $tag ==="
  python rge_scan_v5.py --spectrum "$spec" --outdir "$OUTDIR/$tag" | tee "$OUTDIR/$tag.log"
  python - "$OUTDIR/$tag/best_fit.json" "$tag" >> "$SUMMARY" <<'EOF'
import json, sys
d = json.load(open(sys.argv[1])); b = d["best"]
if b is None:
    print(f"{sys.argv[2]},nan,nan,nan,nan,nan")
else:
    print(f"{sys.argv[2]},{b['log10_MGUT']:.4f},{b['log10_MI']:.4f},{100*b['max_frac']:.5f},{b['alpha_GUT']:.5f},{b['tau_p_eppi0_years']:.3e}")
EOF
done

echo
echo "Summary: $SUMMARY"
cat "$SUMMARY"

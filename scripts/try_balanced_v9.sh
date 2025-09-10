
#!/usr/bin/env bash
set -euo pipefail

# Ensure UTF-8 prints on Windows
export PYTHONIOENCODING=utf-8

echo "=== v9: targeted-lift sweep (two-loop, LR slopes tuned) ==="

# ---- knobs centered on your v8 best ----
B3_LIST=(-6.2 -6.0 -5.8 -5.6)          # less negative → slower α3 running → higher MGUT
B2_LIST=(-1.4 -1.5 -1.6 -1.7)          # slightly more negative → helps lift meet
BBL_LIST=(6.3 6.5 6.7)                  # gentle U(1)_{B-L} slope nudge
SPLIT_LIST=(0.8 0.9 1.0)               # heavier colored vs doublets helps τ
THRESH_PATTERNS=(
  "0.001,0.001,-0.001"
  "0.000,-0.001,0.001"
)

# ---- acceptance cuts ----
RES_MAX_PCT=0.20
LOG10_MGUT_MIN_PHASE1=15.70   # ~5.0e15 GeV → clears SK in your normalization
LOG10_MGUT_MIN_PHASE2=15.88   # ~7.5e15 GeV → ~1e35 yr target
TAU_MIN_PHASE1=1.6e34
TAU_MIN_PHASE2=1.0e35

SPECTRUM_JSON="spectrum_balanced_v8.json"   # your input spectrum
OUTDIR="out/v9"
LOGDIR="$OUTDIR/logs"
CANDDIR="$OUTDIR/candidates"
mkdir -p "$LOGDIR" "$CANDDIR"

SUMMARY="$OUTDIR/scan_v9_summary.csv"
echo "tag,log10MGUT,log10MI,residual_pct,alphaGUT,tau" > "$SUMMARY"

run_one () {
  local b3="$1" b2="$2" bbl="$3" split="$4" da1="$5" da2="$6" da3="$7"
  local tag="b3_${b3}_b2_${b2}_bbl_${bbl}_split_${split}_da_${da1}_${da2}_${da3}"
  local log="$LOGDIR/${tag}.log"

  echo "--- scan with b3LR=$b3, b2L=b2R=$b2, bBL=$bbl, split_r=$split, da=($da1,$da2,$da3) ---"

  # Run your existing scanner (adjust path if needed)
  python rge_scan_v5.py --two-loop --multiplets "$SPECTRUM_JSON"     --b3LR "$b3" --b2L "$b2" --b2R "$b2" --bBL "$bbl"     --split-r "$split" --scan-points 160     --dalpha1 "$da1" --dalpha2 "$da2" --dalpha3 "$da3"     | tee "$log"

  # ---- parse from the log ----
  local log10mgut=$(grep -E "Best MGUT" "$log" | tail -n1 | sed -E 's/.*10\^([0-9.]+).*/\1/')
  local log10mi=$(grep -E "Best MI" "$log"   | tail -n1 | sed -E 's/.*10\^([0-9.]+).*/\1/')
  local resid=$(grep -E "Residual mismatch @ MGUT:" "$log" | tail -n1 | sed -E 's/.*: *([0-9.]+)%.*/\1/')
  local alpha=$(grep -E "alpha_GUT" "$log" | tail -n1 | sed -E 's/.*alpha_GUT[^0-9]*([0-9.]+).*/\1/')
  local tau=$(grep -E "tau_p" "$log" | tail -n1 | sed -E 's/.*≈ *([0-9.e+\-]+) .*/\1/')

  echo "${tag},${log10mgut},${log10mi},${resid},${alpha},${tau}" >> "$SUMMARY"
}

# ---- wide grid sweep ----
for b3 in "${B3_LIST[@]}"; do
  for b2 in "${B2_LIST[@]}"; do
    for bbl in "${BBL_LIST[@]}"; do
      for split in "${SPLIT_LIST[@]}"; do
        for pat in "${THRESH_PATTERNS[@]}"; do
          IFS=',' read -r da1 da2 da3 <<< "$pat"
          run_one "$b3" "$b2" "$bbl" "$split" "$da1" "$da2" "$da3"
        done
      done
    done
  done
done

# ---- selections ----
# Phase 1: clears SK bound
awk -F, -v rmax="$RES_MAX_PCT" -v mcut="$LOG10_MGUT_MIN_PHASE1" -v tcut="$TAU_MIN_PHASE1" '
  NR==1 {print; next}
  { lg=$2+0; resid=$4+0; tau=$6+0; if (lg>=mcut && resid<=rmax && tau>=tcut) print }'   "$SUMMARY" > "$CANDDIR/phase1.csv"

# Phase 2: push toward 1e35 yr
awk -F, -v rmax="$RES_MAX_PCT" -v mcut="$LOG10_MGUT_MIN_PHASE2" -v tcut="$TAU_MIN_PHASE2" '
  NR==1 {print; next}
  { lg=$2+0; resid=$4+0; tau=$6+0; if (lg>=mcut && resid<=rmax && tau>=tcut) print }'   "$SUMMARY" > "$CANDDIR/phase2.csv"

echo
echo "=== v9: FAST autotune (global; strict τ ≥ 1e35) ==="
python auto_tune_fast.py --two-loop --multiplets "$SPECTRUM_JSON"   --tau-min "$TAU_MIN_PHASE2"   --mnu-target 0.05 --y-ref 0.1   --mi-min 10.8 --mi-max 13.2   --dmin -0.6 --dmax 0.6   --trials 20000

# ---- local ring sweep around the winning pattern ----
echo
echo "=== v9: local ring sweep around winning point ==="
for B3 in -5.9 -5.8 -5.7; do
  for B2 in -1.6 -1.65; do
    for BBL in 6.5 6.6; do
      for R in 0.9 1.0; do
        echo "--- ring scan with b3LR=$B3, b2L=b2R=$B2, bBL=$BBL, split_r=$R ---"
        python rge_scan_v5.py --two-loop --multiplets "$SPECTRUM_JSON"           --b3LR "$B3" --b2L "$B2" --b2R "$B2" --bBL "$BBL" --split-r "$R"           --dalpha1 0.001 --dalpha2 0.001 --dalpha3 -0.001           --mi-min 10.8 --scan-points 220
      done
    done
  done
done

echo
echo "Done. Summaries:"
echo " - $SUMMARY"
echo " - $CANDDIR/phase1.csv (clears SK)"
echo " - $CANDDIR/phase2.csv (~1e35 yr zone)"

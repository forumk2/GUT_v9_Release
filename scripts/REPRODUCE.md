# REPRODUCE: Scan Reproducibility Check

This script checks that the scanner produces the same output for the same parameter inputs across runs. It is a reproducibility check, not a physics validator. It does not validate any physics claims.

## What it checks

1. **Reproduce two reference parameter points** and archive logs and JSON output:
   - Point A: `b3=-5.8, b2L=b2R=-1.6, bBL=6.5, split_r=0.9, Δα=(+0.001, +0.001, -0.001)`
     - Expected output: residual ≈ 0.06%, τ\_p ≈ 5.5×10³⁵ yr (tree-level estimate)
   - Point B: `b3=-5.8, b2L=b2R=-1.7, bBL=6.7, split_r=1.0, Δα=(+0.001, +0.001, -0.001)`
     - Expected output: residual ≈ 0.05%, τ\_p ≈ 1.7×10³⁶ yr (tree-level estimate)

   These Δα values are free parameters passed as inputs; the check verifies numerical consistency, not physical correctness.

2. **If `out/v9/scan_v9_summary.csv` is present**:
   - Print top-10 lowest-residual points (residual ≤ 0.2%)
   - Run a τ scaling consistency check relative to Point A

3. **Save a summary** to `results/validate/`.

## What it does not check

- Whether the threshold corrections (Δα) are physically motivated
- Whether the slope coefficients correspond to a specific GUT model
- Whether the τ\_p estimates are accurate (they are tree-level only)
- Whether a "low residual" point constitutes unification in any physical sense

## Quick start

From the repo root (where `rge_scan_v5.py` lives):

```bash
bash scripts/validate.sh
```

Output is archived under `results/validate/<timestamp>/`.

## Interpreting results

A passing run means the scanner is numerically consistent: the same inputs produce the same outputs. It does not mean the scan has produced a physics prediction. The reference points (A and B) are anchors for reproducibility, not claimed results.

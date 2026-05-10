# VALIDATE: Turn v9 into a paper-grade, falsifiable result

This validator automates reproduction of your **two key winners** and captures evidence.

## What it automates
1. **Reproduce the two winners** and archive logs/JSON:
   - A: `b3=-5.8, b2L=b2R=-1.6, bBL=6.5, split_r=0.9, Δα=(+,+,−)` → residual ≈ 0.06%, τ ~ 5.5e35
   - B: `b3=-5.8, b2L=b2R=-1.7, bBL=6.7, split_r=1.0, Δα=(+,+,−)` → residual ≈ 0.05%, τ ~ 1.7e36

2. **If available** (`out/v9/scan_v9_summary.csv`):
   - Print **Top-10** by τ (residual ≤ 0.2%)
   - Run **τ scaling check** vs winner A

3. Save a **summary report** in `results/validate/`.

## Quick start
From your repo root (where `rge_scan_v5.py` lives):

```bash
bash scripts/validate.sh

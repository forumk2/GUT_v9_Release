# Two-Loop RGE Parameter Scanner

This is a parameter-search harness that scans over slope coefficients and threshold corrections (Δα), searching for points of low residual mismatch under a specified multiplet spectrum, using two-loop gauge coupling running and a tree-level proton lifetime estimator. It does not derive unification — the threshold corrections are free parameters that are tuned by the scanner, not computed from the heavy multiplet spectrum.

---

## What this is not

- **Threshold corrections (Δα) are free parameters**, scanned over a grid, not derived from the heavy spectrum. A low-residual point is a fit to the chosen Δα values, not a prediction.
- **Slope coefficients** (b3\_LR, b2L\_LR, b2R\_LR, bBL\_LR, split\_r) are also free parameters in the scan. They are inputs, not outputs.
- **Reported τ\_p values are tree-level estimates.** They do not include hadronic matrix elements, one-loop corrections to decay operators, or proper threshold matching. They should be treated as order-of-magnitude placeholders.
- **Low-residual points are fits, not predictions of unification.** Finding a point where the couplings nearly meet under the scanned parameters does not imply the underlying model predicts unification.
- **Results should not be interpreted as physical predictions** without independent derivation of the threshold corrections from the multiplet spectrum and proper matrix-element computation for proton decay.

---

## Contents

```
scripts/
  rge_scan_v5.py           # main scanner: two-loop RGE over a parameter grid
  auto_tune_fast.py        # fast local search around a seed point
  rge_core.py              # RGE integration core
  two_loop_core.py         # two-loop beta function coefficients
  top10.py                 # display top-N scan points by residual or tau
  try_balanced_v9.sh       # wide grid sweep + autotune + local ring sweep
  one_shots.sh             # single-point runs for specific parameter sets
  validate.sh              # reproducibility check script
  spectrum_balanced_v8.json  # multiplet spectrum definition used for scans
  config.json              # scanner configuration
reports/
  historical/
    v9_victory_report.md   # historical document; see note inside
README.md                  # this file
```

---

## Quick start

### Prerequisites
- Python 3.10+ with numpy
- Git Bash, WSL, or any Bash shell (Windows: use Git Bash)

### Ensure UTF-8 output (Windows)

**PowerShell**
```powershell
$env:PYTHONIOENCODING = 'utf-8'
[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new()
```

**CMD**
```cmd
chcp 65001
set PYTHONIOENCODING=utf-8
```

**Bash (Git Bash/WSL)**
```bash
export PYTHONIOENCODING=utf-8
```

### Run a grid sweep
```bash
cd scripts
chmod +x try_balanced_v9.sh
./try_balanced_v9.sh
```

### Run a single-point scan
```bash
python scripts/rge_scan_v5.py --two-loop --multiplets spectrum_balanced_v8.json \
  --b3LR -5.8 --b2L -1.6 --b2R -1.6 --bBL 6.5 --split-r 0.9 \
  --dalpha1 0.001 --dalpha2 0.001 --dalpha3 -0.001
```

### Run the reproducibility check
```bash
bash scripts/validate.sh
```
See [scripts/REPRODUCE.md](scripts/REPRODUCE.md) for details.

---

## Outputs

The scanner produces:

- `out/v9/scan_v9_summary.csv` — full grid results with columns for each free parameter and computed quantities (residual mismatch, log10(MGUT), α\_GUT, τ\_p)
- `out/v9/candidates/phase1.csv` — points passing a first-pass residual cut
- `out/v9/candidates/phase2.csv` — points in a tighter residual/tau band
- `out/best_fit.json` — parameter set and outputs for the lowest-residual point found
- `results/validate/` — archived logs from reproducibility runs

No output is labeled a "winner." Points that appear in phase2.csv are parameter combinations that produce low residual under the scan conditions.

### Inspecting scan results
```bash
# Show the 10 lowest-residual points from the scan
python scripts/top10.py out/v9/scan_v9_summary.csv -N 10 --by tau

# Tau scaling check relative to a reference point
python scripts/top10.py out/v9/scan_v9_summary.csv -N 10 --by tau \
  --scaling-check --ref-tag b3_-5.8_b2_-1.6_bbl_6.5_split_0.9_da_0.001_0.001_-0.001
```

---

## Reference scan points

These are example parameter points the scanner produces low-residual output for, under `spectrum_balanced_v8.json`, used as reproducibility anchors for the scan.

**Point A**
- b3\_LR = -5.8, b2L\_LR = -1.6, b2R\_LR = -1.6, bBL\_LR = 6.5, split\_r = 0.9
- Δα = (+0.001, +0.001, -0.001) ← free parameters
- two\_loop = True, multiplets\_used = True
- Outputs: residual ≈ 0.06%, log10(MGUT) ≈ 16.07, α\_GUT ≈ 0.0235, τ\_p ≈ 5.5×10³⁵ yr (tree-level estimate)

**Point B**
- b3\_LR = -5.8, b2L\_LR = -1.7, b2R\_LR = -1.7, bBL\_LR = 6.7, split\_r = 1.0
- Δα = (+0.001, +0.001, -0.001) ← free parameters
- two\_loop = True, multiplets\_used = True
- Outputs: residual ≈ 0.05%, log10(MGUT) ≈ 16.07, α\_GUT ≈ 0.0235, τ\_p ≈ 1.7×10³⁶ yr (tree-level estimate)

These points demonstrate that the scanner can locate regions of parameter space with sub-0.1% residual. They are not predictions of physical unification.

---

## Configuring paths

If `rge_scan_v5.py` and `auto_tune_fast.py` live in a different directory:
- Edit the shell scripts and replace `python rge_scan_v5.py` with the full path
- Or set `PYTHONPATH` so imports like `rge_core` resolve

---

## License

Licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE) for terms.

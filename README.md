
# GUT v9 Milestone (Two-Loop, LR Slopes Tuned)

This package contains a small, repo-ready bundle to **reproduce the v9 milestone** you just hit:
- **Proton lifetime** τₚ(p→e⁺π⁰) ≈ **5.5 × 10³⁵ years**
- **Residual mismatch** @ MGUT ≈ **0.06%**
- **MGUT** ≈ **1.17 × 10¹⁶ GeV**, **α_GUT ≈ 0.0235**
- Based on **two-loop running** with LR slopes tuned and tiny Δα thresholds.

> **Note:** This bundle assumes you already have your working codebase checked out with the runners:
> - `rge_scan_v5.py` (scanner)
> - `auto_tune_fast.py` (fast autotuner)
> - `spectrum_balanced_v8.json` (multiplet spectrum JSON you used for v8/v9)
> - Any modules they import: `rge_core.py`, `two_loop_core.py`, etc.
>
> If these are in another directory, see **Configuring paths** below.

---

## Contents

```
scripts/
  try_balanced_v9.sh      # wide v9 sweep + selections + FAST autotune + local ring sweep
  one_shots.sh            # two single-shot commands that jump-started MGUT and τ
reports/
  v9_victory_report.md    # milestone summary (copy for README / release notes)
README.md                 # this file
```

---

## Quick start

### 0) Prereqs
- Python 3.10+ with numpy (and any other deps your codebase needs).
- Git Bash or WSL (or any Bash shell). On Windows PowerShell, install Git Bash.
- Your working codebase with `rge_scan_v5.py`, `auto_tune_fast.py`, and `spectrum_balanced_v8.json` available.

### 1) Place this folder
Put this `gut_v9_release/` folder **inside** the same directory that contains your scanners, or set the environment variable `PYTHONPATH` / use relative paths (see configuring paths).

### 2) Ensure UTF‑8 output (Windows)
Some consoles choke on Unicode (≈, →, π). Two options:

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

### 3) Run the v9 sweep
```bash
cd scripts
chmod +x try_balanced_v9.sh
./try_balanced_v9.sh
```

Outputs:
- `out/v9/scan_v9_summary.csv` (summary of runs)
- `out/v9/candidates/phase1.csv` (clears SK)
- `out/v9/candidates/phase2.csv` (~10^35 yr zone)
- Autotune outputs in `out/v9/`

### 4) (Optional) Run the “one-shot” reproducer
These are the single commands that produced the winning MGUT and τ in our session.
```bash
cd scripts
chmod +x one_shots.sh
./one_shots.sh
```

---

## Configuring paths

If your `rge_scan_v5.py` and `auto_tune_fast.py` live elsewhere, you can:
- Edit the script and replace `python rge_scan_v5.py` with `python /full/path/to/rge_scan_v5.py`
- Or set `PYTHONPATH` so imports like `rge_core` resolve.

## Winning Commands (so far): 

python rge_scan_v5.py --two-loop --multiplets spectrum_balanced_v8.json \
  --b3LR -5.8 --b2L -1.6 --b2R -1.6 --bBL 6.5 --split-r 0.9 \
  --dalpha1 0.001 --dalpha2 0.001 --dalpha3 -0.001

---

## Winning parameters (for reference)

- **b3_LR** = -5.8
- **b2L_LR** = **-1.6** (also b2R_LR = -1.6)
- **bBL_LR** = 6.5
- **split_r** = 0.9
- **Δα thresholds** = (+0.001, +0.001, -0.001)
- **two_loop** = True
- **multiplets_used** = True

**Results:**  
- log10(MI) ≈ 10.563 → MI ≈ 3.65×10^10 GeV  
- log10(MGUT) ≈ 16.069 → MGUT ≈ 1.17×10^16 GeV  
- Residual mismatch @ MGUT: **0.06%**  
- α_GUT ≈ **0.0235**  
- τₚ(p→e⁺π⁰) ≈ **5.5×10^35 years**

---

## Remarks

---

## License
Use freely within your repo/project. If you publish, please include a short acknowledgment of the v9 milestone sweep/scaffold.

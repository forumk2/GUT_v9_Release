
# GUT v9 Milestone (Two-Loop, LR Slopes Tuned)
 
- **Proton lifetime** τₚ(p→e⁺π⁰) ≈ **5.5 × 10³⁵ years**
- **Residual mismatch** @ MGUT ≈ **0.06%**
- **MGUT** ≈ **1.17 × 10¹⁶ GeV**, **α_GUT ≈ 0.0235**
- Based on **two-loop running** with LR slopes tuned and tiny Δα thresholds.

> **Note:** This bundle assumes you already have the working codebase checked out with the runners:
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
- Python 3.10+ with numpy (and any other deps the codebase needs).
- Git Bash or WSL (or any Bash shell). On Windows PowerShell, install Git Bash.
- the working codebase with `rge_scan_v5.py`, `auto_tune_fast.py`, and `spectrum_balanced_v8.json` available.

### 1) Place this folder
Put this `gut_v9_release/` folder **inside** the same directory that contains the scanners, or set the environment variable `PYTHONPATH` / use relative paths (see configuring paths).

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

If the `rge_scan_v5.py` and `auto_tune_fast.py` live elsewhere, you can:
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

# Show the 10 best points by tau (residual ≤ 0.2%)
python scripts/top10.py out/v9/scan_v9_summary.csv -N 10 --by tau

# Same, but check tau scaling vs the winning tag
python scripts/top10.py out/v9/scan_v9_summary.csv -N 10 --by tau \
  --scaling-check --ref-tag b3_-5.8_b2_-1.6_bbl_6.5_split_0.9_da_0.001_0.001_-0.001

---

## Misc

# VERIFY: Evidence Checklist for the v9 Milestone
A rigorous, falsifiable plan to confirm (or refute) that the **v9 constants & slope/threshold pattern** robustly produce tight unification and safe proton decay.

## A. Reproducibility (core)
1. **Exact one-shot winner**  
   - Run the known winning command and capture stdout and JSON.  
   - Save to `results/winner/`.  
   - Verify: residual ≤ **0.06%**, τ ≥ **1e35 yr**, log10(MGUT) ≈ **16.07**.

2. **Full v9 grid (288)**  
   - Run `try_balanced_v9.sh`.  
   - Confirm the same winner (or a neighbor) appears in `phase2.csv`.  
   - Save `out/v9/scan_v9_summary.csv` under version control (or LFS).

## B. Scaling & Mechanism (physics sanity)
3. **τ scaling law check**  
   - For candidate points near the winner, verify:  
     \[ \tau_p \propto \frac{M_{\mathrm{GUT}}^4}{\alpha_{\mathrm{GUT}}^2} \]  
   - Compute predicted τ from (MGUT, α) relative to the winner and compare to reported τ (within ~10–20%).

4. **Slope/threshold ablations**  
   - Vary one knob at a time from the winner: b3, b2L/R, bBL, split_r, Δα pattern.  
   - Expect: easing b3 (less negative) and slightly hardening b2L/R raise MGUT; Δα=(+,+,−) lifts meet by ~0.05–0.1 dex.

## C. Robustness (review-proof)
5. **SM input uncertainties**  
   - Re-run with α3(MZ), sin²θ_W, α_em varied by ±1σ.  
   - Record band of {residual, MGUT, τ}. Winner should remain in the “safe τ / tight residual” zone.

6. **Threshold modeling**  
   - Replace Δα “fudge” with **explicit heavy multiplet thresholds** from the spectrum where possible (or bracket with plausible ranges).  
   - Verify qualitative stability of MGUT and τ.

7. **Integrator / step-size stability**  
   - Halve and double RGE step sizes (or tolerance) at two-loop.  
   - Residual and MGUT should shift only slightly; trends should persist.

8. **Scheme/option toggles** (if supported)  
   - Check 1-loop vs 2-loop to show necessity of 2-loop precision.  
   - If MS̄/DR̄ toggles exist, compare for qualitative stability.

## D. Neutrino sector consistency
9. **ν-friendly autotune**  
   - Run `auto_tune_fast.py` with `--mnu-target ~0.04` and `--y-ref ~0.007`.  
   - Show at least one point with Σmν in the cosmology-friendly range **while** keeping τ ≥ 1e35 yr and residual ≤ 0.2%.

## E. Presentation (“paper-ready”)
10. **Top-10 table**  
    - Produce a top-10 leaderboard under residual ≤0.2% sorted by τ.  
    - Include (log10 MI, log10 MGUT, τ, α_GUT, b’s, split_r, Δα).

11. **Two figures**  
    - (i) Coupling lines meeting at MGUT for winner vs a near-miss.  
    - (ii) τ vs MGUT with points from the grid; overlay τ ∝ M⁴/α² curve.

## F. Pass/Fail Criteria
- **PASS**: There exists a neighborhood of points meeting: residual ≤ 0.2%, τ ≥ 1e35 yr, and a ν-consistent configuration (with smaller y or slightly larger MI). τ scaling holds within tolerance; results stable under A–C variations.  
- **FAIL**: Winner is isolated and collapses under small input/step variations, τ scaling breaks badly, or ν consistency cannot be achieved without destroying τ/residual.

--

## License
Use freely within the repo/project. If you publish, please include a short acknowledgment of the v9 milestone sweep/scaffold.

# REPRODUCE: Reproducibility Check

`validate.sh` checks that the code is numerically consistent: the beta
derivation reproduces known coefficients, and the example spectra reproduce
the reference numbers below. It does not validate any physics claim.

## What it checks

1. **`test_lr_betas.py`** — the generic beta routine must give the textbook
   SM one-loop vector (41/10, −19/6, −7) and the full two-loop matrix exactly
   as rationals; the minimal LR content must give (−7, −3, −7/3, 11/2); the
   threshold routine must have the leading-log sign and normalisation.

2. **Example spectra** (two-loop, default `config.json`, default 400 RK4
   steps, 120-point M_I grid on log10 M_I ∈ [4, 15]):

   | spectrum | log10 M_I | log10 M_GUT | α_GUT | residual | τ_p |
   |---|---|---|---|---|---|
   | lr_minimal | 9.454 | 16.177 | 0.02205 | 0.0000% | 1.68 × 10³⁶ yr |
   | so10_10_126_thresholds | 9.484 | 16.194 | 0.02245 | 0.0000% | 1.90 × 10³⁶ yr |
   | lr_dparity_symmetric | 10.494 | 15.358 | 0.02313 | 0.0000% | 8.11 × 10³² yr |

   These numbers are outputs of the listed multiplet content; nothing was
   tuned to produce them. The zero residuals are the exactly-determined case
   explained in README.md, not evidence of anything beyond the existence of a
   solution.

3. **Summary** saved to `results/validate/<timestamp>/`.

## What it does not check

- Whether the spectrum files describe a realistic model.
- Whether the assumed heavy masses in `so10_10_126_thresholds.json` are right.
- Whether τ_p is accurate (tree-level only).
- Anything about Yukawa effects, M_I thresholds, or finite matching terms.

## Run

```bash
cd scripts
bash validate.sh
```

`results/validate/20250909_203538/` predates this refactor and records the
retired knob-based mode; it is kept only as history.

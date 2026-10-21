# Two-Loop RGE Scan for an Intermediate-Scale Left-Right Model

This repository runs two-loop gauge-coupling evolution through the chain
SM → SU(3)_c × SU(2)_L × SU(2)_R × U(1)_{B-L} → unified, for a multiplet
spectrum supplied as a JSON file, and reports the intermediate scale M_I,
the unification scale M_GUT, α_GUT, a residual mismatch, and a tree-level
proton-lifetime estimate.

**Everything that enters the calculation is either low-energy data
(`config.json`) or multiplet content (a spectrum file).** The beta
coefficients, the GUT-scale threshold corrections and the split of the
hypercharge coupling at M_I are all computed from the spectrum. There are no
slope coefficients, threshold shifts, or split ratios on the command line, and
the scanner refuses them if they are passed.

Earlier versions of this repository (see `reports/historical/` and the git
history) scanned over free LR slope coefficients and free Δα shifts. Those
low-residual points were fits. That mode has been removed, not just
deprecated.

---

## What the code does

For each trial M_I on a grid:

1. Two-loop SM running M_Z → M_I (`two_loop_core.py`, RK4).
2. Tree-level matching SM → LR at M_I (`rge_core.match_SM_to_LR_with_ratio`).
3. M_GUT from the α_2L = α_3 crossing. The one-loop analytic locator in
   `rge_core.py` seeds a secant refinement so that the crossing holds for the
   two-loop, threshold-corrected couplings.
4. Two-loop LR running M_I → M_GUT with one- and two-loop gauge coefficients
   **derived from the spectrum** (`lr_betas.py`).
5. One-loop leading-log thresholds at M_GUT from the spectrum's heavy
   multiplets and their masses (`rge_scan_v5.compute_multiplet_thresholds`):
   Δ(1/α_i) = (b_i^H / 2π) ln(M_H / M_GUT), with b_i^H derived from the
   representation.
6. The split of 1/α_1 into 1/α_2R and 1/α_BL at M_I, which low-energy data
   does not fix, is set by a model property declared in the spectrum file:
   - `"split_at_MI": "unification"` solves it from α_2R(M_GUT) = α_BL(M_GUT)
     (left-right parity broken at M_GUT, as in most SO(10) chains).
   - `"split_at_MI": "lr_parity"` imposes g_L = g_R at M_I (manifest
     left-right symmetry). Nothing is solved for.
7. Residual = max fractional deviation of (1/α_1, 1/α_2L, 1/α_3, 1/α_2R,
   1/α_BL) from their mean at M_GUT (`rge_core.mismatch_metrics_n`).
8. Tree-level τ_p estimate ∝ M_GUT⁴ / (α_GUT² α_H²).

The grid minimum of the residual is refined by a golden-section search in
log10 M_I and reported.

### How to read the residual

- In **unification** mode the unknowns (M_I, M_GUT, α_GUT, split) are as
  many as the unification conditions. A residual of zero therefore means only
  that a self-consistent solution exists inside the scanned range. The
  physical output is the **predicted** M_I, M_GUT, α_GUT and α_2R(M_I). If no
  solution exists in range, the residual is non-zero at a boundary.
- In **lr_parity** mode the system is over-determined and the residual is a
  genuine test of the content. For example, the minimal content below gives a
  1.7% residual at best under g_L = g_R, i.e. it does not unify that way.

---

## Spectrum files

Format `lr_spectrum_v1` (see the `lr_betas.py` docstring). Each field gives
its dimensions under SU(3)_c, SU(2)_L, SU(2)_R, its B-L charge (quarks 1/3,
leptons −1), its type (`weyl`, `dirac`, `complex`, `real`, or `vector` for
heavy gauge multiplets in the threshold list) and a multiplicity. Threshold
entries also carry `mass_GeV`.

```
scripts/spectra/
  lr_minimal.json               fermions (3 gen, incl. nu_R) + Phi(1,2,2,0) + Delta_R(1,1,3,2); no thresholds
  so10_10_126_thresholds.json   same LR content + remaining components of a complex 10_H and 126_H
                                of SO(10) as GUT-scale thresholds, all at an ASSUMED common 1e16 GeV
  lr_dparity_symmetric.json     fermions + Phi + Delta_L + Delta_R, g_L = g_R imposed at M_I
```

Derived LR-stage coefficients for the minimal content (printed by
`python scripts/lr_betas.py scripts/spectra/lr_minimal.json`):

```
b_i   (SU3c, SU2L, SU2R, U1BL) = (-7, -3, -7/3, 11/2)
b_ij  SU3c:  -26   9/2   9/2   1/2
      SU2L:   12     8     3   3/2
      SU2R:   12     3  80/3  27/2
      U1BL:    4   9/2  81/2  61/2
```

`test_lr_betas.py` checks that the same routine reproduces the textbook SM
one- and two-loop coefficients exactly as rationals.

---

## Results for the example spectra (two-loop, default settings)

| spectrum | split | log10 M_I | log10 M_GUT | α_GUT | residual | τ_p (tree-level) |
|---|---|---|---|---|---|---|
| lr_minimal | unification | 9.454 | 16.177 | 0.02205 (1/45.35) | 0.0000% | 1.7 × 10³⁶ yr |
| so10_10_126_thresholds | unification | 9.484 | 16.194 | 0.02245 (1/44.54) | 0.0000% | 1.9 × 10³⁶ yr |
| lr_dparity_symmetric | lr_parity | 10.494 | 15.358 | 0.02313 (1/43.24) | 0.0000% | 8.1 × 10³² yr |

Zero residuals in the first two rows are the exactly-determined case
described above; the content *predicts* those scales. The third row is also
exactly determined because its content is left-right symmetric; its M_GUT is
low enough that the tree-level τ_p falls well below the Super-Kamiokande
bound, which is the familiar problem of D-parity-conserving chains. The
one-loop run of `lr_minimal` (`--one-loop`) gives log10 M_I = 9.02 and
log10 M_GUT = 16.66.

---

## What remains approximate

- **τ_p is a tree-level dimensional estimate.** No hadronic matrix element,
  no operator running, no flavour structure, M_X = M_GUT assumed. Treat as
  order-of-magnitude.
- **Two-loop coefficients are gauge-only.** Yukawa contributions to the gauge
  betas are not included in either stage.
- **Thresholds are one-loop leading-log step functions.** No finite,
  scheme-dependent matching constants.
- **Matching at M_I is tree-level.** No thresholds from W_R, Δ_R, ν_R or the
  heavy half of the bidoublet; exactly one SM Higgs doublet is assumed light
  below M_I.
- **The example spectra make assumptions:** the heavy 10_H/126_H components in
  `so10_10_126_thresholds.json` sit at a single illustrative mass; the SO(10)/LR
  gauge bosons are assumed degenerate with M_GUT; B-L assignments follow the
  standard convention (quarks 1/3, leptons −1, Δ_R ±2, bidoublet 0).
- A low residual, where it appears, is a property of the supplied multiplet
  content and the scan range. It is not a completed unification calculation.

---

## Contents

```
scripts/
  rge_scan_v5.py        scanner (spectrum-derived mode)
  lr_betas.py           group theory + beta derivation + spectrum loader
  rge_core.py           matching, M_GUT locator, residual metric
  two_loop_core.py      two-loop RK4 integrator
  test_lr_betas.py      consistency checks
  sweep_spectra.sh      run every spectrum in spectra/ and tabulate
  one_shots.sh          the three example runs
  validate.sh           reproducibility check (see REPRODUCE.md)
  top10.py              tabulate a summary CSV
  spectra/              example spectrum files
  config.json           low-energy inputs, scan range, M_GUT bracket
  out/                  outputs of the example runs
reports/historical/     old knob-based "victory report"; not a current claim
```

## Quick start

Python 3.10+ with numpy. On Windows use Git Bash and
`export PYTHONIOENCODING=utf-8`.

```bash
cd scripts
python test_lr_betas.py
python rge_scan_v5.py --spectrum spectra/lr_minimal.json
bash sweep_spectra.sh          # all spectra -> out/spectra_summary.csv
bash validate.sh               # reproducibility check
```

Options: `--mi-min/--mi-max/--scan-points` (M_I grid), `--one-loop`,
`--rk-steps` (default 400; results unchanged to 4 significant figures between
200 and 3000), `--outdir`.

## License

Apache License 2.0. See [LICENSE](LICENSE).

# Independent audit — Expedition 001

Audit date: 2026-10-28. Baseline `8d44f8c5ca6eed3a83572471dac6f9229743ec67`;
retired scanner at its parent `b79a41f`. Findings distinguish those revisions.
Historical files and results have not been rewritten. Reproduction evidence:
[`baselines.json`](../../experiments/great_pumpkin_001/results/baselines.json).

## Confirmed findings

| ID | Revision / source | Finding, consequence, disposition |
|---|---|---|
| A1 | Retired `rge_scan_v5.lr_two_loop_defaults`, `auto_tune_fast.run_point` | The LR two-loop matrix is zero; the flag enables SM gauge two-loop and LR one-loop running. Current baseline already fixes this with representation-derived coefficients. Do not describe this as a current defect. |
| A2 | Retired `spectrum_balanced_v8.json` and CLI slopes | Raw coefficient vectors lack representations, spins and multiplicities. Some negative entries cannot describe ordinary scalar/fermion matter; unidentified vectors do not establish a completion. Low residuals are fits, not model predictions. |
| A3 | Retired threshold routine | It adds `b/(12π) ln(MG/MH)`. If `b` is the multiplet's contribution to the one-loop running coefficient, the correction to the inferred full coupling is instead `b/(2π) ln(MH/MG)`, opposite sign and six times the magnitude. The undefined old vectors preclude a physical reinterpretation by merely rescaling them. Current baseline already corrected the convention. |
| A4 | Retired residual / split | It checks only the projected hypercharge combination, L and colour. With LR one-loop running, changing the R/BL split leaves that combination invariant. Small residuals therefore do not establish four-coupling equality. Current baseline includes R and BL in the residual and solves their equality. |
| A5 | Both locators at departure | `find_MGUT_by_alpha2_eq_alpha3` clamps an out-of-bracket solution and returns the endpoint as MG. Current `Model.locate_MGUT` further clamps secant steps and can return without meeting its gap tolerance. A boundary mismatch remains visible, but the claimed crossing is false. Fixed with explicit rejection and a sign-bracketed, convergence-checked solver. |
| A6 | Current `split_at_MI` | A secant iteration can exit on iteration exhaustion/equal residuals without meeting tolerance. Now rejects unconverged/nonfinite results. The corrected scale locator also requires finite positive endpoint couplings. |
| A7 | Current README, spectrum description, scanner messages | Parity was called generically overdetermining, despite a symmetric example with zero residual. L/R-symmetric beta functions preserve equal couplings: that equality is redundant. Corrected the explanation. This changes interpretation, not reference solutions. |
| A8 | `validate.sh`, `REPRODUCE.md` | The script ran examples and printed numbers but never asserted agreement with stored references. Added an explicit numerical comparison and a test that intentionally drifts a reference. |
| A9 | `test_lr_betas.py` | Direct invocation from repository root failed on `spectra/lr_minimal.json`; the documented `cd scripts` invocation passed. Fixed the test's fixture path. Bash wrappers already change to their own directory. This was a portability limitation, not a failure of the documented command. |
| A10 | Retired autotuner output | The final/best-progress MG string uses the last trial's `res`, rather than `best`; stored best JSON remains correct. Its objective also rewards `-log10(tau)` with an arbitrary weight. Historical code is preserved, not resurrected or silently corrected. |

## Reproductions and an adversarial counterexample

The archived A and B points reproduce exactly to the displayed Python floats
using the pinned retired code, 3,000 steps, archived MI and archived parameters.
This is point reproduction; it is not a rerun of an unspecified random search.

| Point | log10 MI | log10 MG | historical fractional residual | raw AR−ABL at MG |
|---|---:|---:|---:|---:|
| A | 10.563025210084 | 16.069169434569 | 0.000577360527188 | 11.4582984540 |
| B | 10.310924369748 | 16.193300918736 | 0.000512131589962 | 18.1078933508 |

Inverse couplings are dimensionless. The raw R/BL differences precede the
arbitrary projected SM threshold shifts; the old spectrum supplies no unique
LR lifting of those shifts. Setting the split ratio to **2** leaves each old
residual unchanged to about 1e-15 while changing the individual R and BL
couplings. This is an explicit identifiability failure, not an integrator error.
The Lean `hypercharge_hidden_direction` theorem records its algebraic core.

All three October reference scans were rerun at their actual default 120-point
grid and 400 RK steps, then golden refinement. Baseline minima:

| Spectrum | log10 MI | log10 MG | fractional residual |
|---|---:|---:|---:|
| minimal | 9.4538121963 | 16.1770908192 | 1.80e-9 |
| 10+126 thresholds at 1e16 GeV | 9.4837416544 | 16.1942445773 | 9.39e-10 |
| parity-symmetric | 10.4940065691 | 15.3582322759 | 2.89e-9 |

The independent adaptive solver reaches the same physical root without a
grid or golden minimization. After the locator fixes, the three printed
references remain within their rounding tolerances. The changed behavior is
rejection of absent/unconverged crossings, not improvement of the minima.
See [`verification.json`](../../experiments/great_pumpkin_001/results/verification.json).

Concrete A5 counterexample: `A2=30, A3=10, b2=-3, b3=-7, MI=1e10`, bracket
`[14.5,17.5]`. The mathematical crossing is log10 μ ≈ 23.644; the baseline
returns 17.5 with a nonzero coupling gap. Equal slopes/equal couplings are
coincident lines, not a uniquely located crossing. The new analytic helper
returns `(None,None)` in both cases and computes in log space to avoid
exponentiating an enormous out-of-range root.

## Equation and representation checks

The current RG convention is correct for α=g²/(4π), t=ln μ:
`dα_i/dt = b_i α_i²/(2π) + α_i² Σ_j B_ij α_j/(8π²)`.
The inverse-coupling prototype independently integrates
`dA_i/dt = −b_i/(2π) − Σ_j B_ij/A_j/(8π²)`.
The generic Weyl/Dirac and real/complex scalar weights agree with the
gauge-sector terms in Luo, Wang and Xiao (2003), Eq. (30) and product-group
rules. They exactly recover the SM matrix. The LR matrix also matches
Bertolini, Di Luzio and Malinský (2009), Table VI, after reordering. These
are independent external checks, not only comparisons to local literals.
Sources and equation locations are in [the bibliography](great-pumpkin-001-sources.md).

`qBL=√(3/8)(B−L)` and `q1=√(3/5)Y` give
`Y=T3R+(B−L)/2` and `A1=(3/5)AR+(2/5)ABL`.
The squared coefficients sum to one, as needed for equal normalized parent
couplings. There is one Abelian factor in each running stage, so no pair of
independent U(1)s needing kinetic mixing in this chosen chain.

The listed complex 10+126 scalars, including their light remnants, have total
dimension 136 and summed index contribution `(12,12,12,12)`. This checks
gauge-index bookkeeping, not all SO(10) Clebsch coefficients, conjugations,
Yukawa interactions or vacuum conditions. Colour conjugation is invisible
to these beta functions. The general dimension-only SU(3) interface is
ambiguous for dimension 15 versus 15′ and truncates noninteger dimension
inputs; the selected fields use only 1,3,6,8. The loader is **not** a general
representation, anomaly or scalar-potential validator. These limitations
remain explicitly outside this expedition's small corrective patch.

## Numerical and physical limitations

- The legacy RK4 integrator has no error estimate, positivity event or general
  Landau-pole handling. It remains a historical/current utility for regular
  trajectories. The new inverse-coupling solver checks a perturbative domain
  `A_i>1`, uses adaptive errors and rejects failed integration. Tests compare
  RK4 step refinement with DOP853, RK45, an analytic solution and Brent roots.
- The current five-entry residual double-counts a derived hypercharge
  combination alongside its R and BL parents. Its zero set is appropriate,
  but its nonzero magnitude is a convention-dependent diagnostic, not a
  chi-squared or likelihood. The prototype reports three low-energy boundary
  residuals and explicitly labels that all three inputs were used to solve.
- The current heavy 10+126 example has no explicit SO(10)-breaking sector.
  The prototype accounts for uneaten 45 scalars but does not solve their
  potential. Assumed masses are not derived simply because beta functions are.
- Both versions omit Yukawa contributions in both running stages. Neither
  is a complete two-loop physical prediction. GUT leading logs are not full
  threshold matching; intermediate nondegeneracy, vector masses, finite
  constants and Yukawas require separate treatment.
- The proton expression is a normalized dimensional ansatz, not a computed
  partial decay width. It includes an αH scaling factor despite prose saying
  "no matrix element": no matrix-element calculation or operator/flavour
  normalization justifies its prefactor. `MX=MG` is an assumption.
- `mν=y²(174 GeV)²/MR` is a one-generation type-I illustration with
  `MR=MI`; in SO(10), MR depends on a Majorana Yukawa and a VEV, and type-II
  terms and flavour matrices can matter. The autotuner's neutrino target is
  part of its objective, never an independent validation datum.
- The historical report used a stale proton bound and a cosmological mass
  sum as if it were a simple single-neutrino check. Preserve the report's
  warning, and use dated, channel-specific experimental limits instead.

## What the fixes do not claim

A sign-bracketed solver does not enumerate tangential/multiple roots or roots
across invalid-coupling regions. Failure to bracket is not a theorem of global
nonexistence. The fixed parser still permits spectra the physics does not.
No historical output was recomputed under a different convention and relabeled
as the same experiment. No code change establishes unification in nature.

# What the one-loop equations actually constrain

Conventions and the spectrum are in [model.md](great-pumpkin-001-model.md).
Let `A=(A1,A2,A3)` be the measured normalized inverse couplings at MZ,
`d=PΔ_G−δ_I` the projected matching correction,
`x=ln(MI/MZ)/(2π)`, `y=ln(MG/MI)/(2π)` and `u=1/αG`.
For the fixed minimal light spectrum, write `a=A+d`. Unification is

```
a1 = u + (41/10)x + (4/5)y
a2 = u − (19/6)x − 3y
a3 = u − 7x − 7y.
```

These are affine equations, not evidence for a dynamical mechanism forcing
unification. Taking differences removes u. Their coefficient matrix is

```
C = [109/15   19/5
       23/6      4],       det C = 29/2.
```

The nonzero determinant establishes uniqueness for fixed corrected inputs:

```
x = (40a1 − 78a2 + 38a3)/145
y = (−115a1 + 333a2 − 218a3)/435
u = (35a1 + 693a2 − 293a3)/435.
```

Positive x and y give ordered scales, and u>1 is a useful endpoint domain.
These inequalities are necessary and sufficient for the **projected affine
system** with those inequalities, not for a stable scalar vacuum or a
positive full four-coupling trajectory. The Python solver checks the latter
separately. Degenerate coefficient matrices in other spectra can yield no
solution or an unidentifiable continuum; a scan cannot remove that rank issue.

## Eliminate the upper interval before scanning

The weight `w=(20,−39,19)` annihilates both the universal vector and the
projected LR beta vector. Thus

`2w·(A+d)=145x`.

The upper scale and the unified coupling cannot be used to evade this
constraint. If `2w·d≤B`, then any proposed `x≥xmin` requires
`2w·A+B≥145xmin`. Its strict violation is a no-solution certificate for
the stipulated model. This is ordinary linear elimination / dual bounding,
not new mathematics; the useful contribution is making it explicit and
checkable for this repository's spectrum.

## A bound on a whole heavy-mass family

For the 16 heavy scalar entries in `model.spectrum()` order, set
`zH=log10(MH/MG)`. The exact integers `cH=2w·P bH` are:

| Entry | cH |
|---|---:|
| 10: T, Tbar | 9, 9 |
| 126: (3,1,1), (bar3,1,1) | 9, 9 |
| 126: (6,1,3), (3,1,3) | 207, 75 |
| 126: (bar6,3,1), (bar3,3,1) | −201, −129 |
| 126: ΔL | −28 |
| 126: (8,2,2), (1,2,2) | 8, −18 |
| 126: (3,2,2), (bar3,2,2) | 14, 14 |
| 45: octet, left triplet, right triplet | 19, −26, 8 |

Hence `Σ|cH|=783`, `ΣcH=−21`, and exactly at one loop with tree-level
intermediate matching and no finite constants,

```
log10 MI = log10 MI(degenerate) + Σ cH zH / 145
|log10 MI − log10 MI(degenerate)| ≤ (783/145) log10 f
                                        = 5.4 log10 f.
```

The sign-aligned corner attains this upper bound in the affine model. At
f=3 it is an ordered, positive solution with all heavy masses above MI.
The repository inputs give a center 9.0158165511 and bounds
`[6.4393617756,11.5922713266]`. Thus MI≥1e13 GeV is incompatible with
**this one-loop factor-three scalar box**, fixed vector masses and degenerate
MI thresholds. It is not excluded for all SO(10) models, and the result is
not promoted to a rigorous two-loop bound.

If all heavy scalar masses instead share one z, the response is only
`−21z/145`; at f=3 the absolute shift is 0.06910031965 decades. Complete
degenerate SO(10) multiplets contribute a universal index; incomplete light
remnants and the uneaten 45 bookkeeping leave this small remainder.
The cancellation is expected from group theory, not a new physical invariant.

An exact affine support calculation implies that reaching log10 MI=13 needs
at least `f=10^((13−9.0158165511)/5.4)≈5.467789`. A separate linear program
also enforces ordered scales and the MG bracket and reaches that bound at
log10 MG=15.1358116865. This is a minimum permitted mass range **within the
relaxed independent-mass model**, not proof that an SO(10) potential realizes
the extremal spectrum. Two-loop running at the same masses instead gives
`(13.51528,14.68674)`, displaying substantial perturbative scale shifts.

## Why arbitrary shifts erase the test

For any chosen x,y,u, define
`d=(u,u,u)+bSM x+bLR_projected y−A`. The equations then hold identically.
Without mass-derived bounds or correlations, exact convergence has no
discriminatory content. A universal shift changes u but not x,y. An old
LR split change `(AR,ABL)→(AR+2k,ABL−3k)` preserves hypercharge exactly;
when running is affine, a projected residual cannot constrain that direction.

## Formal status

[`GreatPumpkin.lean`](../../formal/great_pumpkin/GreatPumpkin.lean) proves nine
statements over an explicit standard ordered field: the invariant, the closed
form solution, uniqueness, a positive-scale equivalence, bounded-threshold
obstruction, the 16-variable mass-box support inequality, universal-shift
invariance, arbitrary-offset fitting, and the invisible hypercharge direction.
The equations are multiplied by 30, so the proof does not approximate their
rational coefficients. A rational instantiation checks non-vacuity.

Lean 4.24.0's standard library suffices; Mathlib is not required. Both `lake
build` and a separate `lean --trust=0 GreatPumpkin.lean` process passed. Axiom
reports contain only Lean's usual `propext`, `Classical.choice`, `Quot.sound`;
there are no added axioms or admitted proofs. These foundations and the Lean
kernel are the trusted mathematical base.

The proof does **not** derive beta functions from quantum field theory,
formalize real logarithms/RG differential equations, certify floating-point
input values, or prove a two-loop enclosure. Exact rational coefficient tests
connect Python's spectrum to the integer support theorem; numerical decimal
bounds are evaluated in Python, not certified real-arithmetic interval proofs.
Formalizing that bridge and a two-loop remainder bound remains open.

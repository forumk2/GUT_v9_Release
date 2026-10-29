# Constrained gauge-sector prototype

## Candidate and assumptions

Choose the nonsupersymmetric chain
`SO(10) → SU(3)c × SU(2)L × SU(2)R × U(1)BL → SM`, with D-parity broken
at the first step. The candidate scalar sector is real 45H, complex 126H
(or its conjugate with corresponding charge conventions) and complex 10H;
three 16F generations contain νR. This is the class studied in the 2009,
2012, 2022, 2023 and 2025 sources in [the bibliography](great-pumpkin-001-sources.md).

**This expedition specifies an effective gauge-sector mass ansatz, not a
completed viable SO(10) Lagrangian point.** In particular, the published
2023 light-Higgs problem is a substantive challenge to this completion;
the 2025 discussion does not constitute a resolved scalar-potential check.
Retaining this model makes the distinction testable rather than avoiding it.

Units: masses in GeV; A=1/α dimensionless; t=ln μ. LR array order is
`(3,L,R,BL)`, SM order `(1,L,3)`. Below MI: three SM generations and one
complex Higgs doublet. Above MI: the following fields, assumed active
throughout the interval:

| Field | LR representation (colour,L,R,B−L) | Type / multiplicity |
|---|---|---|
| QL | (3,2,1,1/3) | Weyl ×3 |
| QR conjugate | (bar3,1,2,−1/3) | Weyl ×3 |
| LL | (1,2,1,−1) | Weyl ×3 |
| LR conjugate | (1,1,2,1) | Weyl ×3 |
| Φ | (1,2,2,0) | complex scalar |
| ΔR | (1,1,3,2) | complex scalar |

These fermions are anomaly-free when conjugations are interpreted as above;
the beta-function data format itself drops colour conjugation and cannot
certify anomalies. A mixed light bidoublet direction is represented by one
copy for gauge running; its existence and realistic fermion couplings are
assumptions to be tested, not outputs of this calculation.

For a field, `S_i = T_i × product(other dimensions)`. One-loop contributions
are `(2/3)S_i` for Weyl fermions, `(1/3)S_i` for complex scalars, half that
for real scalars; gauge fields contribute `−11 C2(G_i)/3`. Abelian indices
use normalized charge squared. These yield

```
bLR = (−7, −3, −7/3, 11/2)
BLR = [ −26   9/2   9/2    1/2
         12     8     3    3/2
         12     3  80/3   27/2
          4   9/2  81/2   61/2 ]
```

The gauge-only two-loop table matches the literature independently of the
SM coefficient test. SM b is `(41/10,−19/6,−7)`.

## Heavy masses and corrections

The 13 heavy entries in `so10_10_126_thresholds.json` are the remaining
10+126 components. The prototype adds the uneaten real 45 scalars:
`(8,1,1,0)`, `(1,3,1,0)`, `(1,1,3,0)`; its real singlet has no gauge effect.
The 45 has 45 real components: 30 Goldstones plus 8+3+3+1 uneaten fields.
Those 30 Goldstones belong to the broken SO(10)/3221 vector multiplets;
they are not also counted as physical heavy scalars.

The central hypothesis sets all these heavy scalars and broken GUT vectors
at MG. GUT scalar variations use `z_H=log10(MH/MG)` and the explicit box
`|z_H|≤log10 f`, primarily **f=3**, with **f=10** as a relaxation. A separate
family ties every scalar z to one common value. These are stipulated bounds
and correlations, not measured uncertainties or a probability distribution.
The vectors remain fixed at MG. All intermediate-scale heavy particles are
assumed degenerate at MI for leading-log matching. No independent MI
threshold knobs are introduced. The accepted prototype also requires every
heavy scalar to remain above MI and `MZ<MI<MG`.

For heavy scalar H,
`Δ_i = (b_i^H/2π) ln(MH/MG)` and `A_full = A_EFT + Δ`.
A direct step-decoupling test checks the sign by integrating with H active
above MH. The correction vanishes at degeneracy. This avoids fitting an
unidentified threshold vector. However, the masses themselves still encode
physical freedom; the independent box is deliberately more permissive than
most scalar potentials, so it can support conservative conditional bounds.

The optional finite-constant comparison implements BDM Eq. (26),
`A_p−C_p/(12π)=A_G−C_G/(12π)`: for SO(10), `C_G=8`, so
`Δ_G=(5,6,6,8)/(12π)` in LR order. At MI it adds `−1/(10π)` to the
hypercharge inverse coupling on the SM side, with L and colour unchanged.
This is a matching-convention sensitivity check; it does not supply omitted
mass-dependent intermediate thresholds or Yukawa effects.

## Solver and acceptance

At one loop, project LR b onto SM as `(4/5,−3,−7)` and solve a 2×2
linear system for the two logarithmic intervals, then the common inverse
coupling. At two loops, integrate **downward** from the common full-theory
coupling U: start at `U−Δ_G`, run four inverse couplings to MI, apply the
projection and matching constants, then run three inverse couplings to MZ.
Solve for `(log10 MI, log10 MG,U)` against the three input couplings.

The solved R/BL split is an output, not an extra low-energy observable or a
fourth independent success. Numerical acceptance requires optimizer success,
maximum input residual <1e-7, ordered scales, masses above MI and positive
perturbative trajectories `A_i>1`. This numerical domain is not a proof that
all scalar/Yukawa couplings are perturbative. Range: `log10 MI<15`,
`14.2≤log10 MG≤17.9`, `5≤U≤100`; a boundary minimum with nonzero residual
is reported rejected. Rejection is not global impossibility.

The three couplings at MZ are all used. There is no held-out gauge datum.
Proton decay and neutrino masses would be useful independent observables,
but this prototype lacks the operator/flavour/VEV calculations to predict
them. Conditional illustrations are reported separately from constraints.

## Explicit alternative

The parity-symmetric diagnostic adds ΔL and enforces L/R equality through
its symmetric spectrum. It is a comparison EFT, not a second vacuum of the
chosen D-parity-breaking 45 ansatz. Its equality is preserved by the RG flow
and therefore does not add a new independent equation. The low GUT scale
raises proton-decay concerns, but the toy lifetime does not establish exclusion.

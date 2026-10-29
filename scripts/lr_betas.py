"""
lr_betas.py -- derive gauge beta-function coefficients from a multiplet spectrum.

This module replaces the old practice of passing the left-right (LR) stage beta
coefficients (b3_LR, b2L_LR, b2R_LR, bBL_LR) on the command line.  Betas are
now *outputs* of a spectrum file that lists representations and quantum
numbers under

    SU(3)_c x SU(2)_L x SU(2)_R x U(1)_{B-L}      ("LR" group)

or, for validation only,

    SU(3)_c x SU(2)_L x U(1)_Y                    ("SM" group).

Conventions
-----------
Couplings run as  (t = ln mu)

    d alpha_i / dt = (b_i / 2pi) alpha_i^2 + (1 / 8pi^2) alpha_i^2 sum_j b_ij alpha_j

which is the form used by two_loop_core.step_rge_gauge_two_loop.  The general
gauge-sector coefficients for a product group with factors G_i are
(Machacek & Vaughn 1983; Jones 1982), per field of representation R with
Dynkin index T_i(R) in factor i, quadratic Casimir C2_i(R), and total
dimension d(R) = prod_k dim_k(R):

    S_i(R)   = T_i(R) * d(R) / dim_i(R)                      ("index with multiplicity")

    b_i      = -(11/3) C2(G_i)
             + (4/3) kappa * sum_F S_i(F)
             + (1/3) eta   * sum_S S_i(S)

    b_ij     = -(34/3) C2(G_i)^2 delta_ij
             + kappa * sum_F S_i(F) [ (20/3) C2(G_i) delta_ij + 4 C2_j(F) ]
             + eta   * sum_S S_i(S) [ ( 2/3) C2(G_i) delta_ij + 4 C2_j(S) ]

with kappa = 1/2 (Weyl fermion) or 1 (Dirac), eta = 1/2 (real scalar) or 1
(complex scalar).  For an abelian factor with normalised charge q the index and
Casimir are both q^2 and C2(G) = 0.

Abelian normalisations (GUT / SO(10) compatible, so that the matching
1/alpha_1 = (3/5)/alpha_2R + (2/5)/alpha_BL used in rge_core holds):

    U(1)_Y   :  q = sqrt(3/5) * Y            (Y = T3 + ... with Q = T3 + Y)
    U(1)_B-L :  q = sqrt(3/8) * (B-L)        (quarks B-L = 1/3, leptons -1)

The *only* two-loop pieces NOT included here are the Yukawa-coupling
contributions to the gauge betas.  There is no free parameter anywhere in this
module: everything follows from the representation content.

Validation: test_lr_betas.py checks that this routine reproduces the standard
one- and two-loop SM coefficients exactly (as rational numbers).

Spectrum file format ("lr_spectrum_v1")
---------------------------------------
{
  "format": "lr_spectrum_v1",
  "name": "...",
  "split_at_MI": "unification" | "lr_parity",   # how 1/alpha_Y splits into
                       # 1/alpha_2R and 1/alpha_BL at M_I: solved from the GUT
                       # boundary condition alpha_2R = alpha_BL at M_GUT, or
                       # fixed by manifest left-right symmetry g_L = g_R at M_I.
                       # A discrete model property, not a scanned parameter.
  "lr_stage": [        # fields assumed massless between M_I and M_GUT
    {"name": "Q_L", "rep": [3, 2, 1, "1/3"], "type": "weyl", "copies": 3}, ...
  ],
  "gut_thresholds": [  # heavy fields near M_GUT; give one-loop threshold logs
    {"name": "...", "rep": [3, 1, 1, "-2/3"], "type": "complex", "copies": 1,
     "mass_GeV": 1e16}, ...
  ]
}

"rep" lists the dimension under each non-abelian factor (sign/"bar" is
irrelevant for beta functions) followed by the un-normalised abelian charge
(B-L for the LR group).  Numbers may be given as JSON numbers or as strings
such as "1/3" to keep exact rational arithmetic.

"type" is one of: weyl, dirac, complex, real, vector.  "vector" is only
meaningful in "gut_thresholds" (a heavy gauge boson together with its eaten
Goldstone: b = -(21/6) T) -- the gauge bosons of the LR group itself are added
automatically to the LR-stage betas and must not be listed.
"""

from fractions import Fraction as Fr
import json
import math

import numpy as np

# --------------------------------------------------------------------------
# Group definitions
# --------------------------------------------------------------------------

# Each factor: ("SU", N) or ("U1", k2) where k2 is the squared normalisation
# applied to the raw charge:  q = sqrt(k2) * charge.
GROUPS = {
    "SM": {
        "factors": [("SU", 3), ("SU", 2), ("U1", Fr(3, 5))],
        "names": ["SU3c", "SU2L", "U1Y"],
        "charge_label": "Y",
    },
    "LR": {
        "factors": [("SU", 3), ("SU", 2), ("SU", 2), ("U1", Fr(3, 8))],
        "names": ["SU3c", "SU2L", "SU2R", "U1BL"],
        "charge_label": "B-L",
    },
}

# SU(3) irreps by dimension: (Dynkin index T, quadratic Casimir C2)
_SU3_TABLE = {
    1: (Fr(0), Fr(0)),
    3: (Fr(1, 2), Fr(4, 3)),
    6: (Fr(5, 2), Fr(10, 3)),
    8: (Fr(3), Fr(3)),
    10: (Fr(15, 2), Fr(6)),
    15: (Fr(10), Fr(16, 3)),
    27: (Fr(27), Fr(8)),
}

VALID_TYPES = ("weyl", "dirac", "complex", "real", "vector")


def su_n_invariants(N, dim):
    """Return (T, C2) for the SU(N) irrep of the given dimension.

    SU(2): any dimension n (spin j = (n-1)/2).
    SU(3): tabulated small irreps.
    Other N: singlet, fundamental and adjoint only.
    """
    dim = abs(int(dim))
    if dim == 1:
        return Fr(0), Fr(0)
    if N == 2:
        n = Fr(dim)
        return n * (n * n - 1) / 12, (n * n - 1) / 4
    if N == 3:
        if dim not in _SU3_TABLE:
            raise ValueError(f"SU(3) irrep of dimension {dim} not tabulated")
        return _SU3_TABLE[dim]
    if dim == N:
        return Fr(1, 2), Fr(N * N - 1, 2 * N)
    if dim == N * N - 1:
        return Fr(N), Fr(N)
    raise ValueError(f"SU({N}) irrep of dimension {dim} not supported")


def _as_fraction(x):
    if isinstance(x, Fr):
        return x
    if isinstance(x, str):
        return Fr(x.strip())
    if isinstance(x, int):
        return Fr(x)
    if isinstance(x, float):
        # keep exact arithmetic where possible (e.g. 0.5, 2.0)
        return Fr(x).limit_denominator(10**6)
    raise TypeError(f"cannot interpret {x!r} as a rational number")


def multiplet_invariants(group_name, rep):
    """Per-factor (T_i, C2_i, dim_i) for a representation list `rep`."""
    g = GROUPS[group_name]
    if len(rep) != len(g["factors"]):
        raise ValueError(
            f"group {group_name} needs {len(g['factors'])} rep entries, got {len(rep)}"
        )
    out = []
    for (kind, par), r in zip(g["factors"], rep):
        if kind == "SU":
            T, C2 = su_n_invariants(par, r)
            out.append((T, C2, Fr(abs(int(r)))))
        else:  # U1
            q2 = par * _as_fraction(r) ** 2  # normalised charge squared
            out.append((q2, q2, Fr(1)))
    return out


def gauge_casimirs(group_name):
    """C2(G_i) for each factor (0 for abelian)."""
    return [Fr(par) if kind == "SU" else Fr(0) for kind, par in GROUPS[group_name]["factors"]]


def field_one_loop(group_name, field):
    """One-loop b_i contribution (list of Fractions) of a single spectrum entry.

    Includes `copies`.  For type "vector" the entry is treated as a heavy
    massive gauge multiplet (gauge boson + eaten Goldstone), b = -(21/6) T,
    which is only meaningful for threshold bookkeeping.
    """
    ftype = field.get("type", "weyl")
    if ftype not in VALID_TYPES:
        raise ValueError(f"field {field.get('name')!r}: unknown type {ftype!r}")
    inv = multiplet_invariants(group_name, field["rep"])
    copies = _as_fraction(field.get("copies", 1))
    d_total = Fr(1)
    for _, _, d in inv:
        d_total *= d
    b = []
    for T, _, d in inv:
        S = T * d_total / d
        if ftype == "weyl":
            b.append(Fr(4, 3) * Fr(1, 2) * S)
        elif ftype == "dirac":
            b.append(Fr(4, 3) * S)
        elif ftype == "complex":
            b.append(Fr(1, 3) * S)
        elif ftype == "real":
            b.append(Fr(1, 3) * Fr(1, 2) * S)
        else:  # vector (massive, with Goldstone)
            b.append(-Fr(21, 6) * S)
    return [copies * x for x in b]


def field_two_loop(group_name, field):
    """Two-loop b_ij contribution (nested list of Fractions) of one entry.

    Vectors are rejected: heavy gauge multiplets are not part of the LR-stage
    running and their two-loop contribution is not needed for thresholds.
    """
    ftype = field.get("type", "weyl")
    if ftype == "vector":
        raise ValueError("two-loop contribution of a heavy vector multiplet is not supported")
    inv = multiplet_invariants(group_name, field["rep"])
    CG = gauge_casimirs(group_name)
    n = len(inv)
    copies = _as_fraction(field.get("copies", 1))
    d_total = Fr(1)
    for _, _, d in inv:
        d_total *= d
    if ftype in ("weyl", "dirac"):
        w = Fr(1, 2) if ftype == "weyl" else Fr(1)
        diag_coef = Fr(20, 3)
    else:
        w = Fr(1) if ftype == "complex" else Fr(1, 2)
        diag_coef = Fr(2, 3)
    bij = [[Fr(0)] * n for _ in range(n)]
    for i in range(n):
        T_i, _, d_i = inv[i]
        S_i = T_i * d_total / d_i
        for j in range(n):
            C2_j = inv[j][1]
            val = w * S_i * (4 * C2_j + (diag_coef * CG[i] if i == j else 0))
            bij[i][j] = copies * val
    return bij


def compute_betas(group_name, fields, include_gauge=True):
    """Total one- and two-loop coefficients for `fields` in the given group.

    Returns (b, bij) as lists of Fractions.  The gauge-boson self-contribution
    of each non-abelian factor is included unless include_gauge=False.
    """
    CG = gauge_casimirs(group_name)
    n = len(CG)
    b = [Fr(0)] * n
    bij = [[Fr(0)] * n for _ in range(n)]
    if include_gauge:
        for i in range(n):
            b[i] += -Fr(11, 3) * CG[i]
            bij[i][i] += -Fr(34, 3) * CG[i] ** 2
    for f in fields:
        if f.get("type") == "vector":
            raise ValueError(
                f"field {f.get('name')!r}: type 'vector' is only allowed in gut_thresholds"
            )
        b1 = field_one_loop(group_name, f)
        b2 = field_two_loop(group_name, f)
        for i in range(n):
            b[i] += b1[i]
            for j in range(n):
                bij[i][j] += b2[i][j]
    return b, bij


def to_float(b, bij):
    return (np.array([float(x) for x in b]),
            np.array([[float(x) for x in row] for row in bij]))


# --------------------------------------------------------------------------
# Spectrum files
# --------------------------------------------------------------------------

def load_spectrum(path):
    """Load and validate an lr_spectrum_v1 file.

    Raises a clear error for the retired raw-beta format (entries carrying a
    'b' vector instead of a representation), which is no longer accepted.
    """
    with open(path, "r", encoding="utf-8") as f:
        spec = json.load(f)
    if isinstance(spec, list) or spec.get("format") != "lr_spectrum_v1":
        raise ValueError(
            f"{path}: not an lr_spectrum_v1 file.  Raw beta vectors ('b': [...]) are "
            "no longer accepted; list representations and quantum numbers instead "
            "(see lr_betas.py docstring)."
        )
    spec.setdefault("group", "LR")
    if spec["group"] not in GROUPS:
        raise ValueError(f"unknown group {spec['group']!r}")
    spec.setdefault("lr_stage", [])
    spec.setdefault("gut_thresholds", [])
    spec.setdefault("split_at_MI", "unification")
    if spec["split_at_MI"] not in ("unification", "lr_parity"):
        raise ValueError(f"{path}: split_at_MI must be 'unification' or 'lr_parity'")
    for sec in ("lr_stage", "gut_thresholds"):
        for fld in spec[sec]:
            if "b" in fld:
                raise ValueError(
                    f"{path}: entry {fld.get('name')!r} carries a raw 'b' vector; "
                    "this knob has been removed -- specify 'rep' and 'type'."
                )
            if "rep" not in fld or "type" not in fld:
                raise ValueError(f"{path}: entry {fld.get('name')!r} needs 'rep' and 'type'")
            if fld["type"] not in VALID_TYPES:
                raise ValueError(f"{path}: entry {fld.get('name')!r} has unknown type {fld['type']!r}")
            multiplet_invariants(spec["group"], fld["rep"])  # validates rep
            if sec == "gut_thresholds":
                m = fld.get("mass_GeV")
                if m is None or m <= 0:
                    raise ValueError(f"{path}: threshold entry {fld.get('name')!r} needs mass_GeV > 0")
            elif fld["type"] == "vector":
                raise ValueError(f"{path}: entry {fld.get('name')!r}: vectors are not allowed in lr_stage")
    return spec


def lr_stage_betas(spec):
    """(b, bij) for the LR running stage, ordered (SU3, SU2L, SU2R, U1BL)."""
    return compute_betas(spec["group"], spec["lr_stage"], include_gauge=True)


def threshold_multiplets(spec):
    """List of {name, b (float array), mass_GeV} for GUT-scale heavy fields.

    b is the one-loop coefficient the multiplet *would* contribute to the
    LR-stage running; it is used only inside the leading-log threshold
    formula of rge_scan_v5.compute_multiplet_thresholds.
    """
    out = []
    for fld in spec["gut_thresholds"]:
        b = field_one_loop(spec["group"], fld)
        out.append({
            "name": fld.get("name", "?"),
            "b": np.array([float(x) for x in b]),
            "b_exact": b,
            "mass_GeV": float(fld["mass_GeV"]),
        })
    return out


def format_betas(group_name, b, bij):
    names = GROUPS[group_name]["names"]
    lines = ["  one-loop b_i:"]
    for nm, x in zip(names, b):
        lines.append(f"    {nm:5s}: {str(x):>8s}  ({float(x):+.6f})")
    lines.append("  two-loop b_ij (rows i, columns j):")
    lines.append("           " + "".join(f"{nm:>12s}" for nm in names))
    for nm, row in zip(names, bij):
        lines.append(f"    {nm:5s}  " + "".join(f"{str(x):>12s}" for x in row))
    return "\n".join(lines)


def describe_spectrum(spec):
    """Human-readable summary of the spectrum and the betas derived from it."""
    g = spec["group"]
    names = GROUPS[g]["names"]
    lines = [f"Spectrum: {spec.get('name', '?')}  (group {g})"]
    if spec.get("description"):
        lines.append(f"  {spec['description']}")
    lines.append("  LR-stage fields (massless between M_I and M_GUT):")
    for fld in spec["lr_stage"]:
        b = field_one_loop(g, fld)
        bs = ", ".join(f"{nm}={str(x)}" for nm, x in zip(names, b))
        lines.append(f"    {fld.get('name','?'):14s} rep={fld['rep']} type={fld['type']} "
                     f"copies={fld.get('copies',1)}  -> b: {bs}")
    b, bij = lr_stage_betas(spec)
    lines.append("  Derived LR-stage coefficients (gauge-only two-loop, no Yukawas):")
    lines.append(format_betas(g, b, bij))
    if spec["gut_thresholds"]:
        lines.append("  GUT-scale threshold fields (leading-log, one-loop):")
        for fld in spec["gut_thresholds"]:
            b = field_one_loop(g, fld)
            bs = ", ".join(f"{nm}={str(x)}" for nm, x in zip(names, b))
            lines.append(f"    {fld.get('name','?'):14s} rep={fld['rep']} type={fld['type']} "
                         f"copies={fld.get('copies',1)} M={fld['mass_GeV']:.3e} GeV -> b^H: {bs}")
    else:
        lines.append("  GUT-scale threshold fields: none (all heavy states assumed degenerate with M_GUT)")
    return "\n".join(lines)


# Standard-Model content, used for validation of the generic routine.
SM_FIELDS = [
    {"name": "Q",   "rep": [3, 2, "1/6"],  "type": "weyl",    "copies": 3},
    {"name": "u^c", "rep": [3, 1, "-2/3"], "type": "weyl",    "copies": 3},
    {"name": "d^c", "rep": [3, 1, "1/3"],  "type": "weyl",    "copies": 3},
    {"name": "L",   "rep": [1, 2, "-1/2"], "type": "weyl",    "copies": 3},
    {"name": "e^c", "rep": [1, 1, "1"],    "type": "weyl",    "copies": 3},
    {"name": "H",   "rep": [1, 2, "1/2"],  "type": "complex", "copies": 1},
]


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        print(describe_spectrum(load_spectrum(sys.argv[1])))
    else:
        b, bij = compute_betas("SM", SM_FIELDS)
        print("Standard Model (n_g=3, one Higgs doublet), GUT-normalised U(1)_Y:")
        print(format_betas("SM", b, bij))

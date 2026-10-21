"""
Consistency checks for lr_betas.py.  Run:  python test_lr_betas.py

1. The generic routine must reproduce the textbook SM one- and two-loop
   gauge coefficients (GUT-normalised U(1)_Y) exactly as rationals, and
   match the literal values hard-coded in rge_scan_v5.sm_two_loop_defaults.
2. The minimal LR content must give the textbook one-loop coefficients
   (-7, -3, -7/3, 11/2) for fermions + one bidoublet + one Delta_R.
3. The threshold routine must have the leading-log sign/normalisation:
   a heavy multiplet below M_GUT lowers 1/alpha_i by (b_i^H/2pi) ln(M_GUT/M_H).
"""
from fractions import Fraction as Fr
import math
import numpy as np

import lr_betas as lb
from rge_scan_v5 import sm_two_loop_defaults, compute_multiplet_thresholds


def test_sm():
    b, bij = lb.compute_betas("SM", lb.SM_FIELDS)
    assert b == [Fr(-7), Fr(-19, 6), Fr(41, 10)], b
    expect = [
        [Fr(-26), Fr(9, 2), Fr(11, 10)],
        [Fr(12), Fr(35, 6), Fr(9, 10)],
        [Fr(44, 5), Fr(27, 10), Fr(199, 50)],
    ]
    assert bij == expect, bij
    # The scanner stores SM in (U1, SU2, SU3) order; compare after reordering.
    b_lit, bij_lit = sm_two_loop_defaults()
    perm = [2, 1, 0]
    bf, bijf = lb.to_float(b, bij)
    assert np.allclose(bf[perm], b_lit)
    assert np.allclose(bijf[np.ix_(perm, perm)], bij_lit)


def test_lr_minimal():
    spec = lb.load_spectrum("spectra/lr_minimal.json")
    b, bij = lb.lr_stage_betas(spec)
    assert b == [Fr(-7), Fr(-3), Fr(-7, 3), Fr(11, 2)], b
    # Fermion-only B-L coefficient is 4 (the old config default).
    ferm = [f for f in spec["lr_stage"] if f["type"] == "weyl"]
    bf, _ = lb.compute_betas("LR", ferm)
    assert bf == [Fr(-7), Fr(-3) - Fr(1, 3), Fr(-3) - Fr(1, 3), Fr(4)], bf
    # SU(2)_L and SU(2)_R two-loop rows must coincide for L<->R symmetric
    # fermion content when scalar content is swapped accordingly.
    assert bij[0][1] == bij[0][2]  # SU3 row: colour couples equally to L and R


def test_threshold_sign():
    MGUT = 1e16
    heavy = [{"name": "x", "b": np.array([1.0, 0.0, 0.0, 0.0]), "mass_GeV": 1e15}]
    D = compute_multiplet_thresholds(heavy, MGUT)
    # colour: Delta_3 = (1/2pi) ln(1e15/1e16) < 0
    assert D[2] < 0 and math.isclose(D[2], math.log(0.1) / (2 * math.pi))
    assert D[0] == 0.0 and D[1] == 0.0
    # projection onto GUT-normalised U(1): 3/5 from SU(2)_R, 2/5 from B-L
    heavy = [{"name": "y", "b": np.array([0.0, 0.0, 1.0, 1.0]), "mass_GeV": 1e15}]
    D = compute_multiplet_thresholds(heavy, MGUT)
    assert math.isclose(D[0], math.log(0.1) / (2 * math.pi))
    # degenerate with M_GUT -> no correction
    heavy = [{"name": "z", "b": np.array([1.0, 1.0, 1.0, 1.0]), "mass_GeV": MGUT}]
    assert np.allclose(compute_multiplet_thresholds(heavy, MGUT), 0.0)


if __name__ == "__main__":
    test_sm()
    test_lr_minimal()
    test_threshold_sign()
    print("all lr_betas checks passed")

"""Independent oracles and failure cases; no fitted residual acceptance tests."""
from copy import deepcopy
from fractions import Fraction as F
import math

import numpy as np
import pytest

import model as m
from rge_core import find_MGUT_by_alpha2_eq_alpha3
from rge_scan_v5 import Model, compute_multiplet_thresholds


def test_published_lr_two_loop_matrix():
    # Bertolini/Di Luzio/Malinsky 0903.4049v2, Table VI, permuted to 3,L,R,BL.
    b, B = m.lb.lr_stage_betas(m.spectrum())
    assert b == [F(-7), F(-3), F(-7,3), F(11,2)]
    assert B == [[F(-26),F(9,2),F(9,2),F(1,2)],
                 [F(12),F(8),F(3),F(3,2)],
                 [F(12),F(3),F(80,3),F(27,2)],
                 [F(4),F(9,2),F(81,2),F(61,2)]]


def test_complete_scalar_representations_have_universal_index():
    spec = m.spectrum(include_45=False)
    scalars = [f for f in spec["lr_stage"] if f["type"] == "complex"] + spec["gut_thresholds"]
    assert sum(math.prod(f["rep"][:3]) for f in scalars) == 136
    contributions = [m.lb.field_one_loop("LR",f) for f in scalars]
    assert [sum(v[i] for v in contributions) for i in range(4)] == [F(12)]*4


def test_threshold_sign_from_explicit_decoupling_not_literal_formula():
    # Run low spectrum up to MH, add heavy b above MH, compare to light-only.
    a, b, bh = 40., -3., 2.
    mi, mh, mg = 1e10, 1e15, 1e16
    step = a-b*math.log(mh/mi)/m.TWOPI-(b+bh)*math.log(mg/mh)/m.TWOPI
    light = a-b*math.log(mg/mi)/m.TWOPI
    delta = compute_multiplet_thresholds([dict(b=np.array([0,bh,0,0]),mass_GeV=mh)],mg,False)
    assert light+delta[1] == pytest.approx(step,abs=1e-12)


def test_invariant_and_closed_form_solution():
    cfg, spec = m.inputs(), m.spectrum()
    a = m.low_inverse(cfg)
    result = m.analytic(cfg,spec)
    assert 2*m.W@a == pytest.approx(145*result["x"],abs=1e-10)
    assert result["determinant"] == pytest.approx(14.5)
    assert result["log10_MI"] == pytest.approx(9.0158165511,abs=1e-9)
    numerical = m.solve(cfg,spec,two_loop=False)
    assert numerical["accepted"]
    assert numerical["log10_MI"] == pytest.approx(result["log10_MI"],abs=1e-8)


def test_exact_box_support_and_attainment():
    cfg,spec=m.inputs(),m.spectrum()
    box=m.mass_box_bounds(cfg,spec,math.log10(3))
    assert sum(abs(x) for x in box["response"][0]) == pytest.approx(5.4)
    # Formal theorem's coefficients: exact representation arithmetic,
    # independent of floating-point threshold and matrix-solve implementations.
    coeff=[]
    for field in spec["gut_thresholds"]:
        b=m.lb.field_one_loop("LR",field)
        coeff.append(2*(20*(F(3,5)*b[2]+F(2,5)*b[3])-39*b[1]+19*b[0]))
    assert coeff == [9,9,9,9,207,75,-201,-129,-28,8,-18,14,14,19,-26,8]
    assert sum(abs(c) for c in coeff) == 783
    corner=m.analytic(cfg,spec,box["max_MI_corner"])
    assert corner["log10_MI"] == pytest.approx(box["upper"][0],abs=1e-12)
    assert box["upper"][0] < 12 < 13
    rng=np.random.default_rng(17)
    for _ in range(20):
        z=rng.uniform(-math.log10(3),math.log10(3),len(spec["gut_thresholds"]))
        li=m.analytic(cfg,spec,z)["log10_MI"]
        assert box["lower"][0] <= li <= box["upper"][0]


def test_parity_reduces_rank_in_symmetric_spectrum():
    spec=m.spectrum(parity=True,include_45=False)
    b,B=m.coefficients(spec)[2:]
    swap=[0,2,1,3]
    assert np.array_equal(b,b[swap])
    assert np.array_equal(B,B[np.ix_(swap,swap)])
    result=m.solve(m.inputs(),spec)
    assert result["accepted"]
    assert result["alpha_inv_LR_at_MI"][1] == pytest.approx(result["alpha_inv_LR_at_MI"][2],abs=1e-9)


def test_adaptive_vs_analytic_and_nonperturbative_rejection():
    b=np.array([-7.,-3.,.8])
    zero=np.zeros((3,3)); a=np.array([30.,40.,50.])
    expected=m.evolve(a,b,zero,15,two_loop=False)
    assert np.allclose(m.evolve(a,b,zero,15),expected,rtol=0,atol=1e-10)
    with pytest.raises(ValueError,match="nonperturbative"):
        m.evolve(np.ones(3)*2,np.ones(3)*100,zero,15)


def test_boundary_solution_not_accepted_as_unification():
    cfg=m.inputs(); cfg["alpha3_MZ"]=.01
    result=m.solve(cfg,m.spectrum())
    assert not result["accepted"]
    assert result["max_abs_input_residual"] > 1e-7


def test_locator_does_not_return_clamped_boundary_or_overflow():
    assert find_MGUT_by_alpha2_eq_alpha3(30,10,-3,-7,1e10) == (None,None)
    assert find_MGUT_by_alpha2_eq_alpha3(1e6,10,-3,-7,1e10) == (None,None)
    assert find_MGUT_by_alpha2_eq_alpha3(30,30,-3,-3,1e10) == (None,None)
    mu, lg = find_MGUT_by_alpha2_eq_alpha3(20,10,-3,-7,1e10)
    assert 10-(4/m.TWOPI)*math.log(mu/1e10) == pytest.approx(0,abs=1e-12)


def test_current_locator_rejects_absent_root():
    cfg=m.inputs(); cfg["root_bracket_log10"]=[14.2,14.3]
    model=Model(cfg,m.spectrum(),rk_steps=100)
    assert model.evaluate(9.45) is None


def test_current_root_agrees_with_independent_boundary_solver():
    cfg,spec=m.inputs(),m.spectrum()
    # Remove nondegenerate absolute-mass thresholds: compare the same EFT.
    spec["gut_thresholds"]=[]
    independent=m.solve(cfg,spec)
    result=Model(cfg,spec,rk_steps=200).evaluate(independent["log10_MI"])
    assert result["log10_MGUT"] == pytest.approx(independent["log10_MGUT"],abs=1e-8)
    assert result["max_frac"] < 1e-8


def test_reference_checker_detects_drift(tmp_path):
    import json
    from check_references import check, REFERENCES
    for name,(mi,mg,a) in REFERENCES.items():
        target=tmp_path/name/"best_fit.json"
        target.parent.mkdir()
        target.write_text(json.dumps(dict(best=dict(log10_MI=mi,log10_MGUT=mg,alpha_GUT=a,max_frac=0))))
    check(tmp_path)
    target=tmp_path/"lr_minimal"/"best_fit.json"
    d=json.loads(target.read_text()); d["best"]["log10_MI"]+=.1
    target.write_text(json.dumps(d))
    with pytest.raises(ValueError,match="differs"):
        check(tmp_path)

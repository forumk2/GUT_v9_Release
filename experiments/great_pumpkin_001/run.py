"""Deterministic evidence generator. Run from any working directory."""
from copy import deepcopy
from itertools import product
from pathlib import Path
import hashlib
import json
import math
import platform
import sys

import numpy as np
import scipy
from scipy.optimize import brentq, linprog, minimize

import model as m
from rge_scan_v5 import Model
from two_loop_core import run_segment_two_loop


def numeric_checks(cfg, spec, central):
    logI = central["log10_MI"]
    plain = m.lb.load_spectrum(m.ROOT / "scripts/spectra/lr_minimal.json")
    checks = []
    for n in [100, 200, 400, 800, 1600]:
        row = Model(cfg, plain, rk_steps=n).evaluate(logI)
        checks.append(dict(steps=n, log10_MGUT=row["log10_MGUT"],
                           max_frac=row["max_frac"],
                           four_inverse=row["alpha_inv_at_MGUT"]))
    # An independently bracketed forward root; no scanner scale locator.
    bs, Bs, bl, Bl = m.coefficients(spec)
    atI_sm = m.evolve(m.low_inverse(cfg), bs, Bs,
                      (logI-math.log10(cfg["MZ_GeV"]))*m.LN10)
    def at(logG):
        def lr(r):
            A = atI_sm[0]
            y = A / (.6*r+.4)
            inv = np.array([atI_sm[2], atI_sm[1], r*y, y])
            return m.evolve(inv, bl, Bl, (logG-logI)*m.LN10)
        r = brentq(lambda r: (lr(r)[2]-lr(r)[3]), .3, 1.5, xtol=1e-13)
        return lr(r)
    root = brentq(lambda lg: at(lg)[1]-at(lg)[0], 15, 17, xtol=1e-12)
    methods = [m.solve(cfg, spec, method="RK45", rtol=1e-10),
               m.solve(cfg, spec, method="DOP853", rtol=1e-12)]
    starts = [m.solve(cfg, spec, start=s) for s in
              [[8, 16, 44], [11, 15.5, 48], [10, 17.5, 43]]]
    return dict(rk4=checks, forward_brent_log10_MGUT=root,
                adaptive=methods, alternative_starts=starts)


def input_sensitivity(spec):
    # PDG 2025 review revised November 2025, sections 10.3.1, 10.4.5, Table 10.6.
    # Separate one-at-a-time and box sensitivity, NOT a joint likelihood.
    cfg = m.inputs()
    cfg.update(alpha_em_MZ=1/127.930, sin2thetaW_MZ=.23122,
               alpha3_MZ=.1177, MZ_GeV=91.1880)
    central = m.solve(cfg, spec)
    rows = []
    for em, sw, strong in product([-1, 1], repeat=3):
        c = dict(cfg, alpha_em_MZ=1/(127.930+em*.008),
                 sin2thetaW_MZ=.23122+sw*.00006, alpha3_MZ=.1177+strong*.0009)
        rows.append(dict(signs=[em, sw, strong], result=m.solve(c, spec)))
    return dict(inputs=cfg, independent_halfwidths=[.008, .00006, .0009],
                interpretation="rectangular sensitivity envelope, no covariance or confidence level",
                central=central, corners=rows)


def thresholds_challenge(cfg, spec, count):
    rho = math.log10(3)
    box = m.mass_box_bounds(cfg, spec, rho)
    matrix = np.array(box["response"])
    rng = np.random.default_rng(20261028)
    n = len(spec["gut_thresholds"])
    zs = rng.uniform(-rho, rho, size=(count, n))
    rows = []
    for i, z in enumerate(zs):
        result = m.solve(cfg, spec, z)
        rows.append(dict(z=z.tolist(), result=result))
    corners = {key: m.solve(cfg, spec, box[key]) for key in ["max_MI_corner", "min_MI_corner"]}
    # Correlated common scalar mass: one free mass ratio, not 16 independent ones.
    common = [dict(z=z, one_loop=m.analytic(cfg, spec, np.full(n,z)),
                   two_loop=m.solve(cfg, spec, np.full(n,z))) for z in [-rho, 0, rho]]
    # Minimal mass splitting needed to attain MI>=1e13 at one loop, while
    # keeping MG in the documented bracket. The objective is rho, not a residual.
    center = np.array(box["center"])
    A = np.vstack([np.column_stack((np.eye(n), -np.ones(n))),
                   np.column_stack((-np.eye(n), -np.ones(n))),
                   np.r_[-matrix[0], 0], np.r_[-matrix[1], 0],
                   np.r_[matrix[1], 0], np.r_[matrix[0]-matrix[1], 0]])
    b = np.r_[np.zeros(2*n), center[0]-13, center[1]-14.2,
              17.9-center[1], center[1]-center[0]]
    lp = linprog(np.r_[np.zeros(n), 1], A_ub=A, b_ub=b,
                 bounds=[(None,None)]*n+[(0,None)], method="highs")
    if not lp.success:
        raise RuntimeError(lp.message)
    lp_z = lp.x[:-1]
    # Attempt to defeat the two-loop corner with bounded minimisation.
    # This is local numerical evidence, not a global two-loop theorem.
    def objective(z):
        result = m.solve(cfg, spec, z)
        return -result["log10_MI"] if result["accepted"] else 1e3
    opt = minimize(objective, box["max_MI_corner"], method="L-BFGS-B",
                   bounds=[(-rho,rho)]*n, options={"maxiter":25, "ftol":1e-10})
    return dict(factor3_exact_one_loop_box=box, random_seed=20261028,
                samples=rows, two_loop_corners=corners, common_mass=common,
                local_two_loop_max=dict(success=bool(opt.success), message=str(opt.message),
                                        z=opt.x.tolist(), result=m.solve(cfg,spec,opt.x)),
                one_loop_minimum_factor_for_MI_1e13=dict(factor=10**lp.x[-1],
                    log10_factor=lp.x[-1], z=lp_z.tolist(),
                    result=m.analytic(cfg,spec,lp_z), two_loop_at_same_masses=m.solve(cfg,spec,lp_z)))


def adversarial(cfg, spec):
    bs, _, bl, _ = m.coefficients(spec)
    A = m.low_inverse(cfg)
    rows = []
    for li, lg, U in [(5.,15.,40.), (10.,16.,45.), (13.,17.,50.)]:
        x = (li-math.log10(cfg["MZ_GeV"]))*m.LN10/m.TWOPI
        y = (lg-li)*m.LN10/m.TWOPI
        raw = A - bs*x - (m.P@bl)*y
        shifts = U-raw
        rows.append(dict(log10_MI=li, log10_MGUT=lg, alpha_inv_GUT=U,
                         free_inverse_shifts=shifts.tolist(),
                         max_abs_residual=float(np.max(np.abs(raw+shifts-U)))))
    return rows


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=64)
    args = ap.parse_args()
    cfg, spec = m.inputs(), m.spectrum()
    center = m.solve(cfg, spec)
    one = m.analytic(cfg, spec)
    print("central", center, flush=True)
    checks = numeric_checks(cfg, spec, center)
    print("numerics checked", flush=True)
    sensitivity = input_sensitivity(spec)
    print("inputs checked", flush=True)
    thresholds = thresholds_challenge(cfg, spec, args.samples)
    print("threshold challenge complete", flush=True)
    parity = m.spectrum(parity=True, include_45=False)
    # A D-parity-preserving diagnostic EFT, not the chosen 45 breaking vacuum.
    parity_result = m.solve(cfg, parity)
    # Seesaw numbers are conditional illustrations, not predictions.
    mi, mg, U = 10**center["log10_MI"], 10**center["log10_MGUT"], center["alpha_inv_GUT"]
    y_for_005 = math.sqrt(.05*mi/1e9)/174
    tau = 1e35*(mg/1e16)**4*(.04*U)**2
    report = dict(environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
        input_config=cfg, spectrum=spec,
        coefficients=dict(one_loop=[str(x) for x in m.lb.lr_stage_betas(spec)[0]],
                          two_loop=[[str(x) for x in row] for row in m.lb.lr_stage_betas(spec)[1]]),
        one_loop=one, two_loop=center,
        finite_matching_constants_only=m.solve(cfg,spec,finite=True),
        numerical_checks=checks, input_sensitivity=sensitivity, thresholds=thresholds,
        adversarial_free_offsets=adversarial(cfg,spec), parity_diagnostic=parity_result,
        conditional_observables=dict(tree_dimensional_tau_years=tau,
           yukawa_for_mnu_005eV_and_MR_eq_MI=y_for_005,
           mnu_eV_if_y_01_and_MR_eq_MI=(.1*174)**2/mi*1e9,
           proton_mass_factor_to_saturate_SK=(2.4e34/tau)**.25))
    out = Path(__file__).parent/"results"
    out.mkdir(exist_ok=True)
    report["source_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in [Path(__file__), Path(m.__file__)]}
    (out/"experiments.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("wrote",out/"experiments.json")


if __name__ == "__main__":
    main()

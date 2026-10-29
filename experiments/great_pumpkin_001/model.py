"""Conditional SO(10) -> 3221 -> SM gauge-sector experiment.

Independent inverse-coupling evolution / boundary solver. Representation
arithmetic is reused only after comparison to published coefficient tables.
No Yukawa, scalar-potential, or flavour fit is claimed. See model.md.
"""
from copy import deepcopy
from pathlib import Path
import json
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import lr_betas as lb

TWOPI = 2 * math.pi
LN10 = math.log(10)
P = np.array([[0, 0, .6, .4], [0, 1, 0, 0], [1, 0, 0, 0]])
DIFF = np.array([[1, -1, 0], [0, 1, -1]])
W = np.array([20., -39., 19.])


def inputs():
    return json.loads((ROOT / "scripts/config.json").read_text())


def spectrum(parity=False, include_45=True):
    """Fixed light fields; heavy uneaten 10+126+45 scalars.

    The 30 broken SO(10)/3221 generators and their Goldstones have mass MG.
    The 45 leaves a real colour octet and real L/R triplets, plus a singlet
    with no gauge effect. This is mass bookkeeping, not a vacuum solution.
    """
    spec = lb.load_spectrum(ROOT / "scripts/spectra/so10_10_126_thresholds.json")
    spec = deepcopy(spec)
    spec["name"] = "45_126_complex10_survival"
    if include_45:
        for name, rep in [("45_H: octet", [8, 1, 1, 0]),
                          ("45_H: left", [1, 3, 1, 0]),
                          ("45_H: right", [1, 1, 3, 0])]:
            spec["gut_thresholds"].append(dict(name=name, rep=rep, type="real", copies=1,
                                                mass_GeV=1e16))
    if parity:
        spec["lr_stage"].append(dict(name="Delta_L", rep=[1, 3, 1, 2], type="complex"))
        spec["gut_thresholds"] = [f for f in spec["gut_thresholds"]
                                   if f["name"] != "126_H: Delta_L"]
    return spec


def low_inverse(cfg):
    em, s, strong = cfg["alpha_em_MZ"], cfg["sin2thetaW_MZ"], cfg["alpha3_MZ"]
    return np.array([.6 * (1 - s) / em, s / em, 1 / strong])


def coefficients(spec):
    bs, Bs = lb.to_float(*lb.compute_betas("SM", lb.SM_FIELDS))
    perm = [2, 1, 0]
    bl, Bl = lb.to_float(*lb.lr_stage_betas(spec))
    return bs[perm], Bs[np.ix_(perm, perm)], bl, Bl


def threshold_matrix(spec):
    return np.array([list(map(float, lb.field_one_loop("LR", f)))
                     for f in spec["gut_thresholds"]]).T.reshape(4, -1)


def thresholds(spec, log_mass_ratios):
    """Delta(full - EFT); z=log10(M_H/M_G), bounded externally."""
    z = np.asarray(log_mass_ratios, dtype=float)
    if z.shape != (len(spec["gut_thresholds"]),) or not np.isfinite(z).all():
        raise ValueError("one finite log10 mass ratio per heavy field required")
    return threshold_matrix(spec) @ z * LN10 / TWOPI


def evolve(a, b, B, log_ratio, two_loop=True, rtol=1e-10, method="DOP853"):
    """Integrate A_i=1/alpha_i; alpha_i < 1 throughout accepted segments."""
    a = np.asarray(a, dtype=float)
    if np.min(a) <= 1 or not np.isfinite(a).all():
        raise ValueError("nonperturbative or nonfinite initial coupling")
    if not two_loop:
        result = a - b * log_ratio / TWOPI
        if np.min(result) <= 1:
            raise ValueError("nonperturbative one-loop trajectory")
        return result

    def rhs(t, inv):
        return -b / TWOPI - (B @ (1 / inv)) / (8 * math.pi**2)

    def boundary(t, inv):
        return np.min(inv) - 1
    boundary.terminal = True
    sol = solve_ivp(rhs, (0, log_ratio), a, method=method, rtol=rtol,
                    atol=rtol * .01, events=boundary)
    if not sol.success or len(sol.t_events[0]) or not np.isfinite(sol.y).all():
        raise ValueError("nonperturbative or failed integration")
    return sol.y[:, -1]


def matching_constants(finite=False):
    # BDM 2009 Eq. (26): A_p-C_p/(12pi)=A_G-C_G/(12pi).
    gut = np.array([5, 6, 6, 8]) / (12 * math.pi) if finite else np.zeros(4)
    intermediate = np.array([-1 / (10 * math.pi), 0, 0]) if finite else np.zeros(3)
    return gut, intermediate


def analytic(cfg, spec, z=None, finite=False):
    """Exact affine solution in x=ln(MI/MZ)/(2pi), y=ln(MG/MI)/(2pi)."""
    if z is None:
        z = np.zeros(len(spec["gut_thresholds"]))
    bs, _, bl, _ = coefficients(spec)
    dg, di = matching_constants(finite)
    delta = P @ (thresholds(spec, z) + dg) - di
    C = np.column_stack((DIFF @ bs, DIFF @ (P @ bl)))
    if abs(np.linalg.det(C)) < 1e-12:
        raise ValueError("singular scale system: intersection not identifiable")
    x, y = np.linalg.solve(C, DIFF @ (low_inverse(cfg) + delta))
    U = (low_inverse(cfg) + delta - bs * x - (P @ bl) * y).mean()
    logI = math.log10(cfg["MZ_GeV"]) + TWOPI / LN10 * x
    logG = logI + TWOPI / LN10 * y
    return dict(log10_MI=logI, log10_MGUT=logG, alpha_inv_GUT=U,
                x=x, y=y, determinant=float(np.linalg.det(C)),
                ordered=bool(x > 0 and y > 0 and U > 1))


def solve(cfg, spec, z=None, two_loop=True, finite=False, rtol=1e-10,
          method="DOP853", start=None):
    """Shoot downward from a common GUT boundary, fit three measured inputs.

    Scales and U are outputs conditional on the spectrum. A small residual
    establishes existence only, not an independent low-energy prediction.
    """
    if z is None:
        z = np.zeros(len(spec["gut_thresholds"]))
    dg, di = matching_constants(finite)
    delta = thresholds(spec, z) + dg
    bs, Bs, bl, Bl = coefficients(spec)
    mz = math.log10(cfg["MZ_GeV"])

    def predict(p, detail=False):
        li, lg, U = p
        at_I = evolve(U - delta, bl, Bl, (li-lg)*LN10, two_loop, rtol, method)
        at_Z = evolve(P @ at_I + di, bs, Bs, (mz-li)*LN10, two_loop, rtol, method)
        return (at_Z, at_I) if detail else at_Z

    def residual(p):
        try:
            return predict(p) - low_inverse(cfg)
        except ValueError:
            return np.full(3, 1e6)

    if start is None:
        first = analytic(cfg, spec, z, finite)
        start = [first["log10_MI"], first["log10_MGUT"], first["alpha_inv_GUT"]]
    lower, upper = [mz+.001, 14.2, 5.], [15., 17.9, 100.]
    fit = least_squares(residual, np.clip(start, np.array(lower)+1e-6, np.array(upper)-1e-6),
                        bounds=(lower, upper), xtol=1e-12, ftol=1e-12, gtol=1e-12,
                        max_nfev=150, x_scale="jac")
    li, lg, U = fit.x
    err = float(np.max(np.abs(residual(fit.x))))
    accepted = bool(fit.success and err < 1e-7 and mz < li < lg and U > 1
                    and np.all(lg + np.asarray(z) > li))
    try:
        _, at_I = predict(fit.x, True)
    except ValueError:
        return dict(log10_MI=li, log10_MGUT=lg, alpha_inv_GUT=U,
                    max_abs_input_residual=err, accepted=False,
                    alpha_inv_LR_at_MI=None, split_r=None, nfev=fit.nfev,
                    method=method, two_loop=two_loop, finite_constants=finite)
    return dict(log10_MI=li, log10_MGUT=lg, alpha_inv_GUT=U,
                max_abs_input_residual=err, accepted=accepted,
                alpha_inv_LR_at_MI=at_I.tolist(), split_r=at_I[2]/at_I[3],
                nfev=fit.nfev, method=method, two_loop=two_loop, finite_constants=finite)


def mass_box_bounds(cfg, spec, decades):
    """Exact coordinate extrema over independent |log10(MH/MG)|<=decades.

    These extrema are attained at (generally different) corners of the box.
    They are NOT simultaneous extrema or confidence intervals.
    """
    bs, _, bl, _ = coefficients(spec)
    C = np.column_stack((DIFF @ bs, DIFF @ (P @ bl)))
    response_xy = np.linalg.solve(C, DIFF @ P @ threshold_matrix(spec))
    response = np.vstack([response_xy[0], response_xy.sum(axis=0)])
    base = analytic(cfg, spec)
    center = np.array([base["log10_MI"], base["log10_MGUT"]])
    halfwidth = decades * np.abs(response).sum(axis=1)
    return dict(center=center.tolist(), halfwidth=halfwidth.tolist(),
                lower=(center-halfwidth).tolist(), upper=(center+halfwidth).tolist(),
                response=response.tolist(),
                max_MI_corner=(decades*np.sign(response[0])).tolist(),
                min_MI_corner=(-decades*np.sign(response[0])).tolist())

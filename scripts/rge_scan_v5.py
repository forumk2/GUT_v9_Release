"""
rge_scan_v5.py -- two-loop gauge-coupling scan for an intermediate-scale
left-right model, in "GUT-scan mode": every input is either low-energy data
(config.json) or multiplet content (a spectrum file).  There are no tunable
slope coefficients and no manual threshold shifts.

Pipeline for each trial intermediate scale M_I on a grid:

  1. SM two-loop running M_Z -> M_I            (two_loop_core.run_segment_two_loop)
  2. tree-level matching SM -> LR at M_I        (rge_core.match_SM_to_LR_with_ratio)
  3. locate M_GUT from alpha_2L = alpha_3: a bracketed, safeguarded secant
     enforces the crossing with the two-loop, threshold-corrected couplings
                                               (Model.locate_MGUT)
  4. LR two-loop running M_I -> M_GUT with betas DERIVED from the spectrum
                                               (lr_betas.lr_stage_betas)
  5. one-loop leading-log thresholds at M_GUT from the spectrum's heavy
     multiplets and their masses              (compute_multiplet_thresholds)
  6. the split of U(1)_Y into SU(2)_R x U(1)_B-L at M_I is a model property
     declared in the spectrum file ("split_at_MI"), never scanned:
        "unification": solved from alpha_2R(M_GUT) = alpha_BL(M_GUT)
        "lr_parity":   g_L = g_R at M_I (manifest left-right symmetry)
                                               (Model.split_at_MI)
  7. residual metric over (alpha_1, alpha_2L, alpha_3) at M_GUT
                                               (rge_core.mismatch_metrics)
  8. tree-level proton-lifetime estimate       (proton_decay_tau_years)

The reported point is the grid minimum of the residual, refined by a local
golden-section search in log10(M_I).

HOW TO READ THE RESIDUAL.  In "unification" mode the unknowns (M_I, M_GUT,
alpha_GUT, split) are exactly as many as the four unification conditions, so
a near-zero residual means only that a self-consistent solution exists inside
the scanned range; the physical content of the result is the *predicted*
M_I, M_GUT, alpha_GUT and the implied alpha_2R(M_I).  In "lr_parity" mode the
split is fixed by symmetry. For symmetric spectra this condition is redundant
under running; only independent extra conditions over-determine the system.

WHAT REMAINS APPROXIMATE (see also README.md):
  * tau_p is a tree-level, dimensional estimate: no hadronic matrix element,
    no operator running, no flavour structure; M_X = M_GUT is assumed.
  * Two-loop LR coefficients are gauge-only; Yukawa contributions (top,
    bottom, tau, neutrino Dirac/Majorana) are not included in either stage.
  * Thresholds are leading-log one-loop step functions.  Scheme-dependent
    finite constants (e.g. the MS-bar gauge-boson term) are not included.
  * Matching at M_I is tree-level: no thresholds for W_R, Delta_R, nu_R or the
    heavy half of the bidoublet; one SM Higgs doublet is assumed light.
  * M_GUT requires a converged crossing within the configured bracket;
    a bracket boundary without a crossing is rejected.
  * The quantum numbers of every field are whatever the spectrum file says;
    see the spectrum 'description' for the assumptions made in the examples.
"""
import argparse
import csv
import json
import math
import os

import numpy as np

from rge_core import (gut_normalized_alphas_at_MZ, match_SM_to_LR_with_ratio,
                      mismatch_metrics, mismatch_metrics_n)
from two_loop_core import run_segment_two_loop, run_segment_one_loop
import lr_betas


def load_config(path):
    with open(path, "r") as f:
        return json.load(f)


def proton_decay_tau_years(MX_GeV, alpha_GUT, alpha_H_GeV3):
    """Tree-level dimensional estimate, tau ~ M_X^4 / (alpha_GUT^2 alpha_H^2).
    Order-of-magnitude placeholder only (see module docstring)."""
    return (1e35 * (MX_GeV / 1e16) ** 4 * (0.04 / alpha_GUT) ** 2 * (0.012 / alpha_H_GeV3) ** 2)


def seesaw_mnu_eV(y, v_GeV, MR_GeV):
    return (y * v_GeV) ** 2 / MR_GeV * 1e9


def sm_two_loop_defaults():
    """SM gauge coefficients, order (U1_Y GUT-normalised, SU2_L, SU3_c),
    n_g = 3, one Higgs doublet.  test_lr_betas.py checks that lr_betas derives
    exactly these numbers from the SM field content."""
    b = np.array([41.0 / 10.0, -19.0 / 6.0, -7.0])
    bij = np.array([
        [199.0 / 50.0, 27.0 / 10.0, 44.0 / 5.0],
        [9.0 / 10.0, 35.0 / 6.0, 12.0],
        [11.0 / 10.0, 9.0 / 2.0, -26.0],
    ])
    return b, bij


def compute_multiplet_thresholds(heavy, MGUT, project=True):
    """One-loop leading-log threshold corrections at M_GUT.

    heavy : list of {"name", "b": 4-vector (SU3, SU2L, SU2R, U1BL), "mass_GeV"}
            as produced by lr_betas.threshold_multiplets.  b is the one-loop
            coefficient the multiplet would contribute to the running; it is
            derived from the representation, never typed in.

    A heavy multiplet H of mass M_H changes the full-theory inverse coupling
    at the matching scale relative to running with the light content alone by

        Delta_i = (b_i^H / 2 pi) * ln(M_H / M_GUT)

    (valid for M_H on either side of M_GUT; it vanishes for M_H = M_GUT).

    Returns the 4-vector in LR order, or with project=True the 3-vector in
    the (alpha_1, alpha_2L, alpha_3) basis used by the residual metric, with
    Delta_1 = (3/5) Delta_2R + (2/5) Delta_BL (same projection as the
    coupling itself).

    NOTE: the retired version of this function used +b/(12 pi) ln(M_GUT/M_H)
    with hand-typed b vectors; the coefficient and sign have been corrected to
    the standard leading-log result and the inputs are now derived.
    """
    D = np.zeros(4)
    for H in heavy:
        D += (np.asarray(H["b"], dtype=float) / (2.0 * math.pi)) * math.log(H["mass_GeV"] / MGUT)
    if not project:
        return D
    return np.array([0.6 * D[2] + 0.4 * D[3], D[1], D[0]])


class Model:
    """Everything fixed by config + spectrum; evaluates one trial M_I."""

    def __init__(self, cfg, spec, two_loop=True, rk_steps=400):
        self.cfg = cfg
        self.spec = spec
        self.two_loop = two_loop
        self.rk_steps = rk_steps
        self.MZ = cfg["MZ_GeV"]
        a1, a2, a3 = gut_normalized_alphas_at_MZ(cfg["alpha_em_MZ"], cfg["sin2thetaW_MZ"], cfg["alpha3_MZ"])
        self.alphas_MZ = np.array([a1, a2, a3])
        self.bSM_1, self.bSM_2 = sm_two_loop_defaults()
        b, bij = lr_betas.lr_stage_betas(spec)
        self.bLR_exact, self.bijLR_exact = b, bij
        self.bLR_1, self.bLR_2 = lr_betas.to_float(b, bij)       # order (3, 2L, 2R, BL)
        self.heavy = lr_betas.threshold_multiplets(spec)
        self.bracket = tuple(cfg.get("root_bracket_log10", (14.2, 17.9)))
        self.split_mode = spec.get("split_at_MI", "unification")

    # -- stage 1: SM running ------------------------------------------------
    def sm_to_MI(self, MI):
        if self.two_loop:
            a = run_segment_two_loop(self.alphas_MZ, self.bSM_1, self.bSM_2, self.MZ, MI, nsteps=self.rk_steps)
            return 1.0 / a[0], 1.0 / a[1], 1.0 / a[2]
        a_inv = run_segment_one_loop(1.0 / self.alphas_MZ, self.bSM_1, self.MZ, MI)
        return tuple(a_inv.tolist())

    # -- stage 4: LR running for a given split ratio ----------------------------
    def lr_to_MGUT(self, a1inv_MI, a2inv_MI, a3inv_MI, r, MI, MGUT):
        a3inv, aLinv, aRinv, aBLinv = match_SM_to_LR_with_ratio(a1inv_MI, a2inv_MI, a3inv_MI, r_ratio=r)
        if aRinv <= 0 or aBLinv <= 0:
            return None
        if self.two_loop:
            alphas_MI = np.array([1.0 / a3inv, 1.0 / aLinv, 1.0 / aRinv, 1.0 / aBLinv])
            a = run_segment_two_loop(alphas_MI, self.bLR_1, self.bLR_2, MI, MGUT, nsteps=self.rk_steps)
            return 1.0 / a  # inverse couplings at M_GUT, LR order
        a_inv = np.array([a3inv, aLinv, aRinv, aBLinv])
        return run_segment_one_loop(a_inv, self.bLR_1, MI, MGUT)

    # -- stage 6: the U(1)_Y -> SU(2)_R x U(1)_B-L split at M_I -----------------------
    def split_at_MI(self, a1inv_MI, a2inv_MI, a3inv_MI, MI, MGUT, Delta4, tol=1e-10, max_iter=25):
        """Return (r, a_inv4_full): the split ratio r = (1/alpha_2R)/(1/alpha_BL) at
        M_I and the full-theory inverse couplings at M_GUT (thresholds applied).

        Two model-defined options (spectrum key "split_at_MI"):
          "unification" (default): r is solved so that alpha_2R(M_GUT) = alpha_BL(M_GUT),
              i.e. the GUT boundary condition.  Appropriate when left-right (D-)parity
              is broken at M_GUT so g_L != g_R below it.  The one-loop solution is
              analytic and seeds a secant iteration for the two-loop case.
          "lr_parity": manifest left-right symmetry, g_L = g_R at M_I, so
              1/alpha_2R(M_I) = 1/alpha_2L(M_I) and 1/alpha_BL follows from
              1/alpha_1.  Nothing is solved for; the residual is then a genuine
              test of the content only if it is an independent constraint.
        Returns None if no physical (positive-coupling) choice exists.
        """
        if self.split_mode == "lr_parity":
            x = a2inv_MI
            y = (a1inv_MI - 0.6 * x) / 0.4
            if x <= 0 or y <= 0:
                return None
            r = x / y
            a = self.lr_to_MGUT(a1inv_MI, a2inv_MI, a3inv_MI, r, MI, MGUT)
            return None if a is None else (r, a + Delta4)

        L = math.log(MGUT / MI)
        b2R, bBL = self.bLR_1[2], self.bLR_1[3]
        # one-loop: x - (b2R/2pi) L + Delta_R = y - (bBL/2pi) L + Delta_BL,  A = 0.6 x + 0.4 y
        D = (b2R - bBL) * L / (2.0 * math.pi) + (Delta4[3] - Delta4[2])
        y = a1inv_MI - 0.6 * D
        x = y + D
        if x <= 0 or y <= 0:
            return None
        r = x / y

        def f(rv):
            a = self.lr_to_MGUT(a1inv_MI, a2inv_MI, a3inv_MI, rv, MI, MGUT)
            if a is None:
                return None, None
            full = a + Delta4
            return full[2] - full[3], full

        f0, full0 = f(r)
        if f0 is None:
            return None
        if not self.two_loop or abs(f0) < tol:
            return r, full0
        r1 = r * (1.0 + 1e-3)
        f1, full1 = f(r1)
        if f1 is None:
            return None
        for _ in range(max_iter):
            if abs(f1) < tol or f1 == f0:
                break
            r2 = r1 - f1 * (r1 - r) / (f1 - f0)
            r, f0 = r1, f1
            r1 = r2
            f1, full1 = f(r1)
            if f1 is None:
                return None
        if not math.isfinite(f1) or abs(f1) >= tol:
            return None
        return r1, full1

    # -- stages 3-6 at fixed M_I ---------------------------------------------------
    def _at_MGUT(self, a1inv_MI, a2inv_MI, a3inv_MI, MI, logMGUT):
        MGUT = 10 ** logMGUT
        Delta4 = compute_multiplet_thresholds(self.heavy, MGUT, project=False)
        sol = self.split_at_MI(a1inv_MI, a2inv_MI, a3inv_MI, MI, MGUT, Delta4)
        if sol is None:
            return None
        r, full = sol
        return {"MGUT": MGUT, "logMGUT": logMGUT, "Delta4": Delta4, "r": r, "full": full,
                "gap23": full[1] - full[0]}  # 1/alpha_2L - 1/alpha_3 at M_GUT

    def locate_MGUT(self, a1inv_MI, a2inv_MI, a3inv_MI, MI, tol=1e-10, max_iter=25):
        """Safeguarded secant crossing, with a sign bracket and residual check.

        Unlike the pre-expedition implementation, an out-of-range analytic
        seed does not become an alleged crossing at a clamped boundary. A
        two-loop root may still exist when the one-loop seed is out of range.
        No claim is made to enumerate tangential or multiple crossings.
        """
        lo, hi = self.bracket
        lo = max(lo, math.log10(MI))
        if lo >= hi:
            return None

        def at(x):
            p = self._at_MGUT(a1inv_MI, a2inv_MI, a3inv_MI, MI, x)
            if p is None or not np.isfinite(p["full"]).all() or np.min(p["full"]) <= 0:
                return None
            return p

        left, right = at(lo), at(hi)
        if left is None or right is None:
            return None
        for p in (left, right):
            if abs(p["gap23"]) < tol:
                return p
        if left["gap23"] * right["gap23"] > 0:
            return None
        # Interpolation is fast for smooth gauge running; every fifth step
        # bisects so stagnation cannot masquerade as convergence.
        for iteration in range(max(60, max_iter)):
            fl, fr = left["gap23"], right["gap23"]
            x = (lo + hi) / 2 if iteration % 5 == 4 else hi - fr*(hi-lo)/(fr-fl)
            if not lo < x < hi:
                x = (lo + hi) / 2
            p = at(x)
            if p is None:
                return None
            if abs(p["gap23"]) < tol:
                return p
            if fl*p["gap23"] > 0:
                lo, left = x, p
            else:
                hi, right = x, p
        return None

    # -- one trial M_I -------------------------------------------------------
    def evaluate(self, logMI):
        MI = 10 ** logMI
        a1inv_MI, a2inv_MI, a3inv_MI = self.sm_to_MI(MI)
        pt = self.locate_MGUT(a1inv_MI, a2inv_MI, a3inv_MI, MI)
        if pt is None or pt["MGUT"] <= MI:
            return None
        MGUT, logMGUT, Delta4, r, full = pt["MGUT"], pt["logMGUT"], pt["Delta4"], pt["r"], pt["full"]
        a3inv_G, a2inv_G, aRinv_G, aBLinv_G = full
        a1inv_G = 0.6 * aRinv_G + 0.4 * aBLinv_G
        # Residual over every coupling that must coincide at M_GUT.  In
        # "unification" mode alpha_2R = alpha_BL = alpha_1 by construction and
        # this reduces to the original (alpha_1, alpha_2L, alpha_3) metric; in
        # "lr_parity" mode the alpha_2R/alpha_BL mismatch is also counted.
        max_frac, spread, mean = mismatch_metrics_n([a1inv_G, a2inv_G, a3inv_G, aRinv_G, aBLinv_G])
        max_frac_3, _, _ = mismatch_metrics(a1inv_G, a2inv_G, a3inv_G)
        alpha_GUT = 1.0 / mean
        tau = proton_decay_tau_years(MGUT, alpha_GUT, self.cfg["alpha_H_GeV3"])
        _, _, aRinv_MI, aBLinv_MI = match_SM_to_LR_with_ratio(a1inv_MI, a2inv_MI, a3inv_MI, r_ratio=r)
        return {
            "log10_MI": logMI, "MI_GeV": MI,
            "log10_MGUT": logMGUT, "MGUT_GeV": MGUT,
            "split_r_solved": r, "split_mode": self.split_mode,
            "alpha_inv_at_MI": {"alpha_1": a1inv_MI, "alpha_2L": a2inv_MI, "alpha_3": a3inv_MI,
                                "alpha_2R": aRinv_MI, "alpha_BL": aBLinv_MI},
            "max_frac": max_frac, "max_frac_3": max_frac_3, "spread": spread,
            "alpha_GUT": alpha_GUT,
            "alpha_inv_at_MGUT": {"alpha_1": a1inv_G, "alpha_2L": a2inv_G, "alpha_3": a3inv_G,
                                  "alpha_2R": aRinv_G, "alpha_BL": aBLinv_G},
            "thresholds_Delta_inv_alpha": {"alpha_1": 0.6 * Delta4[2] + 0.4 * Delta4[3],
                                           "alpha_2L": Delta4[1], "alpha_3": Delta4[0],
                                           "alpha_2R": Delta4[2], "alpha_BL": Delta4[3]},
            "tau_p_eppi0_years": tau,
        }


def golden_refine(model, lo, hi, iters=24):
    """Golden-section minimisation of the residual in log10(M_I) on [lo, hi]."""
    g = (math.sqrt(5.0) - 1.0) / 2.0

    def val(x):
        res = model.evaluate(x)
        return (1e9, None) if res is None else (res["max_frac"], res)

    c = hi - g * (hi - lo)
    d = lo + g * (hi - lo)
    fc, rc = val(c)
    fd, rd = val(d)
    for _ in range(iters):
        if fc < fd:
            hi, d, fd, rd = d, c, fc, rc
            c = hi - g * (hi - lo)
            fc, rc = val(c)
        else:
            lo, c, fc, rc = c, d, fd, rd
            d = lo + g * (hi - lo)
            fd, rd = val(d)
    return rc if fc < fd else rd


def run_scan(cfg, spec, args):
    model = Model(cfg, spec, two_loop=not args.one_loop, rk_steps=args.rk_steps)

    log10_min = args.mi_min if args.mi_min is not None else cfg.get("scan_log10_MI_min", 4.0)
    log10_max = args.mi_max if args.mi_max is not None else cfg.get("scan_log10_MI_max", 15.0)
    npts = args.scan_points if args.scan_points is not None else cfg.get("scan_points", 120)

    os.makedirs(args.outdir, exist_ok=True)
    csv_path = os.path.join(args.outdir, "scan.csv")
    rows, best = [], None
    grid = np.linspace(log10_min, log10_max, npts)
    for logMI in grid:
        res = model.evaluate(float(logMI))
        if res is None:
            continue
        rows.append(res)
        if best is None or res["max_frac"] < best["max_frac"]:
            best = res

    if best is not None and npts > 1:
        h = (log10_max - log10_min) / (npts - 1)
        refined = golden_refine(model, max(log10_min, best["log10_MI"] - h), min(log10_max, best["log10_MI"] + h))
        if refined is not None and refined["max_frac"] < best["max_frac"]:
            best = refined

    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["log10_MI", "log10_MGUT", "split_r_solved", "max_frac_mismatch", "alpha_GUT", "tau_p_eppi0_years"])
        for r in rows:
            w.writerow([r["log10_MI"], r["log10_MGUT"], r["split_r_solved"], r["max_frac"], r["alpha_GUT"], r["tau_p_eppi0_years"]])

    names = lr_betas.GROUPS[spec["group"]]["names"]
    summary = {
        "mode": "spectrum_derived",
        "spectrum_file": os.path.abspath(args.spectrum),
        "spectrum_name": spec.get("name"),
        "two_loop": not args.one_loop,
        "rk_steps": args.rk_steps,
        "derived_lr_betas_one_loop": {n: str(x) for n, x in zip(names, model.bLR_exact)},
        "derived_lr_betas_two_loop": {n: [str(x) for x in row] for n, row in zip(names, model.bijLR_exact)},
        "gut_threshold_multiplets": [{"name": H["name"], "b": [str(x) for x in H["b_exact"]], "mass_GeV": H["mass_GeV"]} for H in model.heavy],
        "free_parameters": [],
        "split_at_MI": model.split_mode,
        "scanned": {"log10_MI": [log10_min, log10_max, npts]},
        "best": best,
        "how_to_read_residual": (
            "unification mode: the U(1) split is solved from alpha_2R = alpha_BL at M_GUT, so the system is exactly "
            "determined; a near-zero residual means a self-consistent solution exists in range and the outputs are the "
            "predicted M_I, M_GUT, alpha_GUT and alpha_2R(M_I)."
            if model.split_mode == "unification" else
            "lr_parity mode: g_L = g_R at M_I is imposed. For L/R-symmetric spectra and matching, equality "
            "persists under running and is redundant; over-determination occurs only when it adds an independent condition."
        ),
        "caveats": [
            "tau_p is a tree-level dimensional estimate (no matrix elements, no operator running); M_X = M_GUT assumed",
            "two-loop LR coefficients are gauge-only (no Yukawa contributions)",
            "thresholds are one-loop leading-log step functions; no finite scheme constants",
            "matching at M_I is tree-level (no W_R / Delta_R / nu_R / heavy-bidoublet thresholds)",
            "M_GUT requires a converged alpha_2L = alpha_3 crossing inside root_bracket_log10; absent crossings are rejected",
            "quantum numbers and heavy masses are exactly those listed in the spectrum file",
        ],
    }
    with open(os.path.join(args.outdir, "best_fit.json"), "w") as f:
        json.dump(summary, f, indent=2)
    return model, summary, csv_path


def main():
    ap = argparse.ArgumentParser(
        description="LR gauge-coupling scan with spectrum-derived betas and thresholds (no free b's, no manual dalpha).")
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--spectrum", required=True, help="lr_spectrum_v1 JSON file (representations + masses)")
    ap.add_argument("--mi-min", type=float, default=None, help="log10 of lowest M_I on the grid")
    ap.add_argument("--mi-max", type=float, default=None, help="log10 of highest M_I on the grid")
    ap.add_argument("--scan-points", type=int, default=None)
    ap.add_argument("--one-loop", action="store_true", help="disable two-loop running (default: two-loop)")
    ap.add_argument("--rk-steps", type=int, default=400, help="RK4 steps per running segment (results are unchanged at 4 significant figures between 200 and 3000)")
    ap.add_argument("--outdir", default="out")
    # Retired knobs: refuse loudly rather than silently ignore.
    import sys
    retired = ("--b3LR", "--b2L", "--b2R", "--bBL", "--split-r", "--dalpha1", "--dalpha2", "--dalpha3", "--multiplets")
    used = [a for a in sys.argv[1:] if a.split("=", 1)[0] in retired]
    if used:
        ap.error(f"{', '.join(used)}: removed. Betas and thresholds are derived from --spectrum and the U(1) "
                 "split is fixed by the spectrum's split_at_MI setting; there are no tunable knobs.")
    args = ap.parse_args()

    cfg = load_config(args.config)
    spec = lr_betas.load_spectrum(args.spectrum)
    print(lr_betas.describe_spectrum(spec))
    print()

    model, summary, csv_path = run_scan(cfg, spec, args)
    best = summary["best"]

    print("=== Scan complete (v5, spectrum-derived mode) ===")
    if best is None:
        print("No physical point found: alpha_2L = alpha_3 crossing or positive-coupling split does not exist in range.")
    else:
        print(f"Best M_I   ≈ 10^{best['log10_MI']:.3f} GeV")
        print(f"Best M_GUT ≈ 10^{best['log10_MGUT']:.3f} GeV")
        print(f"Residual mismatch @ M_GUT: {best['max_frac']*100:.4f}%")
        print(f"alpha_GUT ≈ {best['alpha_GUT']:.5f}   (1/alpha_GUT ≈ {1/best['alpha_GUT']:.3f})")
        print(f"tau_p (p→ e^+π^0) ≈ {best['tau_p_eppi0_years']:.2e} years  [tree-level estimate]")
        how = "solved from alpha_2R = alpha_BL at M_GUT" if best["split_mode"] == "unification" else "fixed by g_L = g_R at M_I"
        print(f"split 1/alpha_2R : 1/alpha_BL at M_I = {best['split_r_solved']:.5f}  ({how})")
        am = best["alpha_inv_at_MI"]
        print("1/alpha at M_I  : " + ", ".join(f"{k}={v:.4f}" for k, v in am.items()))
        ai = best["alpha_inv_at_MGUT"]
        print("1/alpha at M_GUT (thresholds applied): " + ", ".join(f"{k}={v:.4f}" for k, v in ai.items()))
        th = best["thresholds_Delta_inv_alpha"]
        print("threshold Delta(1/alpha) from spectrum: " + ", ".join(f"{k}={v:+.4f}" for k, v in th.items()))
        MR = best["MI_GeV"]
        for y in [1.0, 0.3, 0.1, 0.03]:
            print(f"Seesaw mν for y={y:.2f}: ~{seesaw_mnu_eV(y, 174.0, MR):.3e} eV (M_R ≈ M_I)")
    print("\nFree parameters tuned in this run: none.")
    print("Reading the residual: " + summary["how_to_read_residual"])
    print("Approximations: " + "; ".join(summary["caveats"]))
    print(f"\nCSV written to : {csv_path}")
    print(f"Summary JSON   : {os.path.join(args.outdir, 'best_fit.json')}")


if __name__ == "__main__":
    main()

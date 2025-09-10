import json, math, csv, os, argparse
import numpy as np
from rge_core import gut_normalized_alphas_at_MZ, match_SM_to_LR_with_ratio, find_MGUT_by_alpha2_eq_alpha3, mismatch_metrics
from two_loop_core import run_segment_two_loop, run_segment_one_loop

def load_config(path):
    with open(path, "r") as f:
        return json.load(f)

def proton_decay_tau_years(MX_GeV, alpha_GUT, alpha_H_GeV3):
    return (1e35 * (MX_GeV/1e16)**4 * (0.04/alpha_GUT)**2 * (0.012/alpha_H_GeV3)**2)

def seesaw_mnu_eV(y, v_GeV, MR_GeV):
    return (y*v_GeV)**2 / MR_GeV * 1e9

def sm_two_loop_defaults():
    # b vector (U1, SU2, SU3) and b_ij matrix in GUT-normalized basis for SM (n_g=3, n_H=1)
    b = np.array([41.0/10.0, -19.0/6.0, -7.0])
    bij = np.array([
        [199.0/50.0, 27.0/10.0, 44.0/5.0],
        [9.0/10.0,   35.0/6.0,  12.0     ],
        [11.0/10.0,  9.0/2.0,  -26.0     ]
    ])
    return b, bij

def lr_two_loop_defaults():
    # Placeholder LR (alpha3, alpha2L, alpha2R, alphaBL) two-loop coefficients -> zeros by default
    b = None; bij = None
    return b, bij

def compute_multiplet_thresholds(multiplets, MGUT, group='SM', n_gauges=3):
    """Compute Δ_i = (1/12π) sum_H b_i^(H) ln(MGUT / M_H) for given multiplets.
    Each entry: {"name":..., "group":"SM" or "LR", "b": [..], "mass_GeV": ...}
    Returns numpy array Δ of length n_gauges (missing entries treated as zeros).
    """
    import math, numpy as np
    two_pi = 2*math.pi
    Delta = np.zeros(n_gauges, dtype=float)
    if not multiplets:
        return Delta
    for H in multiplets:
        if H.get("group","SM") != group:
            continue
        bH = H.get("b", [])
        mH = H.get("mass_GeV", MGUT)
        if len(bH) != n_gauges or mH <= 0:
            continue
        Delta += (np.array(bH, dtype=float) / (12.0*math.pi)) * math.log(MGUT / mH)
    return Delta

def run_scan(cfg, args):
    # inputs at MZ
    MZ = cfg["MZ_GeV"]
    alpha1_MZ, alpha2_MZ, alpha3_MZ = gut_normalized_alphas_at_MZ(cfg["alpha_em_MZ"], cfg["sin2thetaW_MZ"], cfg["alpha3_MZ"])
    a1inv_MZ, a2inv_MZ, a3inv_MZ = 1/alpha1_MZ, 1/alpha2_MZ, 1/alpha3_MZ

    # SM 1-loop b's
    b1_SM, b2_SM, b3_SM = cfg["b1_SM"], cfg["b2_SM"], cfg["b3_SM"]

    # LR b's
    b3_LR = args.b3LR if args.b3LR is not None else cfg.get("b3_LR", -7.0)
    b2L_LR = args.b2L if args.b2L is not None else cfg.get("b2L_LR", -3.0)
    b2R_LR = args.b2R if args.b2R is not None else cfg.get("b2R_LR", -3.0)
    bBL_LR = args.bBL if args.bBL is not None else cfg.get("bBL_LR", 4.0)

    # Two-loop flags
    use_two_loop = args.two_loop

    # Load multiplets if provided
    multiplets = None
    if args.multiplets:
        with open(args.multiplets, "r") as f:
            multiplets = json.load(f)

    # MI scan bounds
    log10_min = args.mi_min if args.mi_min is not None else cfg.get("scan_log10_MI_min", 9.0)
    log10_max = args.mi_max if args.mi_max is not None else cfg.get("scan_log10_MI_max", 15.0)
    npts = args.scan_points if args.scan_points is not None else cfg.get("scan_points", 120)

    # Threshold manual shifts
    d1 = args.dalpha1 or 0.0
    d2 = args.dalpha2 or 0.0
    d3 = args.dalpha3 or 0.0

    # matching split ratio
    r_ratio = args.split_r if args.split_r is not None else 1.0

    os.makedirs("out", exist_ok=True)
    csv_path = "out/scan.csv"
    rows = []
    best = {"max_frac": 1e9}

    # Two-loop coefficients
    import numpy as np
    bSM_1, bSM_2 = sm_two_loop_defaults()
    bLR_1, bLR_2 = lr_two_loop_defaults()
    if bLR_1 is None:
        bLR_1 = np.array([b3_LR, b2L_LR, b2R_LR, bBL_LR])
        bLR_2 = np.zeros((4,4))

    for logMI in np.linspace(log10_min, log10_max, npts):
        MI = 10**logMI

        if not use_two_loop:
            # 1-loop SM to MI
            a_inv_SM = np.array([a1inv_MZ, a2inv_MZ, a3inv_MZ])
            a_inv_MI = run_segment_one_loop(a_inv_SM, np.array([b1_SM,b2_SM,b3_SM]), MZ, MI)
            a1inv_MI, a2inv_MI, a3inv_MI = a_inv_MI.tolist()
        else:
            # two-loop SM to MI
            alphas_SM = np.array([alpha1_MZ, alpha2_MZ, alpha3_MZ])
            alphas_MI = run_segment_two_loop(alphas_SM, bSM_1, bSM_2, MZ, MI, nsteps=3000)
            a1inv_MI, a2inv_MI, a3inv_MI = (1/alphas_MI[0], 1/alphas_MI[1], 1/alphas_MI[2])

        # Match to LR at MI (split ratio)
        a3inv_LR_MI, aLinv_LR_MI, aRinv_LR_MI, aBLinv_LR_MI = match_SM_to_LR_with_ratio(a1inv_MI, a2inv_MI, a3inv_MI, r_ratio=r_ratio)

        # Find MGUT via alpha2==alpha3 (1-loop relation in inverse space is fine as locator)
        MGUT, logMGUT = find_MGUT_by_alpha2_eq_alpha3(aLinv_LR_MI, a3inv_LR_MI, b2L_LR, b3_LR, MI, (14.2, 17.9))
        if MGUT is None:
            continue

        if not use_two_loop:
            # 1-loop LR to MGUT
            aRinv_MGUT = run_segment_one_loop(aRinv_LR_MI, np.array(b2R_LR), MI, MGUT)
            aBLinv_MGUT = run_segment_one_loop(aBLinv_LR_MI, np.array(bBL_LR), MI, MGUT)
            a2inv_MGUT  = run_segment_one_loop(aLinv_LR_MI, np.array(b2L_LR), MI, MGUT)
            a3inv_MGUT  = run_segment_one_loop(a3inv_LR_MI, np.array(b3_LR), MI, MGUT)
            a1inv_MGUT  = (3/5)*aRinv_MGUT + (2/5)*aBLinv_MGUT
        else:
            # two-loop LR to MGUT (alphas formulation)
            # reconstruct alphas at MI from inverses (approx)
            alphaR_MI = 1.0/aRinv_LR_MI
            alphaBL_MI = 1.0/aBLinv_LR_MI
            alphaL_MI = 1.0/aLinv_LR_MI
            alpha3_MI = 1.0/a3inv_LR_MI
            alphas_MI_LR = np.array([alpha3_MI, alphaL_MI, alphaL_MI*0 + alphaR_MI, alphaBL_MI]) # (3,2L,2R,BL)
            alphas_MGUT = run_segment_two_loop(alphas_MI_LR, bLR_1, bLR_2, MI, MGUT, nsteps=3000)
            alpha3_MG = alphas_MGUT[0]; alpha2L_MG = alphas_MGUT[1]; alpha2R_MG = alphas_MGUT[2]; alphaBL_MG = alphas_MGUT[3]
            a3inv_MGUT = 1.0/alpha3_MG
            a2inv_MGUT = 1.0/alpha2L_MG
            a1inv_MGUT = (3/5)*(1.0/alpha2R_MG) + (2/5)*(1.0/alphaBL_MG)

        # Apply physical multiplet thresholds mapped to SM (1,2,3)
        Delta_phys = np.zeros(3)
        if multiplets:
            from two_loop_core import np as _np
            Delta_phys = compute_multiplet_thresholds(multiplets, MGUT, group='SM', n_gauges=3)
        # manual Δ's
        a1inv_MGUT += (d1 + float(Delta_phys[0]))
        a2inv_MGUT += (d2 + float(Delta_phys[1]))
        a3inv_MGUT += (d3 + float(Delta_phys[2]))

        # Mismatch & outputs
        max_frac, spread, mean = mismatch_metrics(a1inv_MGUT, a2inv_MGUT, a3inv_MGUT)
        alpha_GUT = 1.0/mean
        tau_years = proton_decay_tau_years(MX_GeV=MGUT, alpha_GUT=alpha_GUT, alpha_H_GeV3=cfg["alpha_H_GeV3"])

        rows.append({
            "log10_MI": logMI,
            "log10_MGUT": logMGUT,
            "max_frac_mismatch": max_frac,
            "alpha_GUT": alpha_GUT,
            "tau_p_eppi0_years": tau_years
        })

        if max_frac < best["max_frac"]:
            best = {
                "log10_MI": logMI,
                "MI_GeV": 10**logMI,
                "log10_MGUT": logMGUT,
                "MGUT_GeV": MGUT,
                "max_frac": max_frac,
                "alpha_GUT": alpha_GUT,
                "tau_p_eppi0_years": tau_years,
                "params": {
                    "b3_LR": b3_LR, "b2L_LR": b2L_LR, "b2R_LR": b2R_LR, "bBL_LR": bBL_LR,
                    "split_r": r_ratio, "dalpha1": d1, "dalpha2": d2, "dalpha3": d3,
                    "two_loop": use_two_loop, "multiplets_used": bool(multiplets)
                }
            }

    # write CSV and JSON
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["log10_MI","log10_MGUT","max_frac_mismatch","alpha_GUT","tau_p_eppi0_years"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open("out/best_fit.json","w") as f:
        json.dump(best, f, indent=2)

    return best, csv_path

def main():
    ap = argparse.ArgumentParser(description="v5: LR scan with optional two-loop running and physical multiplet thresholds")
    ap.add_argument("--config", default="config.json")
    # LR b's and matching
    ap.add_argument("--b3LR", type=float, default=None)
    ap.add_argument("--b2L", type=float, default=None)
    ap.add_argument("--b2R", type=float, default=None)
    ap.add_argument("--bBL", type=float, default=None)
    ap.add_argument("--split-r", dest="split_r", type=float, default=None)
    # Scan control
    ap.add_argument("--mi-min", type=float, default=None)
    ap.add_argument("--mi-max", type=float, default=None)
    ap.add_argument("--scan-points", type=int, default=None)
    # Thresholds and extras
    ap.add_argument("--dalpha1", type=float, default=0.0)
    ap.add_argument("--dalpha2", type=float, default=0.0)
    ap.add_argument("--dalpha3", type=float, default=0.0)
    ap.add_argument("--multiplets", type=str, default=None, help="JSON file of heavy multiplets with beta contributions and masses")
    ap.add_argument("--two-loop", dest="two_loop", action="store_true", help="Enable two-loop running (SM and LR gauge-only)")

    args = ap.parse_args()
    cfg = load_config(args.config)
    best, csv_path = run_scan(cfg, args)

    print("\n=== Scan complete (v5) ===")
    if best.get("max_frac", 1e9) < 1e8:
        print(f"Best MI  \u2248 10^{best['log10_MI']:.3f} GeV")
        print(f"Best MGUT\u2248 10^{best['log10_MGUT']:.3f} GeV")
        print(f"Residual mismatch @ MGUT: {best['max_frac']*100:.2f}%")
        print(f"alpha_GUT \u2248 {best['alpha_GUT']:.4f}")
        print(f"tau_p (p\u2192 e^+\u03c0^0) \u2248 {best['tau_p_eppi0_years']:.2e} years")
        p = best.get("params", {})
        if p:
            print("\nSettings:")
            for k,v in p.items():
                print(f"  {k}: {v}")
    else:
        print("No solution found; adjust knobs/coefficients or ranges.")

    # Seesaw band at best MI
    MR = 10**best["log10_MI"]
    v = 174.0
    for y in [1.0, 0.3, 0.1, 0.03]:
        mnu = (y*v)**2 / MR * 1e9
        print(f"Seesaw m\u03bd for y={y:.2f}: ~{mnu:.3e} eV (MR\u2248MI\u2248 10^{best['log10_MI']:.2f} GeV)")

    print(f"\nCSV written to: {csv_path}")
    print("Best fit JSON : out/best_fit.json")

if __name__ == '__main__':
    main()

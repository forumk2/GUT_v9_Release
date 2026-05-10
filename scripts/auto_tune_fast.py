import argparse, json, math, os, time, random
import numpy as np
from rge_core import gut_normalized_alphas_at_MZ, match_SM_to_LR_with_ratio, find_MGUT_by_alpha2_eq_alpha3, mismatch_metrics
from two_loop_core import run_segment_two_loop, run_segment_one_loop
from rge_scan_v5 import sm_two_loop_defaults, lr_two_loop_defaults, compute_multiplet_thresholds

def proton_decay_tau_years(MX_GeV, alpha_GUT, alpha_H_GeV3):
    return (1e35 * (MX_GeV/1e16)**4 * (0.04/alpha_GUT)**2 * (0.012/alpha_H_GeV3)**2)
def load_config(path): return json.load(open(path))

def run_point(cfg, bLR, split_r, logMI, deltas, use_two_loop, multiplets, rk_steps):
    MZ = cfg["MZ_GeV"]
    alpha1_MZ, alpha2_MZ, alpha3_MZ = gut_normalized_alphas_at_MZ(cfg["alpha_em_MZ"], cfg["sin2thetaW_MZ"], cfg["alpha3_MZ"])
    a1inv_MZ, a2inv_MZ, a3inv_MZ = 1/alpha1_MZ, 1/alpha2_MZ, 1/alpha3_MZ

    MI = 10**logMI
    if not use_two_loop:
        a_inv_MI = run_segment_one_loop(np.array([a1inv_MZ,a2inv_MZ,a3inv_MZ]), np.array([cfg["b1_SM"],cfg["b2_SM"],cfg["b3_SM"]]), MZ, MI)
        a1inv_MI, a2inv_MI, a3inv_MI = a_inv_MI.tolist()
    else:
        bSM_1, bSM_2 = sm_two_loop_defaults()
        alphas_MI = run_segment_two_loop(np.array([alpha1_MZ,alpha2_MZ,alpha3_MZ]), bSM_1, bSM_2, MZ, MI, nsteps=rk_steps)
        a1inv_MI, a2inv_MI, a3inv_MI = (1/alphas_MI[0], 1/alphas_MI[1], 1/alphas_MI[2])

    a3inv_LR_MI, aLinv_LR_MI, aRinv_LR_MI, aBLinv_LR_MI = match_SM_to_LR_with_ratio(a1inv_MI, a2inv_MI, a3inv_MI, r_ratio=split_r)
    MGUT, logMGUT = find_MGUT_by_alpha2_eq_alpha3(aLinv_LR_MI, a3inv_LR_MI, bLR["b2L"], bLR["b3"], MI, (14.2,17.9))
    if MGUT is None: return None

    if not use_two_loop:
        aRinv_MGUT = run_segment_one_loop(aRinv_LR_MI, np.array(bLR["b2R"]), MI, MGUT)
        aBLinv_MGUT = run_segment_one_loop(aBLinv_LR_MI, np.array(bLR["bBL"]), MI, MGUT)
        a2inv_MGUT  = run_segment_one_loop(aLinv_LR_MI, np.array(bLR["b2L"]), MI, MGUT)
        a3inv_MGUT  = run_segment_one_loop(a3inv_LR_MI, np.array(bLR["b3"]), MI, MGUT)
        a1inv_MGUT  = (3/5)*aRinv_MGUT + (2/5)*aBLinv_MGUT
    else:
        bLR_1, bLR_2 = lr_two_loop_defaults()
        if bLR_1 is None:
            bLR_1 = np.array([bLR["b3"], bLR["b2L"], bLR["b2R"], bLR["bBL"]])
            bLR_2 = np.zeros((4,4))
        alphaR_MI = 1.0/aRinv_LR_MI; alphaBL_MI = 1.0/aBLinv_LR_MI; alphaL_MI = 1.0/aLinv_LR_MI; alpha3_MI = 1.0/a3inv_LR_MI
        alphas_MGUT = run_segment_two_loop(np.array([alpha3_MI, alphaL_MI, alphaR_MI, alphaBL_MI]), bLR_1, bLR_2, MI, MGUT, nsteps=rk_steps)
        a3inv_MGUT = 1.0/alphas_MGUT[0]; a2inv_MGUT = 1.0/alphas_MGUT[1]
        a1inv_MGUT = (3/5)*(1.0/alphas_MGUT[2]) + (2/5)*(1.0/alphas_MGUT[3])

    Delta_phys = np.zeros(3)
    if multiplets:
        Delta_phys = compute_multiplet_thresholds(multiplets, MGUT, group='SM', n_gauges=3)
    d1,d2,d3 = deltas
    a1inv_MGUT += d1 + float(Delta_phys[0])
    a2inv_MGUT += d2 + float(Delta_phys[1])
    a3inv_MGUT += d3 + float(Delta_phys[2])

    max_frac, spread, mean = mismatch_metrics(a1inv_MGUT, a2inv_MGUT, a3inv_MGUT)
    alpha_GUT = 1.0/mean
    tau = proton_decay_tau_years(MX_GeV=MGUT, alpha_GUT=alpha_GUT, alpha_H_GeV3=cfg["alpha_H_GeV3"])
    return {"log10_MI": logMI, "log10_MGUT": logMGUT, "alpha_GUT": alpha_GUT, "mismatch": max_frac, "tau": tau, "MGUT_GeV": MGUT}

def main():
    ap = argparse.ArgumentParser(description="Fast random-sampling autotune with progress")
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--two-loop", dest="two_loop", action="store_true")
    ap.add_argument("--multiplets", type=str, default=None)
    # LR + split
    ap.add_argument("--b3LR", type=float, default=-7.0)
    ap.add_argument("--b2L", type=float, default=-2.0)
    ap.add_argument("--b2R", type=float, default=-2.0)
    ap.add_argument("--bBL", type=float, default=5.0)
    ap.add_argument("--split-r", dest="split_r", type=float, default=1.0)
    # Random sampling box
    ap.add_argument("--mi-min", type=float, default=10.6)
    ap.add_argument("--mi-max", type=float, default=12.8)
    ap.add_argument("--dmin", type=float, default=-0.6)
    ap.add_argument("--dmax", type=float, default=0.6)
    ap.add_argument("--trials", type=int, default=8000)
    ap.add_argument("--seed", type=int, default=42)
    # Targets & weights
    ap.add_argument("--tau-min", type=float, default=1.0e35)
    ap.add_argument("--mnu-target", type=float, default=0.05)
    ap.add_argument("--y-ref", type=float, default=0.1)
    ap.add_argument("--mismatch-wt", type=float, default=10.0)
    ap.add_argument("--mnu-wt", type=float, default=1.0)
    ap.add_argument("--tau-wt", type=float, default=1.0)
    # RK steps (speed/accuracy knob)
    ap.add_argument("--rk-steps", type=int, default=1200)
    args = ap.parse_args()

    random.seed(args.seed); np.random.seed(args.seed)
    cfg = load_config(args.config)
    multiplets = json.load(open(args.multiplets)) if args.multiplets else None
    bLR = {"b3": args.b3LR, "b2L": args.b2L, "b2R": args.b2R, "bBL": args.bBL}

    os.makedirs("out/auto", exist_ok=True)
    best = None; best_score = 1e99
    last_print = time.time()

    for t in range(1, args.trials+1):
        logMI = random.uniform(args.mi_min, args.mi_max)
        d1 = random.uniform(args.dmin, args.dmax)
        d2 = random.uniform(args.dmin, args.dmax)
        d3 = random.uniform(args.dmin, args.dmax)

        # neutrino objective from MI (no need to recompute inside run_point)
        MR = 10**logMI; v = 174.0
        mnu_ref = (args.y_ref*v)**2 / MR * 1e9
        nu_obj = abs(math.log10(max(mnu_ref,1e-30)/args.mnu_target))

        res = run_point(cfg, bLR, args.split_r, logMI, (d1,d2,d3), args.two_loop, multiplets, args.rk_steps)
        if res is None: 
            continue
        if res["tau"] < args.tau_min:
            # Discard but occasionally print to show activity
            if time.time() - last_print > 2.0 and t % 200 == 0:
                print(f"[{t}/{args.trials}] τ too small ({res['tau']:.2e}); exploring...")
                last_print = time.time()
            continue

        score = args.mismatch_wt*res["mismatch"] + args.mnu_wt*nu_obj + args.tau_wt*(-math.log10(res["tau"]))
        if score < best_score:
            best_score = score
            best = {**res, "MI_GeV": 10**logMI, "dalpha1": d1, "dalpha2": d2, "dalpha3": d3,
                    "mnu_ref_eV": mnu_ref, "y_ref": args.y_ref,
                    "two_loop": args.two_loop, "multiplets_used": bool(multiplets)}
            json.dump(best, open("out/auto/best_autotune_fast.json","w"), indent=2)

        # progress print every ~2s or every 500 trials
        if time.time() - last_print > 2.0 or (t % 500 == 0):
            if best:
                print(f"[{t}/{args.trials}] best mismatch={best['mismatch']*100:.2f}%, "
                      f"log10(MI)={math.log10(best['MI_GeV']):.2f}, "
                      f"log10(MGUT)={res['log10_MGUT']:.2f}, τ≈{best['tau']:.2e}")
            else:
                print(f"[{t}/{args.trials}] searching...")
            last_print = time.time()

    if best:
        print("\n=== Auto-tune FAST complete ===")
        print(f"Best MI   ≈ 10^{math.log10(best['MI_GeV']):.3f} GeV")
        print(f"Best MGUT ≈ 10^{res['log10_MGUT']:.3f} GeV")
        print(f"Residual mismatch @ MGUT: {best['mismatch']*100:.2f}%")
        print(f"alpha_GUT ≈ {best['alpha_GUT']:.4f}")
        print(f"Proton lifetime τ_p ≈ {best['tau']:.2e} years (≥ {args.tau_min:.2e})")
        print(f"Thresholds Δ(1/α): ({best['dalpha1']:+.3f}, {best['dalpha2']:+.3f}, {best['dalpha3']:+.3f})")
        print(f"Seesaw mν (y={best['y_ref']}) at MI: ~{best['mnu_ref_eV']:.3e} eV (target {args.mnu_target} eV)")
        print("Saved: out/auto/best_autotune_fast.json")
    else:
        print("No feasible point found in fast mode. Try more trials or relax constraints.")

if __name__ == '__main__':
    main()

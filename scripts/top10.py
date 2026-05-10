#!/usr/bin/env python3
import csv, argparse, math, sys
from pathlib import Path

def read_rows(path):
    with open(path, newline='', encoding='utf-8') as f:
        r = csv.DictReader(f)
        for row in r:
            row['log10MGUT'] = float(row['log10MGUT'])
            row['log10MI'] = float(row['log10MI'])
            row['residual_pct'] = float(row['residual_pct'])
            row['alphaGUT'] = float(row['alphaGUT'])
            try:
                row['tau'] = float(row['tau'])
            except:
                row['tau'] = float(eval(row['tau']))
            yield row

def main():
    ap = argparse.ArgumentParser(description="Print top-N points from scan_v9_summary.csv")
    ap.add_argument("csv_path", nargs="?", default="out/v9/scan_v9_summary.csv")
    ap.add_argument("--max-residual", type=float, default=0.2, help="Percent")
    ap.add_argument("-N", type=int, default=10)
    ap.add_argument("--by", choices=["tau","mgut","residual"], default="tau")
    ap.add_argument("--ref-tag", default=None, help="Tag of reference point for tau scaling check")
    ap.add_argument("--scaling-check", action="store_true", help="Compare tau to (M/Mref)^4 * (alpha_ref/alpha)^2 scaling")
    args = ap.parse_args()

    rows = [r for r in read_rows(args.csv_path) if r['residual_pct'] <= args.max_residual]
    if not rows:
        print("No rows pass residual threshold.", file=sys.stderr); sys.exit(1)

    key = {"tau": lambda r: r['tau'],
           "mgut": lambda r: r['log10MGUT'],
           "residual": lambda r: -r['residual_pct']}[args.by]
    rows.sort(key=key, reverse=True)

    print(f"# Top {min(args.N,len(rows))} (residual ≤ {args.max_residual:.3f}%) sorted by {args.by}")
    print("tag, log10MI, log10MGUT, residual%, alphaGUT, tau")
    for r in rows[:args.N]:
        print(f"{r['tag']}, {r['log10MI']:.3f}, {r['log10MGUT']:.3f}, {r['residual_pct']:.3f}, {r['alphaGUT']:.4f}, {r['tau']:.3e}")

    if args.scaling_check:
        if not args.ref_tag:
            print("\n[scaling] --ref-tag required for scaling check", file=sys.stderr); sys.exit(2)
        ref = next((r for r in rows if r['tag']==args.ref_tag), None)
        if not ref:
            print(f"\n[scaling] reference tag not found: {args.ref_tag}", file=sys.stderr); sys.exit(3)
        Mref = 10**ref['log10MGUT']; aref = ref['alphaGUT']; tref = ref['tau']
        print(f"\n[scaling] reference: {ref['tag']}  Mref={Mref:.3e}, alpha_ref={aref:.5f}, tau_ref={tref:.3e}")
        print("tag, predicted_tau_from_scaling, actual_tau, ratio(actual/pred)")
        for r in rows[:args.N]:
            M = 10**r['log10MGUT']; a = r['alphaGUT']
            tau_pred = tref * (M/Mref)**4 * (aref/a)**2
            print(f"{r['tag']}, {tau_pred:.3e}, {r['tau']:.3e}, {r['tau']/tau_pred:.3f}")

if __name__ == "__main__":
    main()

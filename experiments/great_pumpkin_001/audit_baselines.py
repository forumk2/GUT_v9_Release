"""Reproduce archived points with pinned historical code in a temporary dir."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
BASE = "8d44f8c5ca6eed3a83572471dac6f9229743ec67"


def historical():
    names = ["rge_core.py", "two_loop_core.py", "rge_scan_v5.py", "auto_tune_fast.py",
             "config.json", "spectrum_balanced_v8.json"]
    with tempfile.TemporaryDirectory() as temp:
        folder = Path(temp)
        hashes = {}
        for name in names:
            data = subprocess.check_output(["git", "show", f"{BASE}^:scripts/{name}"], cwd=ROOT)
            (folder/name).write_bytes(data)
            hashes[name] = hashlib.sha256(data).hexdigest()
        archives = [json.loads((ROOT / f"scripts/results/validate/20250909_203538/winner_{s}_best_fit.json").read_text())
                    for s in ["A", "B"]]
        (folder/"references.json").write_text(json.dumps(archives))
        code = '''import json,numpy as np
from auto_tune_fast import run_point
from rge_core import gut_normalized_alphas_at_MZ,match_SM_to_LR_with_ratio
from rge_scan_v5 import sm_two_loop_defaults
from two_loop_core import run_segment_two_loop,run_segment_one_loop
cfg=json.load(open('config.json')); heavy=json.load(open('spectrum_balanced_v8.json')); rows=[]
for ref in json.load(open('references.json')):
 p=ref['params']; b={'b3':p['b3_LR'],'b2L':p['b2L_LR'],'b2R':p['b2R_LR'],'bBL':p['bBL_LR']}
 shifts=[p['dalpha1'],p['dalpha2'],p['dalpha3']]
 out=run_point(cfg,b,p['split_r'],ref['log10_MI'],shifts,True,heavy,3000)
 alt=run_point(cfg,b,2.,ref['log10_MI'],shifts,True,heavy,3000)
 sm=run_segment_two_loop(np.array(gut_normalized_alphas_at_MZ(cfg['alpha_em_MZ'],cfg['sin2thetaW_MZ'],cfg['alpha3_MZ'])),*sm_two_loop_defaults(),cfg['MZ_GeV'],ref['MI_GeV'],nsteps=3000)
 initial=np.array(match_SM_to_LR_with_ratio(*(1/sm),p['split_r']))
 full=run_segment_one_loop(initial,np.array(list(b.values())),ref['MI_GeV'],out['MGUT_GeV'])
 rows.append({'reference':ref,'reproduced':out,'split_r_2':alt,'raw_LR_inverse_at_GUT':full.tolist(),'raw_R_minus_BL':float(full[2]-full[3])})
print(json.dumps(rows))
'''
        (folder/"reproduce.py").write_text(code)
        rows = json.loads(subprocess.check_output([sys.executable, "reproduce.py"], cwd=folder))
    return {"revision": BASE + "^", "source_sha256": hashes, "points": rows}


def current():
    # Pin baseline code too: later fixes cannot rewrite the before evidence.
    names = ["rge_core.py", "two_loop_core.py", "rge_scan_v5.py", "lr_betas.py",
             "config.json", "test_lr_betas.py"]
    specs = ["lr_minimal", "so10_10_126_thresholds", "lr_dparity_symmetric"]
    with tempfile.TemporaryDirectory() as temp:
        folder = Path(temp)
        for name in names + [f"spectra/{s}.json" for s in specs]:
            target = folder / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(subprocess.check_output(["git", "show", f"{BASE}:scripts/{name}"], cwd=ROOT))
        env = dict(os.environ, PYTHONIOENCODING="utf-8")
        good = subprocess.run([sys.executable, "test_lr_betas.py"], cwd=folder, capture_output=True, env=env)
        rows = []
        for s in specs:
            print("baseline scan:", s, flush=True)
            proc = subprocess.run([sys.executable, "rge_scan_v5.py", "--spectrum", f"spectra/{s}.json",
                                   "--outdir", f"out/{s}"], cwd=folder, env=env, capture_output=True, check=True)
            result = json.loads((folder/f"out/{s}/best_fit.json").read_text())
            result["spectrum_file"] = f"scripts/spectra/{s}.json"
            rows.append(result)
    return {"revision": BASE, "documented_cwd_test_exit": good.returncode,
            "documented_cwd_test_log": good.stdout.decode(), "scans": rows}


if __name__ == "__main__":
    out = Path(__file__).parent / "results"
    out.mkdir(exist_ok=True)
    report = {"historical": historical(), "current": current()}
    (out/"baselines.json").write_text(json.dumps(report, indent=2) + "\n")
    print("wrote", out/"baselines.json")

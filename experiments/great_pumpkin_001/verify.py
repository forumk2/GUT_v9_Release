"""Run portable tests and the real Bash validation entry point from repo root.

The Bash output lives in ignored scratch; only its result summary is evidence.
Use --quick for Python tests alone. Full mode requires Git Bash / bash.
"""
from pathlib import Path
import argparse
import json
import os
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--quick",action="store_true")
    args=ap.parse_args()
    env=dict(os.environ,PYTHONIOENCODING="utf-8",PYTEST_DISABLE_PLUGIN_AUTOLOAD="1")
    records=[]
    commands=[([sys.executable,"-m","pytest","scripts/test_lr_betas.py",
                "experiments/great_pumpkin_001/test_expedition.py","-q"],ROOT),
              ([sys.executable,"scripts/test_lr_betas.py"],ROOT),
              ([sys.executable,"test_lr_betas.py"],ROOT/"scripts")]
    for cmd,cwd in commands:
        proc=subprocess.run(cmd,cwd=cwd,env=env,capture_output=True,text=True,encoding="utf-8")
        records.append(dict(command=cmd[1:],cwd="root" if cwd==ROOT else "scripts",
                            exit_code=proc.returncode,stdout=proc.stdout,stderr=proc.stderr))
        print(proc.stdout,flush=True)
        if proc.returncode:
            raise RuntimeError(proc.stderr)
    if not args.quick:
        gitbash=Path("C:/Program Files/Git/bin/bash.exe")
        bash=str(gitbash) if os.name=="nt" and gitbash.exists() else shutil.which("bash")
        if bash is None:
            raise RuntimeError("Full verification requires Bash; use --quick for Python tests")
        scratch=HERE/"scratch"/"validation"
        env["VALIDATE_OUTDIR"]=scratch.as_posix()
        print("Running actual validate.sh from repository root...",flush=True)
        proc=subprocess.run([bash,"scripts/validate.sh"],cwd=ROOT,env=env,
                            capture_output=True,text=True,encoding="utf-8")
        (HERE/"scratch").mkdir(exist_ok=True)
        (HERE/"scratch"/"validate.log").write_text(proc.stdout+proc.stderr,encoding="utf-8")
        records.append(dict(command="bash scripts/validate.sh",cwd="root",exit_code=proc.returncode,
                            numerical_comparator_passed="All three numerical references match" in proc.stdout))
        if proc.returncode:
            raise RuntimeError(proc.stdout[-4000:]+proc.stderr)
        for name in ["lr_minimal","so10_10_126_thresholds","lr_dparity_symmetric"]:
            d=json.loads((scratch/name/"best_fit.json").read_text())
            records.append(dict(spectrum=name,best=d["best"]))
    out=HERE/"results"; out.mkdir(exist_ok=True)
    (out/("verification_quick.json" if args.quick else "verification.json")).write_text(
        json.dumps(records,indent=2)+"\n",encoding="utf-8")
    print("Verification passed.")


if __name__=="__main__":
    main()

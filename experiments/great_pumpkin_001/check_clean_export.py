"""Verify a fresh Git export with no local scratch or build products.

Committed historical evidence is present; validation writes new scratch output.
"""
from pathlib import Path
import io
import json
import os
import subprocess
import sys
import tempfile
import zipfile

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent


def main():
    revision=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()
    archive=subprocess.check_output(["git","archive","--format=zip",revision],cwd=ROOT)
    with tempfile.TemporaryDirectory(prefix="great-pumpkin-clean-") as directory:
        clean=Path(directory)
        zipfile.ZipFile(io.BytesIO(archive)).extractall(clean)
        env=dict(os.environ,PYTHONIOENCODING="utf-8",PYTEST_DISABLE_PLUGIN_AUTOLOAD="1")
        proc=subprocess.run([sys.executable,"experiments/great_pumpkin_001/verify.py"],
                            cwd=clean,env=env,capture_output=True,text=True,encoding="utf-8")
        print(proc.stdout,flush=True)
        if proc.returncode:
            raise RuntimeError(proc.stderr)
        evidence=json.loads((clean/"experiments/great_pumpkin_001/results/verification.json").read_text())
    out=HERE/"results/clean_export.json"
    out.write_text(json.dumps(dict(revision=revision,method="git archive to a new temporary directory",
                                  exit_code=proc.returncode,checks=evidence),indent=2)+"\n",
                   encoding="utf-8",newline="\n")
    print("Clean export passed:",revision)


if __name__=="__main__":
    main()

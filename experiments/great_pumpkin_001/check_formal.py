"""Fresh Lean kernel check plus explicit axiom and source audit."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess

ROOT=Path(__file__).resolve().parents[2]
PROJECT=ROOT/"formal/great_pumpkin"


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--lean",default="lean")
    args=ap.parse_args()
    lean=str(Path(args.lean).resolve()) if Path(args.lean).exists() else args.lean
    source=PROJECT/"GreatPumpkin.lean"
    text=source.read_text(encoding="utf-8")
    if re.search(r"\b(sorry|admit|axiom)\b",text):
        raise RuntimeError("unverified proof construct in source")
    version=subprocess.check_output([lean,"--version"],text=True,encoding="utf-8").strip()
    if "version 4.24.0," not in version:
        raise RuntimeError("pinned Lean 4.24.0 required: "+version)
    proc=subprocess.run([lean,"--trust=0","GreatPumpkin.lean"],cwd=PROJECT,
                        text=True,encoding="utf-8",capture_output=True,check=True)
    reports=re.findall(r"depends on axioms: \[(.*?)\]",proc.stdout)
    if len(reports)!=9:
        raise RuntimeError("expected nine axiom reports")
    for report in reports:
        if set(report.split(", "))-{"propext","Classical.choice","Quot.sound"}:
            raise RuntimeError("unexpected axiom: "+report)
    result=dict(version=version,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                command="lean --trust=0 GreatPumpkin.lean",exit_code=proc.returncode,
                stdout=proc.stdout,stderr=proc.stderr,axiom_reports_checked=len(reports))
    out=Path(__file__).parent/"results/formal.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(proc.stdout)
    print("Formal check passed.")


if __name__=="__main__":
    main()

# Expedition 001 experiments

Start with the [scientific report](../../docs/expeditions/great-pumpkin-001.md).
This directory is experimental; committed `results/` are reproducibility
evidence, not established physical predictions. `scratch/` is ignored.

From repository root, with Python 3.11+ and the pinned requirements installed:

```text
python experiments/great_pumpkin_001/verify.py --quick
python experiments/great_pumpkin_001/verify.py
python experiments/great_pumpkin_001/audit_baselines.py
python experiments/great_pumpkin_001/run.py
python experiments/great_pumpkin_001/plot.py
```

Full verification uses Bash. On Windows it prefers `C:/Program Files/Git/bin/bash.exe`.
`audit_baselines.py` retrieves pinned files using Git into temporary directories,
never switches this checkout, and needs the full historical commits. The other
commands work in an exported clean source tree. `run.py --samples N` controls
the deterministic random mass-box sample count; the committed evidence uses
64, seed 20261028. Sample ranges are not global bounds.

Files: `model.py` derives representation coefficients, solves the analytic
system and shoots inverse couplings downward with SciPy; `run.py` generates
the science experiments; `test_expedition.py` covers independent checks and
failure cases; `verify.py` invokes the actual repository validation script;
`check_formal.py` records a pinned fresh Lean invocation. Plotting is optional
for the numerical results and uses the saved JSON only.

Artifacts:

- `baselines.json`: pre-expedition code and exact historical point reproductions.
- `experiments.json`: central cases, numerical/input challenges, mass boxes,
  adversarial offsets and conditional observable illustrations.
- `verification.json`: tests and corrected full-reference scans.
- `formal.json`: compiler identity, source hash and axiom reports.
- `threshold_geometry.png` / `.svg`: standalone scientific figure.

Reproduction can change irrelevant environment strings and sub-rounding floats.
Inspect numerical tolerances and scientific conventions rather than demanding
byte identity of generated JSON. Source hashes identify the core experiment
and formal proof files; Git identifies the complete dependency context.

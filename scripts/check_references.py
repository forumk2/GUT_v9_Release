"""Fail on numerical reference drift; not a physics validation."""
import json
from pathlib import Path
import sys

REFERENCES = {
    "lr_minimal": (9.454, 16.177, .02205),
    "so10_10_126_thresholds": (9.484, 16.194, .02245),
    "lr_dparity_symmetric": (10.494, 15.358, .02313),
}


def check(directory):
    for name, reference in REFERENCES.items():
        best = json.loads((Path(directory)/name/"best_fit.json").read_text())["best"]
        if best is None:
            raise ValueError(f"{name}: no solution")
        for key, expect, tol in zip(["log10_MI","log10_MGUT","alpha_GUT"], reference,
                                   [.0006,.0006,.000006]):
            if not abs(best[key]-expect) <= tol:
                raise ValueError(f"{name}: {key}={best[key]} differs from {expect} beyond {tol}")
        if not best["max_frac"] < 1e-6:
            raise ValueError(f"{name}: residual too large")
    print("All three numerical references match within rounding tolerances.")


if __name__ == "__main__":
    check(sys.argv[1])

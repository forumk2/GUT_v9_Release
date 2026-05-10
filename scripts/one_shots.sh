
#!/usr/bin/env bash
set -euo pipefail
export PYTHONIOENCODING=utf-8

# A) Ease b3, harden b2’s, keep the “good” threshold pattern
python rge_scan_v5.py --two-loop --multiplets spectrum_balanced_v8.json   --b3LR -5.8 --b2L -1.6 --b2R -1.6 --bBL 6.5 --split-r 0.9   --dalpha1 0.001 --dalpha2 0.001 --dalpha3 -0.001

# B) Slightly stronger variant (may push MGUT further)
python rge_scan_v5.py --two-loop --multiplets spectrum_balanced_v8.json   --b3LR -5.8 --b2L -1.7 --b2R -1.7 --bBL 6.7 --split-r 1.0   --dalpha1 0.001 --dalpha2 0.001 --dalpha3 -0.001

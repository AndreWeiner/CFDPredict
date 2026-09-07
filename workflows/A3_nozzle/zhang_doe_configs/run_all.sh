#!/bin/bash
# Run all Zhang DOE configs sequentially via the headless workflow runner.
# Reihenfolge: erst Verifikations-Optimum A4B1C2D4, dann L16-Array.
# Ueberspringt Cases die bereits ein `command_finished` Marker haben.
set -eo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../../.." && pwd)
RUNS="$REPO/runs/zhang_doe"
mkdir -p "$RUNS"

# Optimum zuerst -> falls man fruehzeitig abbricht, hat man den wichtigsten Datenpunkt
ORDER=(A4B1C2D4 \
       A1B1C1D1 A1B2C2D2 A1B3C3D3 A1B4C4D4 \
       A2B1C2D3 A2B2C1D4 A2B3C4D1 A2B4C3D2 \
       A3B1C3D4 A3B2C4D3 A3B3C1D2 A3B4C2D1 \
       A4B1C4D2 A4B2C3D1 A4B3C2D4 A4B4C1D3)

for code in "${ORDER[@]}"; do
    series="$RUNS/$code"
    if [ -f "$series/command_finished" ]; then
        echo "SKIP $code (already finished)"
        continue
    fi
    echo "=== Zhang DOE $code ==="
    mkdir -p "$series"
    python "$REPO/tools/run_workflow.py" A3_nozzle "$HERE/$code.json" \
        --workdir "$series" \
        --no-tail \
        2>&1 | tail -5
    echo "=== $code DONE ==="
done
echo "All Zhang configs through."

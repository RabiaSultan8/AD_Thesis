#!/bin/bash

HADDOCK_DIR=~/haddock_work/HADDOCK/haddock_runs
LOG_DIR=~/haddock_work/HADDOCK/haddock_logs
HADDOCK_RUN="/home/zxkj006/anaconda3/bin/python /home/zxkj006/haddock_work/haddock2.5-2025-08/haddock/run_haddock.py"
mkdir -p ${LOG_DIR}

export HADDOCK=/home/zxkj006/haddock_work/haddock2.5-2025-08
export HADDOCKTOOLS=${HADDOCK}/tools
export PYTHONPATH=${HADDOCK}:${PYTHONPATH}

# ── Phase 1: sequential run1 creation for all 40 ────────────────────────────
echo "=== Creating run1 directories for all 40 pairs (sequential) ==="
for pair_dir in ${HADDOCK_DIR}/*/; do
    pair_name=$(basename ${pair_dir})
    log=${LOG_DIR}/${pair_name}.log

    cd ${pair_dir}
    rm -rf run1 *.job *.inp 2>/dev/null

    echo "[$(date '+%H:%M:%S')] SETUP: ${pair_name}"
    ${HADDOCK_RUN} >> ${log} 2>&1

    if [ -d "run1" ]; then
        echo "[$(date '+%H:%M:%S')] OK:    ${pair_name}"
    else
        echo "[$(date '+%H:%M:%S')] FAIL:  ${pair_name} — check ${log}"
    fi
done

n_ok=$(ls -d ${HADDOCK_DIR}/*/run1 2>/dev/null | wc -l)
echo ""
echo "Setup complete: ${n_ok}/40 run1 directories ready"
echo ""

if [ "${n_ok}" -eq 0 ]; then
    echo "No run1 directories created — aborting."
    exit 1
fi

# ── Build sorted list of ready pairs ────────────────────────────────────────
mapfile -t ALL_PAIRS < <(ls -d ${HADDOCK_DIR}/*/run1 2>/dev/null | sort)
TOTAL=${#ALL_PAIRS[@]}
HALF=$(( TOTAL / 2 ))

BATCH1=("${ALL_PAIRS[@]:0:${HALF}}")
BATCH2=("${ALL_PAIRS[@]:${HALF}}")

echo "Total ready: ${TOTAL}  |  Batch1: ${#BATCH1[@]}  |  Batch2: ${#BATCH2[@]}"
echo ""

# ── Docking function ─────────────────────────────────────────────────────────
dock_pair() {
    local run1_dir=$1
    local pair_name=$(basename $(dirname ${run1_dir}))
    local log=${LOG_DIR}/${pair_name}.log

    export HADDOCK=/home/zxkj006/haddock_work/haddock2.5-2025-08
    export HADDOCKTOOLS=${HADDOCK}/tools
    export PYTHONPATH=${HADDOCK}:${PYTHONPATH}
    local HADDOCK_RUN="/home/zxkj006/anaconda3/bin/python /home/zxkj006/haddock_work/haddock2.5-2025-08/haddock/run_haddock.py"

    echo "[$(date '+%H:%M:%S')] DOCKING START: ${pair_name}"
    cd ${run1_dir}
    ${HADDOCK_RUN} >> ${log} 2>&1
    echo "[$(date '+%H:%M:%S')] DOCKING DONE:  ${pair_name}"
}

export -f dock_pair
export LOG_DIR

# ── Batch 1: first 20 ────────────────────────────────────────────────────────
echo "=== BATCH 1: Docking pairs 1-${HALF} in parallel ==="
for run1_dir in "${BATCH1[@]}"; do
    dock_pair "${run1_dir}" &
done
wait
echo ""
echo "[$(date '+%H:%M:%S')] Batch 1 complete."
echo ""

# ── Batch 2: next 20 ─────────────────────────────────────────────────────────
echo "=== BATCH 2: Docking pairs $((HALF+1))-${TOTAL} in parallel ==="
for run1_dir in "${BATCH2[@]}"; do
    dock_pair "${run1_dir}" &
done
wait
echo ""
echo "[$(date '+%H:%M:%S')] Batch 2 complete."

echo ""
echo "============================================"
echo "All ${TOTAL} pairs completed!"
echo "Logs: ${LOG_DIR}"
echo "============================================"

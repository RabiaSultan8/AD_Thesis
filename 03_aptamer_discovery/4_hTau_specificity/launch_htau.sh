#!/bin/bash

RUNS_DIR=/home/zxkj006/haddock_work/hTau_dockings/htau_runs
LOG_DIR=/home/zxkj006/haddock_work/hTau_dockings/htau_logs
HADDOCK_RUN="/home/zxkj006/anaconda3/bin/python /home/zxkj006/haddock_work/haddock2.5-2025-08/haddock/run_haddock.py"

export HADDOCK=/home/zxkj006/haddock_work/haddock2.5-2025-08
export HADDOCKTOOLS=${HADDOCK}/tools
export PYTHONPATH=${HADDOCK}:${PYTHONPATH}

mkdir -p ${LOG_DIR}

# ── Step 1: Sequential setup ─────────────────────────────────────────────────
echo "=== Creating run1 directories (sequential) ==="
for pair_dir in ${RUNS_DIR}/*/; do
    pair_name=$(basename ${pair_dir})
    log=${LOG_DIR}/${pair_name}.log

    cd ${pair_dir}
    rm -rf run1 *.job *.inp 2>/dev/null

    echo "[$(date '+%H:%M:%S')] SETUP: ${pair_name}"
    ${HADDOCK_RUN} >> ${log} 2>&1

    if [ -d "run1" ]; then
        echo "[$(date '+%H:%M:%S')] OK:    ${pair_name}"
    else
        echo "[$(date '+%H:%M:%S')] FAIL:  ${pair_name}"
    fi
done

n_ok=$(ls -d ${RUNS_DIR}/*/run1 2>/dev/null | wc -l)
echo ""
echo "Setup complete: ${n_ok}/20 run1 directories ready"

if [ "${n_ok}" -eq 0 ]; then
    echo "No run1 directories created — aborting."
    exit 1
fi

# ── Step 2: Parallel docking — all 20 at once ────────────────────────────────
echo ""
echo "=== Docking all 20 in parallel ==="

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

for pair_dir in ${RUNS_DIR}/*/; do
    [ -d "${pair_dir}/run1" ] && dock_pair "${pair_dir}/run1" &
done

wait

echo ""
echo "============================================"
echo "All 20 hTau control dockings completed."
echo "============================================"

#!/usr/bin/env python3
import os
import glob
from pathlib import Path

# --- USER SETTINGS ---
RESULTS_DIR_BASE  = Path("./docking_results_wave1/")
CREATEPL_EXEC     = "createpl"
TOP_N_MODELS      = 10
MAX_PARALLEL_JOBS = 4

def models_already_exist(aptamer_dir):
    return aptamer_dir.is_dir() and len(list(aptamer_dir.glob("model_*.pdb"))) > 0

def main():
    Path("logs").mkdir(exist_ok=True)
    hdock_out_files = sorted(glob.glob(f"{RESULTS_DIR_BASE}/**/*.out", recursive=True))

    if not hdock_out_files:
        print(f"ERROR: No .out files found under {RESULTS_DIR_BASE}")
        return

    print(f"Found {len(hdock_out_files)} HDOCK .out files")

    commands_to_run = []
    skipped = 0

    for out_filepath_str in hdock_out_files:
        out_filepath = Path(out_filepath_str).resolve()  # absolute path
        variant_dir  = out_filepath.parent
        ligand_name  = out_filepath.stem.replace('_dock', '')
        aptamer_dir  = (variant_dir / ligand_name).resolve()

        if models_already_exist(aptamer_dir):
            skipped += 1
            continue

        aptamer_dir.mkdir(exist_ok=True)

        # KEY FIX: cd INTO the unique aptamer_dir, pass absolute path to .out
        # createpl generates model_*.pdb in current directory → goes straight to aptamer_dir
        # No mv needed. No race condition. Each job has its own isolated workspace!
        command = (
            f"cd {aptamer_dir} && "
            f"{CREATEPL_EXEC} "
            f"{out_filepath} "                          # absolute path to .out file
            f"{ligand_name}_top{TOP_N_MODELS}.pdb "
            f"-nmax {TOP_N_MODELS} -complex -models "
            f"2>> {Path('logs/createpl_errors.log').resolve()} && echo '.'"
        )
        commands_to_run.append(command)

    pending = len(commands_to_run)
    print(f"Skipped (already done) : {skipped}")
    print(f"Pending                : {pending}")

    if not commands_to_run:
        print("All createpl jobs already complete!")
        return

    with open("commands_createpl.txt", "w") as f:
        for cmd in commands_to_run: f.write(cmd + "\n")

    with open("run_parallel_createpl.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"echo 'Starting {pending} createpl jobs on {MAX_PARALLEL_JOBS} cores...'\n")
        f.write(
            f"cat commands_createpl.txt | "
            f"xargs -P {MAX_PARALLEL_JOBS} -I % bash -c '%' | "
            f"tqdm --total {pending} --desc 'createpl Progress' > /dev/null\n"
        )
        f.write("echo 'All createpl jobs complete.'\n")

    os.chmod("run_parallel_createpl.sh", 0o755)
    print(f"\nRun: bash run_parallel_createpl.sh")

if __name__ == "__main__":
    main()

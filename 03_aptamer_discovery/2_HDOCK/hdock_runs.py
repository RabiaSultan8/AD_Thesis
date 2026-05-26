#!/usr/bin/env python3
import os
from pathlib import Path

# --- USER SETTINGS ---
RECEPTORS_DIR    = Path("./Final_pTau_Variants/")
LIGANDS_DIR      = Path("./Final_Aptamers/")
RESULTS_DIR_BASE = Path("./docking_results_wave1/")

HDOCK_EXEC       = "hdock"
MAX_PARALLEL_JOBS = 64   

NUM_DECOYS       = 100      
ANGLE_SAMPLING   = 15       

def get_pdb_files(directory):
    return sorted(Path(directory).glob("*.pdb"))

def job_is_done(out_file):
    # CHANGED: Ensure the file isn't just an empty error file
    return out_file.exists() and out_file.stat().st_size > 1000

def main():
    receptors = get_pdb_files(RECEPTORS_DIR)
    ligands   = get_pdb_files(LIGANDS_DIR)
    
    if not receptors:
        print("ERROR: No receptors found! Check the 'Final_pTau_Variants' folder.")
        return
    if not ligands:
        print("ERROR: No ligands found! Check the 'Final_Aptamers' folder.")
        return

    total_jobs = len(receptors) * len(ligands)

    commands_to_run = []
    skipped = 0
    
    for receptor in receptors:
        variant_dir = RESULTS_DIR_BASE / receptor.stem
        variant_dir.mkdir(parents=True, exist_ok=True)

        for ligand in ligands:
            out_file = variant_dir / f"{ligand.stem}_dock.out"
            
            # If the file exists and is reasonably large, skip it
            if job_is_done(out_file): 
                skipped += 1
                continue

            # Command executes, and if successful, prints a dot for tqdm
            cmd = (
                f"cd {variant_dir.resolve()} && "
                f"{HDOCK_EXEC} "
                f"{receptor.resolve()} "
                f"{ligand.resolve()} "
                f"-out {out_file.name} "
                f"-angle {ANGLE_SAMPLING} "
                f"-num {NUM_DECOYS} > /dev/null 2>&1 && echo '.'" 
            )
            commands_to_run.append(cmd)

    pending_jobs = len(commands_to_run)

    if pending_jobs == 0:
        print(f"All {total_jobs} jobs already complete!")
        return

    with open("commands_hdock.txt", "w") as f:
        for cmd in commands_to_run: f.write(cmd + "\n")

    # Write the Bash script
    with open("run_parallel_hdock.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"echo 'Starting {pending_jobs} HDOCK jobs on {MAX_PARALLEL_JOBS} cores...'\n")
        f.write(
            f"cat commands_hdock.txt | "
            f"xargs -P {MAX_PARALLEL_JOBS} -I % bash -c '%' | "
            f"tqdm --total {pending_jobs} --desc 'HDOCK Progress' > /dev/null\n"
        )
        f.write("echo '\nAll HDOCK jobs complete.'\n")

    os.chmod("run_parallel_hdock.sh", 0o755)
    
    print(f"Total Jobs : {total_jobs}")
    print(f"Skipped    : {skipped}")
    print(f"Pending    : {pending_jobs}")
    print(f"Run 'bash run_parallel_hdock.sh' to start with the live progress bar!")

if __name__ == "__main__":
    main()

import os
import glob
import subprocess
import numpy as np
from Bio.PDB import PDBParser, Selection
import freesasa

PDB_DIR = "trRosettaRNA_Aptamers"
OUTPUT_DIR = "Final_Aptamers"
USALIGN_BIN = "USalign"
CLASH_CUTOFF = 0.5
TM_CUTOFF = 0.85

def check_model(pdb_file):
    parser = PDBParser(QUIET=True)
    try:
        # 1. Biopython parsing & clash check
        structure = parser.get_structure("RNA", pdb_file)
        atoms = list(Selection.unfold_entities(structure, 'A'))
        
        if len(atoms) < 10: 
            return False, "Invalid PDB: Too few atoms"
            
        coords = np.array([a.get_coord() for a in atoms])
        diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
        dist_sq = np.sum(diff**2, axis=-1)
        np.fill_diagonal(dist_sq, 100)
        
        min_dist = np.sqrt(np.min(dist_sq))
        if min_dist < CLASH_CUTOFF:
            return False, f"Clash: atoms overlapping at {min_dist:.2f} A"

        # 2. FreeSASA accessibility check
        sasa_struct = freesasa.Structure(pdb_file)
        sasa_result = freesasa.calc(sasa_struct)
        area = sasa_result.totalArea()
        
        if area < 500:
            return False, f"SASA too low: {area:.2f} A^2"
            
        return True, "Pass"
    except Exception as e:
        return False, f"Parsing Error: {str(e)}"

def run_usalign(pdb1, pdb2):
    cmd = [USALIGN_BIN, pdb1, pdb2, "-mol", "RNA"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        for line in result.stdout.split('\n'):
            if line.startswith("TM-score="):
                return float(line.split()[1])
        return 0.0
    except:
        return 0.0

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    pdb_files = sorted(glob.glob(f"{PDB_DIR}/*.pdb"))
    if not pdb_files:
        print(f"No PDB files found in {PDB_DIR}/")
        return

    print(f"Starting V2 3D QC Pipeline for {len(pdb_files)} models...\n")

    valid_models = []
    failure_reasons = {}

    # Step 1 & 2
    for pdb in pdb_files:
        passed, msg = check_model(pdb)
        if passed:
            valid_models.append(pdb)
        else:
            failure_reasons[msg] = failure_reasons.get(msg, 0) + 1

    print(f"Models passing Filter 1 & 2: {len(valid_models)} / {len(pdb_files)}")
    if failure_reasons:
        print("Failure breakdown:")
        for reason, count in failure_reasons.items():
            print(f"  - {reason}: {count} models")
    print("")

    if not valid_models:
        print("Stopping pipeline: No valid models to cluster.")
        return

    # Step 3
    print("Running Filter 3 (US-align 3D Clustering)...")
    final_representatives = []
    
    for i, pdb1 in enumerate(valid_models):
        is_redundant = False
        for rep in final_representatives:
            if run_usalign(pdb1, rep) >= TM_CUTOFF:
                is_redundant = True
                break
        
        if not is_redundant:
            final_representatives.append(pdb1)
            os.system(f"cp {pdb1} {OUTPUT_DIR}/")
            
        if (i+1) % 50 == 0:
            print(f"Processed {i+1} models... currently kept {len(final_representatives)}.")

    print(f"\n============================================================")
    print(f"Pipeline Complete!")
    print(f"Initial Models: {len(pdb_files)}")
    print(f"Final High-Quality, Diverse Models: {len(final_representatives)}")
    print(f"Ready for HDOCK. Files saved to: ./{OUTPUT_DIR}/")
    print(f"============================================================")

if __name__ == "__main__":
    main()

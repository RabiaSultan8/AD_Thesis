import mdtraj as md
import pandas as pd
from pathlib import Path

# ===== USER PARAMETERS (Updated for your naming scheme) =====

# 1. Directory containing all your PDB files
PDB_DIRECTORY = Path("./protein_structures/") 

# 2. List of the 12 key residue numbers you are investigating
#    *** IMPORTANT: FILL THIS LIST WITH YOUR 12 RESIDUE NUMBERS ***
PHOSPHO_SITE_RESIDUE_IDS = [
    181, 202, 205, 212, 217, 231, 262, 356, 396, 404, 409, 422 # <-- REPLACE WITH YOUR LIST
]

# 3. Filename identifiers for your models
HEALTHY_IDENTIFIER = "hTau"
TARGET_IDENTIFIERS = ["pTau", "mpTau"]

# 4. Name of the output CSV file
OUTPUT_CSV_FILE = "epitope_SASA_analysis.csv"

# Mapping of phosphorylated residue names (as they appear in PDB files)
# MDTraj uses standard names, but we still need to know which are phosphorylated.
PHOSPHO_RESIDUE_NAMES = ['TPO', 'SEP', 'PTR']

# ===== SCRIPT LOGIC (No need to edit below this line) =====

def analyze_structure(pdb_file):
    """Analyzes a single PDB file using MDTraj to get SASA and residue info."""
    try:
        # 1. Load the PDB file
        traj = md.load_pdb(pdb_file)
        
        # 2. Compute SASA for each atom and sum per-residue
        # mode='residue' correctly sums atomic SASA values for each residue
        sasa = md.shrake_rupley(traj, mode='residue')[0] # Result is in nm^2, we'll convert to Å^2

        # 3. Get topology information into a pandas DataFrame
        topo_df, _ = traj.topology.to_dataframe()
        
        # 4. We only need the residue-level information, so drop duplicates
        residue_df = topo_df.drop_duplicates(subset=['resSeq', 'resName', 'chainID'])
        
        # 5. Add SASA data (convert nm^2 to Å^2 by multiplying by 100)
        residue_df['sasa_A2'] = sasa * 100
        
        # 6. Add filename for tracking
        residue_df['filename'] = pdb_file.name
        
        return residue_df
        
    except Exception as e:
        print(f"Error processing {pdb_file.name} with MDTraj: {e}")
        return None

def main():
    print("--- Starting Epitope Mapping (MDTraj Engine) ---")
    
    if not PDB_DIRECTORY.is_dir():
        print(f"Error: PDB directory '{PDB_DIRECTORY}' not found.")
        return

    all_pdb_files = list(PDB_DIRECTORY.glob("*.pdb"))
    if not all_pdb_files:
        print(f"Error: No .pdb files found in '{PDB_DIRECTORY}'.")
        return
        
    print(f"Found {len(all_pdb_files)} PDB files to analyze.")
    
    # --- Step 1: Analyze all structures and collect data ---
    all_results_df = pd.DataFrame()
    for pdb_file in all_pdb_files:
        print(f"Processing: {pdb_file.name}")
        residue_data = analyze_structure(pdb_file)
        if residue_data is not None:
            all_results_df = pd.concat([all_results_df, residue_data], ignore_index=True)

    # --- Step 2: Filter for our key residues and classify files ---
    key_residues_df = all_results_df[all_results_df['resSeq'].isin(PHOSPHO_SITE_RESIDUE_IDS)].copy()

    def get_file_type(filename):
        if HEALTHY_IDENTIFIER in filename:
            return "healthy"
        for identifier in TARGET_IDENTIFIERS:
            if identifier in filename:
                return "target"
        return "unknown"

    key_residues_df['file_type'] = key_residues_df['filename'].apply(get_file_type)
    key_residues_df['is_phosphorylated'] = key_residues_df['resName'].isin(PHOSPHO_RESIDUE_NAMES)

    # --- Step 3: Calculate baseline SASA from healthy controls ---
    healthy_df = key_residues_df[key_residues_df['file_type'] == 'healthy']
    if healthy_df.empty:
        print("\nError: No 'healthy' files found. Check HEALTHY_IDENTIFIER and filenames.")
        return
    baseline_sasa = healthy_df.groupby('resSeq')['sasa_A2'].mean().to_dict()
    print("\nCalculated baseline SASA (Å^2) from healthy controls:")
    print(baseline_sasa)

    # --- Step 4: Calculate ΔSASA for all phosphorylated target models ---
    target_phospho_df = key_residues_df[
        (key_residues_df['file_type'] == 'target') & 
        (key_residues_df['is_phosphorylated'])
    ].copy()

    if target_phospho_df.empty:
        print("\nWarning: No phosphorylated residues found in any 'target' files.")
        return

    target_phospho_df['baseline_sasa_A2'] = target_phospho_df['resSeq'].map(baseline_sasa).fillna(0)
    target_phospho_df['delta_sasa_A2'] = target_phospho_df['sasa_A2'] - target_phospho_df['baseline_sasa_A2']
    
    # --- Step 5: Save the ranked results ---
    final_ranked_df = target_phospho_df.sort_values(by='delta_sasa_A2', ascending=False)
    
    # Select and rename columns for clarity
    output_columns = ['filename', 'resSeq', 'resName', 'sasa_A2', 'baseline_sasa_A2', 'delta_sasa_A2']
    final_ranked_df = final_ranked_df[output_columns].rename(columns={
        'resSeq': 'residue_id',
        'resName': 'residue_name'
    })

    final_ranked_df.to_csv(OUTPUT_CSV_FILE, index=False)
    
    print(f"\n--- Epitope Mapping Complete! ---")
    print(f"Results saved to: {OUTPUT_CSV_FILE}")
    print("\nTop 10 Most Promising Individual Epitopes (Highest ΔSASA):")
    print(final_ranked_df.head(10))

if __name__ == "__main__":
    main()
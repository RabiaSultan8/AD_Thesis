# Phase II-3: HADDOCK Flexible Refinement

This module executes the flexible refinement of the top 40 aptamer-pTau pairs identified during the HDOCK blind docking phase. HADDOCK 2.5 is utilized to capture induced-fit conformational adaptations at the binding interface, driven by Ambiguous Interaction Restraints (AIRs).

## Execution Pipeline
The scripts must be executed in the following order to prepare the structures, configure the HADDOCK environment for RNA, and extract the final thermodynamic rankings.

### 1. Structure & Restraint Preparation
* **`prep_ptau.py`**: Cleans the pTau variants for HADDOCK compatibility (sets Chain A, renames TPO to TOP, removes H3T atoms).
* **`prep_aptamers.py`**: Cleans the RNA aptamer models (sets Chain B, verifies O2' and OP1/OP2 nomenclature, adds END records).
* **`make_air_tbl.py`**: Converts the active/passive residue lists generated in the HDOCK phase into HADDOCK-formatted `.tbl` restraint files.

### 2. Environment Setup
* **`patch_master_runcns.py`**: *CRITICAL.* Modifies the HADDOCK `run.cns` master configuration file specifically for RNA-protein docking. It updates the topology/parameter files to `dna-rna-allatom-hj-opls-1.3`, sets dielectric constants (`epsilon_0` and `epsilon_1`) to 78.0 to reflect the highly charged aqueous environment of nucleic acids, and enables DNA/RNA restraints. 
* **`setup_haddock_runs.py`**: Automates the creation of isolated run directories for all 40 pairs, injecting the paired PDBs, the `.tbl` files, and the `run.param` file.
* **`launch_haddock_runs.sh`**: Automates the sequential initialization (run1 setup) and the massively parallel execution of the HADDOCK pipeline (rigid-body, semi-flexible, and explicit solvent refinement) across all 40 pairs.

### 3. Scoring & Ranking
* *(HADDOCK 2.5 is executed natively within each run directory).*
* **`score_water.py`**: Parses the `file.list` and `cluster.out` files from the final explicit solvent refinement stage (`it1/water`). It calculates the `Top4Avg` score (the average HADDOCK score of the top four members of each cluster) and outputs the final ranked shortlist.

## Output
The final output is `haddock_water_scores.csv`, containing the Top 20 aptamer-pTau pairs ranked by their `Top4Avg` cluster score. These 20 pairs advance to the final phospho-specificity gate against unphosphorylated hTau controls.

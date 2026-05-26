# Phase II-4: Thermodynamic Specificity Gate (hTau Controls)

This module executes the critical negative-control docking phase. The top 20 aptamer candidates identified against phosphorylated Tau (pTau) are docked against unphosphorylated healthy Tau (hTau) variants to isolate sequences with strict target specificity.

> **Methodological Note on Restraints (AIRs):**
> To ensure a perfectly controlled thermodynamic comparison, the exact same Ambiguous Interaction Restraint (`ambig.tbl`) files generated for the pTau dockings in 2_HDOCK phase were reused for these hTau negative controls. This forces the aptamers to target the exact same spatial epitope on the healthy protein, allowing us to accurately isolate the energetic contribution of the phosphate groups.

## Execution Pipeline

* **`setup_htau_runs.py`**: Generates isolated HADDOCK run directories for the top 20 aptamers against their corresponding hTau conformers, injecting the PDBs, `run.param` files, and the reused Phase 1 `.tbl` files.
* **`launch_htau.sh`**: Automates the sequential execution of HADDOCK Phase 0 (topology generation) and parallelizes Phase 1 (rigid-body), Phase 2 (semi-flexible), and Phase 3 (explicit solvent refinement) across available cores.
* **`score_htau_water.py`**: Parses the output of the explicit solvent refinement, calculates the `Top4Avg` cluster scores for the hTau complexes, and outputs `htau_water_scores.csv`.
* **`compute_specificity.py`**: The definitive thermodynamic filter. It cross-references the pTau and hTau water-refined scores to calculate $\Delta\Delta G_{\text{specificity}}$ (Score_hTau - Score_pTau). This isolates candidates whose binding is strictly driven by the phosphorylation modifications.

## Output
The final output is `haddock_specificity_scores.csv`, containing the final ranking of aptamer candidates prioritized by maximum $\Delta\Delta G$, identifying the precise diagnostic probes advancing to Phase III molecular dynamics validation.

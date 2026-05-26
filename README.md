# Targeting Pathogenic Phospho-Tau Epitopes and Associated Pathways: A Novel RNA Aptamer and Repurposed Drugs for Alzheimer's Disease

This repository houses the computational framework supporting the Master's Dissertation. The pipeline integrates systems biology, de novo RNA aptamer design, and atomistic molecular dynamics to develop precise diagnostic and therapeutic probes targeting the intrinsically disordered Tau protein (MAPT).

## Project Contributions & Workflow
* **Transcriptomic Target Validation:** Cross-dataset weighted meta-analysis of microarray (GSE138260) and RNA-seq (GSE118553) cohorts to isolate reproducible neuroinflammatory signatures and validate MAPT/GSK3B as primary targets.
* **IDP Conformational Sampling:** D-I-TASSER modeling and 500 ns GROMACS molecular dynamics (MD) relaxation to capture the conformational ensemble of intrinsically disordered Tau441.
* **Phosphorylation Mapping:** Systematic construction of 16 phosphorylated Tau (pTau) variants via CHARMM-GUI, followed by per-residue solvent-accessible surface area (SASA) profiling to map exposed diagnostic epitopes.
* **De Novo Aptamer Generation:** Computational design of a 5,000-sequence G-rich RNA library, funneled through sequential 1D/2D biophysical filtering (ViennaRNA) and 3D structural clustering (CD-HIT-EST, RNApdist, US-align).
* **Hierarchical Virtual Screening:** 3,904 rigid-body blind docking runs (HDOCK) evaluated via a custom epitope-centric composite score, advancing top candidates to explicit-solvent flexible refinement (HADDOCK 2.5).
* **Thermodynamic Specificity Gate:** Negative-control docking against unphosphorylated healthy Tau (hTau) to calculate empirical $\Delta\Delta$G, isolating strictly phospho-specific sequences.
* **Explicit-Solvent MD Validation:** 200 ns production trajectories coupled with MM/PBSA binding free energy calculations (ff14SB/RNA.OL3) to definitively quantify complex stability and validate the final diagnostic candidates.

## Repository Architecture

The codebase is modularized sequentially according to the methodology defined in the manuscript. Each directory contains a dedicated `readme.md` detailing the required inputs, execution environment, and output metrics for that specific phase.

```text
AD_Thesis/
├── 01_transcriptomics/       # Phase Ia: 11-script R pipeline for cross-dataset meta-analysis
├── 02_structure_prep/        # Phase Ib: GROMACS MDPs, SASA scripts, and automation
├── 03_aptamer_discovery/     # Phase II: Hierarchical screening funnel
│   ├── 1_aptamer_library/   # Aptamers generation and filtering. RNAfold, clustering, and 3D QC scripts
│   ├── 2_HDOCK/             # Parallel HDOCK blind docking and Composite Selection Scoring
│   ├── 3_HADDOCK/           # Explicit solvent RNA-protein HADDOCK flexible refinement
│   └── 4_hTau_specificity/  # Negative control docking and ΔΔG extraction
└── 04_md_validation/         # Phase III: 200ns MD automation, PBC correction, MM/PBSA
```

## Execution & Dependencies
Due to the integration of sequence-level bioinformatics and atomistic molecular dynamics, this pipeline requires a heavily distributed computational environment. Core dependencies include:

* **Systems Biology:** R 4.3.2+ (DESeq2, limma, clusterProfiler, ComplexHeatmap)
* **Structural Preparation:** GROMACS 2024.3+, MDTraj, FreeSASA, US-align
* **RNA Tools:** ViennaRNA (RNAfold, RNApdist), trRosettaRNA
* **Docking & Scoring:** HDOCK, HADDOCK 2.5 (patched for RNA topologies)
* **Thermodynamics:** MPI-enabled `gmx_MMPBSA`

Please refer to the localized `README.md` files within each module for precise execution commands.

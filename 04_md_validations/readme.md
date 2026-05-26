# Phase III: Molecular Dynamics & Thermodynamic Validation

This module executes the explicit-solvent molecular dynamics (MD) simulations and subsequent binding free energy calculations for the top 5 aptamer-pTau complexes generated in Phase II.

Because highly charged RNA-protein complexes are inherently prone to topological instability and box-scaling errors during initialization, this protocol employs a rigorous, multi-stage equilibration funnel prior to the 200 ns production phase.

## Execution Pipeline

* **`gromacs_automation.sh`**: The master execution script. It automates complex preparation (splitting chains, handling HADDOCK topological nomenclature like `TOP` -> `TPO`), solvation, ion addition, and the sequential GROMACS execution.
* **Simulation Parameter Files (`.mdp`)**:
    * `em.mdp`: Steepest descent energy minimization (max force < 1000.0 kJ/mol/nm).
    * `ions.mdp`: Single-point topology generation for `genion`.
    * `nvt.mdp`: 0.5 ns canonical ensemble equilibration (V-rescale thermostat, position restraints).
    * `npt_1.mdp`: 2 ns isothermal-isobaric equilibration Phase 1 (C-rescale barostat, position restraints) for gentle initial volume stabilization.
    * `npt_2.mdp`: 3 ns isothermal-isobaric equilibration Phase 2 (Parrinello-Rahman barostat, unrestrained) to establish true atmospheric pressure conditions.
    * `md.mdp`: 200 ns production trajectory execution.
* **`mmpbsa.in`**: Configuration for `gmx_MMPBSA`. Executes Poisson-Boltzmann surface area (MM/PBSA) end-state free energy calculations utilizing the `ff14SB` (protein) and `RNA.OL3` (nucleic acid) forcefields to quantify complex stability.
* **`md_analysis.sh`**: Post-processing automation. It generates an isolated `Analysis/` subdirectory and executes a strict periodic boundary condition (PBC) correction protocol (whole, cluster, center, fit) to remove artifactual global tumbling. It subsequently extracts fundamental kinematics (RMSD, RMSF, H-bond occupancy, Rg) and triggers the final MM/PBSA calculation via MPI.

## Usage
To process a target complex, execute the pipeline sequentially from the root of the designated run directory:

### 1. Execute the Simulation Protocol
```bash
./gromacs_automation.sh <complex_input.pdb> [optional_run_name]
```

### 2. Execute Kinematic & Thermodynamic Analysis
*(Ensure the production `md.xtc` and `md.tpr` files are successfully generated before running this step. Also ensure `mmpbsa.in` configuration file is present in the working directory).*
```bash
./md_analysis.sh
```

## Analytical Note
The trajectories generated here serve to validate the thermodynamic stability of the induced-fit conformations predicted by HADDOCK. The MM/PBSA calculations, while employing an implicit solvent approximation for the free energy extraction, provide the definitive energetic ranking that isolates the most viable clinical diagnostic candidates.

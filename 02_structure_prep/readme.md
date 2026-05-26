# Phase Ib: IDP-Aware Structural Preparation

This module contains the configuration files and automation scripts used to generate the accessible conformational ensemble of intrinsically disordered Tau441 and prepare the phosphorylated variant library.

## Files
* `minim.mdp`, `nvt.mdp`, `npt.mdp`, `md.mdp`: GROMACS parameters for energy minimization, equilibration, and the 500 ns production molecular dynamics run used for conformational sampling.
* `gmx_automation.sh`: Bash pipeline automating the execution of the GROMACS simulation and subsequent trajectory analysis (PBC correction, RMSD, RMSF, Rg, SASA, Clustering).
* `calculate_sasa.py`: MDTraj-based script calculating per-residue solvent-accessible surface area to evaluate phosphosite exposure.
* `plot_pTau_sasa.py`: Generates the publication-quality visualization of ΔSASA distributions.

## Methodological Note on Simulation Times
The `md.mdp` file provided here is configured for the 500 ns conformational space sampling of unphosphorylated Tau441. 

For the post-phosphorylation relaxation runs required to relieve local steric strain after CHARMM-GUI modification, the exact same thermodynamic parameters were applied, but the simulation times were adjusted as follows:
* **pTau variants (single/dual site):** 20 ns production run.
* **mpTau variants (hyperphosphorylated, 12 sites):** 50 ns production run.

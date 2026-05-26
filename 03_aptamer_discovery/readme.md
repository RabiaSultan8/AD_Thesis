# Phase II: RNA Aptamer Discovery & Hierarchical Virtual Screening

This module contains the computational framework for the de novo design, biophysical filtering, and structural screening of a G-rich RNA aptamer library targeted against pathogenic phosphorylated Tau (pTau) epitopes. 

Because Tau is an intrinsically disordered protein (IDP), the binding landscape is massive and highly flexible. To address this, the discovery pipeline employs a strict hierarchical funnel: starting with sequence-level thermodynamics, progressing to rigid-body blind docking, and culminating in explicit-solvent flexible refinement.

## Pipeline Architecture

The workflow is strictly sequential and distributed across three submodules:

### [1. `1_aptamer_library/`](./1_aptamer_library) (Generation & 3D QC)
* **Input:** Rational design seeds and structural constraints (Class A-D).
* **Process:** Generates a 5,000-sequence library, applies sequential 1D/2D biophysical filters (ViennaRNA MFE, loop sizes), removes redundancy (CD-HIT-EST, RNApdist), and performs 3D quality control on trRosettaRNA models (clash detection, FreeSASA, US-align).
* **Output:** A refined, structurally diverse library of 244 high-quality 3D RNA aptamer models.

### [2. `2_HDOCK/`](./2_HDOCK) (Blind Docking & Epitope Selection)
* **Input:** 244 aptamer models and 16 prepared pTau conformational variants.
* **Process:** Executes 3,904 rigid-body blind docking runs to rapidly sample the global binding landscape. Parses the resulting complexes using a custom, epitope-centric Composite Selection Score (CSS) that heavily weights interactions with the specific phosphorylated residues.
* **Output:** The Top 40 (aptamer, pTau) structural pairs and their corresponding Ambiguous Interaction Restraint (AIR) tables.

### [3. `3_HADDOCK/`](./3_HADDOCK) (Flexible Refinement)
* **Input:** The Top 40 rigid-body pairs and AIR tables.
* **Process:** Adapts the HADDOCK 2.5 thermodynamic parameters for highly charged RNA environments (modifying dielectric constants) and executes flexible refinement in explicit solvent to capture induced-fit interface adaptations.
* **Output:** The final Top 20 short-listed aptamer candidates ranked by their `Top4Avg` water-refined cluster scores, ready for rigorous specificity docking with healthy tau.

## Execution
To reproduce Phase II, navigate into the subdirectories sequentially (`1` -> `2` -> `3`) and follow the specific execution instructions detailed in their respective README files.

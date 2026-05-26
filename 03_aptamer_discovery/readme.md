# Phase II: RNA Aptamer Discovery & Hierarchical Virtual Screening

This module contains the computational framework for the de novo design, biophysical filtering, and structural screening of a G-rich RNA aptamer library targeted against pathogenic phosphorylated Tau (pTau) epitopes. 

Because Tau is an intrinsically disordered protein (IDP), the discovery pipeline employs a strict hierarchical funnel: starting with sequence-level thermodynamics, progressing to rigid-body blind docking, culminating in explicit-solvent flexible refinement, and finishing with a rigorous negative-control specificity gate.

## Pipeline Architecture

* **[1. `1_aptamer_library/`](./1_aptamer_library) (Generation & 3D QC)**
    * **Input:** Rational design seeds and structural constraints (Class A-D).
    * **Process:** Generates a 5,000-sequence library, applies 1D/2D biophysical filters, removes redundancy, and performs 3D quality control on trRosettaRNA models (clash detection, FreeSASA, US-align).
    * **Output:** A structurally diverse library of 244 high-quality 3D RNA aptamer models.
* **[2. `2_HDOCK/`](./2_HDOCK) (Blind Docking & Epitope Selection)**
    * **Input:** 244 aptamer models and 16 prepared pTau conformational variants.
    * **Process:** Executes 3,904 rigid-body blind docking runs. Parses complexes using a custom Composite Selection Score (CSS) that heavily weights interactions with specific phosphorylated residues.
    * **Output:** The Top 40 (aptamer, pTau) structural pairs and their corresponding Ambiguous Interaction Restraint (AIR) tables.
* **[3. `3_HADDOCK/`](./3_HADDOCK) (Flexible Refinement)**
    * **Input:** The Top 40 rigid-body pairs and AIR tables.
    * **Process:** Adapts HADDOCK 2.5 thermodynamic parameters for RNA environments and executes flexible refinement in explicit solvent to capture induced-fit adaptations.
    * **Output:** The Top 20 aptamer candidates ranked by `Top4Avg` water-refined cluster scores.
* **[4. `4_hTau_specificity/`](./4_hTau_specificity) (Thermodynamic Specificity Gate)**
    * **Input:** The Top 20 aptamer candidates and unphosphorylated hTau structural controls.
    * **Process:** Executes HADDOCK refinement against hTau and calculates $\Delta\Delta$G_specificity (Score_hTau - Score_pTau) to penalize candidates that bind the unphosphorylated backbone.
    * **Output:** The final prioritized diagnostic probes advancing to MD validation.

## Execution
To reproduce Phase II, navigate into the subdirectories sequentially (`01` -> `02` -> `03` -> `04`) and follow the execution instructions detailed in their respective README files.

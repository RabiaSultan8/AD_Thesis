# Phase II-2: HDOCK Blind Docking & Epitope-Centric Selection

This module executes the first wave of hierarchical docking. Because the precise binding geometry on the intrinsically disordered Tau surface is unknown a priori, blind docking is used to allow the aptamers to explore the entire accessible surface of the 16 pTau variants.

## Execution Pipeline
Execute these scripts sequentially to reproduce the 3,904 docking runs and the custom scoring pipeline:

### 1. Docking & Pose Generation
* **`hdock_runs.py`**: Automates 3,904 blind docking runs across the 244 aptamer candidates and 16 pTau targets utilizing multiprocessing.
* **`createpl_runs.py`**: Extracts the top 10 3D complex models from the HDOCK `.out` files using the `createpl` utility in isolated workspaces.

### 2. Custom Epitope-Centric Scoring
* **`phospho_sites.json`**: Dictionary mapping each pTau variant to its specific phosphorylated residue numbers.
* **`hdock_shortlisting.py`**: The core analytical engine. It parses all poses and calculates a custom `PoseScore` (weighting Base, Sugar, Backbone, and Guanine contacts against the phosphate atoms). It applies a two-tier selection strategy (Broad Binders vs. Specific Binders) using a Composite Selection Score (CSS) to rank the pairs.
* **`prepare_haddock_pairs_v2.py`**: Filters the global shortlist down to the final 40 highest-confidence (aptamer, pTau) pairs and verifies the successful generation of Ambiguous Interaction Restraints (AIR) tables for downstream flexible refinement.

## Output
This pipeline outputs a prioritized list of 40 (aptamer, pTau) pairs along with their respective `*__actives.txt` and `*__passives.txt` AIR files, ready to be passed to HADDOCK.

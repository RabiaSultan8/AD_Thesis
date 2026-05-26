#!/bin/bash
set -e

# ============================================================
# Usage: ./gromacs_automation.sh <complex_input.pdb> [run_name]
# Example: ./gromacs_automation.sh pTau_2.9__Aptamer_1152.pdb Apt1152_pTau2.9
# ============================================================

INPUT_PDB="${1:?ERROR: No input PDB provided. Usage: $0 <input.pdb> [run_name]}"
RUN_NAME="${2:-$(basename "$INPUT_PDB" .pdb)}"

echo "===== INPUT: $INPUT_PDB  |  RUN: $RUN_NAME ====="

# ============================================================
# STEP 0: Preprocessing
# ============================================================
echo "===== STEP 0: Preprocessing ====="

# Rename TOP->TPO and SOP->SEP if present (safety net for some HADDOCK outputs)
# Column-specific to avoid touching unrelated fields
sed -i 's/^\(.\{17\}\)TOP/\1TPO/' "$INPUT_PDB"
sed -i 's/^\(.\{17\}\)SOP/\1SEP/' "$INPUT_PDB"

# ============================================================
# STEP 1: Split protein (chain A) and RNA (chain B)
# ============================================================
echo "===== STEP 1: Split protein and RNA ====="

awk '$5=="A"' "$INPUT_PDB" > protein_raw.pdb
echo "TER" >> protein_raw.pdb
echo "END" >> protein_raw.pdb

awk '$5=="B"' "$INPUT_PDB" > rna_raw.pdb
echo "TER" >> rna_raw.pdb
echo "END" >> rna_raw.pdb

# ============================================================
# STEP 2: Fix RNA phosphate oxygen naming
# HADDOCK outputs OP1/OP2; CHARMM36-jul2022 RNA RTP expects O1P/O2P
# This applies to rna_raw.pdb only — protein TPO already has O1P/O2P from CHARMM-GUI
# ============================================================
echo "===== STEP 2: Fix RNA atom names ====="

sed -i 's/ OP1 / O1P /g' rna_raw.pdb
sed -i 's/ OP2 / O2P /g' rna_raw.pdb

# ============================================================
# STEP 3: pdb2gmx — generate topologies
# Confirmed interactive selections for CHARMM36-jul2022:
#   FF=8 (CHARMM36), Water=1 (TIP3P-CHARMM)
#   Protein: N-term=1 (NH3+), C-term=0 (COO-)  [MET1 option 0 fails]
#   RNA:     5'=4 (5TER), 3'=6 (3TER)
# ============================================================
echo "===== STEP 3: pdb2gmx ====="

echo -e "10\n1\n1\n0" | gmx pdb2gmx \
    -f protein_raw.pdb \
    -o protein.gro \
    -p protein.top \
    -i protein_posre.itp \
    -ter -ignh

echo -e "10\n1\n4\n6" | gmx pdb2gmx \
    -f rna_raw.pdb \
    -o rna.gro \
    -p rna.top \
    -i rna_posre.itp \
    -ter -ignh

# ============================================================
# STEP 4: Merge GRO files
# ============================================================
echo "===== STEP 4: Merge GRO ====="

PROT_ATOMS=$(sed -n '2p' protein.gro | tr -d ' ')
RNA_ATOMS=$(sed -n '2p' rna.gro | tr -d ' ')
TOTAL=$((PROT_ATOMS + RNA_ATOMS))

{
    echo "${RUN_NAME} complex"
    echo "$TOTAL"
    tail -n +3 protein.gro | head -n "$PROT_ATOMS"
    tail -n +3 rna.gro | head -n "$RNA_ATOMS"
    tail -n 1 protein.gro
} > complex.gro

echo "Merged GRO: $PROT_ATOMS protein atoms + $RNA_ATOMS RNA atoms = $TOTAL total"

# ============================================================
# STEP 5: Build clean modular topology
# ============================================================
echo "===== STEP 5: Build CLEAN topology ====="

# ------------------------------------------------------------
# 5.1 Generate protein.itp
# ------------------------------------------------------------
awk '
BEGIN {copy=0}

/\[ moleculetype \]/ {copy=1}

copy {
    if ($0 ~ /^\[ system \]/) exit

    # Remove includes
    if ($0 ~ /^#include/) next

    # Remove POSRES block
    if ($0 ~ /^#ifdef POSRES/) {skip=1; next}
    if ($0 ~ /^#endif/ && skip==1) {skip=0; next}
    if (skip==1) next

    print
}
' protein.top > protein.itp

# Validate
if ! grep -q "\[ moleculetype \]" protein.itp; then
    echo "ERROR: protein.itp generation failed"
    exit 1
fi

# ------------------------------------------------------------
# 5.2 Generate rna.itp
# ------------------------------------------------------------
awk '
BEGIN {copy=0}

/\[ moleculetype \]/ {copy=1}

copy {
    if ($0 ~ /^\[ system \]/) exit

    # Remove includes
    if ($0 ~ /^#include/) next

    # Remove POSRES block
    if ($0 ~ /^#ifdef POSRES/) {skip=1; next}
    if ($0 ~ /^#endif/ && skip==1) {skip=0; next}
    if (skip==1) next

    print
}
' rna.top > rna.itp

# Validate
if ! grep -q "\[ moleculetype \]" rna.itp; then
    echo "ERROR: rna.itp generation failed"
    exit 1
fi

# ------------------------------------------------------------
# 5.3 Extract molecule names
# ------------------------------------------------------------
PROT_MOL=$(awk '
/\[ moleculetype \]/ {
    while (getline) {
        if ($0 !~ /^;/ && NF > 0) {
            print $1
            exit
        }
    }
}
' protein.itp)

RNA_MOL=$(awk '
/\[ moleculetype \]/ {
    while (getline) {
        if ($0 !~ /^;/ && NF > 0) {
            print $1
            exit
        }
    }
}
' rna.itp)

echo "Detected Protein molecule: $PROT_MOL"
echo "Detected RNA molecule: $RNA_MOL"

# ------------------------------------------------------------
# 5.4 Build clean complex.top
# ------------------------------------------------------------
cat > complex.top << EOF
; ============================================================
; Protein-RNA Complex Topology
; ============================================================

; ---- Force field ----
#include "charmm36-jul2022.ff/forcefield.itp"

; ---- Molecule definitions ----
#include "protein.itp"
#ifdef POSRES
#include "protein_posre.itp"
#endif

#include "rna.itp"
#ifdef POSRES
#include "rna_posre.itp"
#endif

; ---- Water ----
#include "charmm36-jul2022.ff/tip3p.itp"

#ifdef POSRES_WATER
[ position_restraints ]
;  i funct  fcx  fcy  fcz
   1    1   1000 1000 1000
#endif

; ---- Ions ----
#include "charmm36-jul2022.ff/ions.itp"

; ============================================================
; System
; ============================================================

[ system ]
Protein-RNA Complex

[ molecules ]
$PROT_MOL     1
$RNA_MOL      1
EOF

# ------------------------------------------------------------
# 5.5 Final validation checks
# ------------------------------------------------------------
echo "----- VALIDATION -----"

if ! grep -q "$PROT_MOL" complex.top; then
    echo "ERROR: Protein molecule missing in topology"
    exit 1
fi

if ! grep -q "$RNA_MOL" complex.top; then
    echo "ERROR: RNA molecule missing in topology"
    exit 1
fi

if [ ! -s protein.itp ] || [ ! -s rna.itp ]; then
    echo "ERROR: One of the ITP files is empty"
    exit 1
fi

echo "Topology built successfully."

# ============================================================
# STEP 6: Define box and solvate
# ============================================================
echo "===== STEP 6: Box and Solvate ====="

gmx editconf \
    -f complex.gro \
    -o boxed.gro \
    -c -d 1.2 -bt dodecahedron

gmx solvate \
    -cp boxed.gro \
    -cs spc216.gro \
    -p complex.top \
    -o solv.gro

# ============================================================
# STEP 7: Add ions
# ============================================================
echo "===== STEP 7: Ions ====="

gmx grompp \
    -f ions.mdp \
    -c solv.gro \
    -p complex.top \
    -o ions.tpr

echo "SOL" | gmx genion \
    -s ions.tpr \
    -o solv_ions.gro \
    -p complex.top \
    -pname NA \
    -nname CL \
    -neutral -conc 0.15

# ============================================================
# STEP 8: Energy minimisation
# ============================================================
echo "===== STEP 8: Energy Minimisation ====="

gmx grompp \
    -f em.mdp \
    -c solv_ions.gro \
    -p complex.top \
    -o em.tpr

gmx mdrun -v -deffnm em


# Verify EM converged
if grep -q "Steepest Descents converged" em.log; then
    echo "EM converged"
else
    echo "WARNING: EM may not have converged"
fi

echo "Potential" | gmx energy -f em.edr -o em_potential.xvg

# ============================================================
# STEP 9: NVT equilibration (0.5 ns, position restrained)
# ============================================================
echo "===== STEP 9: NVT Equilibration ====="

gmx grompp \
    -f nvt.mdp \
    -c em.gro \
    -r em.gro \
    -p complex.top \
    -o nvt.tpr

if [ -f nvt.cpt ]; then
    gmx mdrun -v -deffnm nvt -cpi nvt.cpt -append -nb gpu -pme gpu -update gpu -bonded gpu -pin on -ntomp 8
else
    gmx mdrun -v -deffnm nvt -nb gpu -pme gpu -update gpu -bonded gpu -pin on -ntomp 8
fi

# ============================================================
# STEP 10: NPT equilibration Phase 1 (2 ns, position restrained, Berendsen)
# ============================================================
echo "===== STEP 10: NPT Equilibration Phase 1 (restrained) ====="

gmx grompp \
    -f npt_1.mdp \
    -c nvt.gro \
    -r nvt.gro \
    -t nvt.cpt \
    -p complex.top \
    -o npt1.tpr

if [ -f npt1.cpt ]; then
    gmx mdrun -v -deffnm npt1 -cpi npt1.cpt -append -nb gpu -pme gpu -update gpu -bonded gpu -pin on -ntomp 8
else
    gmx mdrun -v -deffnm npt1 -nb gpu -pme gpu -update gpu -bonded gpu -pin on -ntomp 8
fi

# ============================================================
# STEP 11: NPT equilibration Phase 2 (3 ns, unrestrained, Parrinello-Rahman)
# ============================================================
echo "===== STEP 11: NPT Equilibration Phase 2 (unrestrained) ====="

gmx grompp \
    -f npt_2.mdp \
    -c npt1.gro \
    -t npt1.cpt \
    -p complex.top \
    -o npt2.tpr

if [ -f npt2.cpt ]; then
    gmx mdrun -v -deffnm npt2 -cpi npt2.cpt -append -nb gpu -pme gpu -update gpu -bonded gpu -pin on -ntomp 8
else
    gmx mdrun -v -deffnm npt2 -nb gpu -pme gpu -update gpu -bonded gpu -pin on -ntomp 8
fi

# ============================================================
# STEP 12: Production MD (200 ns)
# ============================================================
echo "===== STEP 12: Production MD (200 ns) ====="

gmx grompp \
    -f md.mdp \
    -c npt2.gro \
    -t npt2.cpt \
    -p complex.top \
    -o md.tpr

if [ -f md.cpt ]; then
    gmx mdrun -v -deffnm md -cpi md.cpt -append -nb gpu -pme gpu -update gpu -bonded gpu -pin on -ntomp 8
else
    gmx mdrun -v -deffnm md -nb gpu -pme gpu -update gpu -bonded gpu -pin on -ntomp 8
fi

echo "===== DONE: ${RUN_NAME} ====="
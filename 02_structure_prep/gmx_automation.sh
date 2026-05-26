#!/bin/bash

# Exit on any error
set -e

# === USER VARIABLES ===
FOLDER_NAME=$(basename "$PWD")
INPUT_PDB="${FOLDER_NAME}.pdb"
FORCEFIELD_TGZ="charmm36-jul2022.ff.tgz"
FORCEFIELD_DIR="charmm36-jul2022.ff"
WATER_MODEL="spc216.gro"

if [ ! -d "$FORCEFIELD_DIR" ]; then
    tar -xzf "$FORCEFIELD_TGZ"
fi

echo -e "1\n1\n1\n0" | gmx pdb2gmx -f "$INPUT_PDB" -o processed.gro -ter -ignh

gmx editconf -f processed.gro -o boxed.gro -c -d 1.5 -bt dodecahedron

gmx solvate -cp boxed.gro -cs "$WATER_MODEL" -p topol.top -o solvated.gro

gmx grompp -f minim.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1

echo "13" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr

gmx mdrun -v -deffnm em

echo q | gmx make_ndx -f em.gro -o index.ndx

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

gmx mdrun -v -deffnm nvt -nb gpu -pme gpu -bonded gpu -update gpu -pin on -ntomp 32

echo "4" | gmx genrestr -f nvt.gro -o posre_BB.itp -n

if grep -q 'POSRES_BB' topol.top; then
    echo "[INFO] POSRES_BB block already exists in topol.top. Skipping insertion."
else
    awk '
        /#endif/ && !done {
            print $0
            print ""
            print "#ifdef POSRES_BB"
            print "#include \"posre_BB.itp\""
            print "#endif"
            done=1
            next
        }
        { print }
    ' topol.top > topol.tmp && mv topol.tmp topol.top
    echo "[INFO] POSRES_BB block successfully inserted into topol.top."
fi

gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

gmx mdrun -v -deffnm npt -nb gpu -pme gpu -bonded gpu -update gpu -pin on -ntomp 32

gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md.tpr

gmx mdrun -v -deffnm md -nb gpu -pme gpu -bonded gpu -update gpu -pin on -ntomp 32

echo "✅ Your Simulation is Completed and Starting Analysis..."

mkdir Analysis

cd Analysis

echo "===== Post Simulation Analysis ====="

echo "===== STEP 1: PBC Correction ====="
echo "System" | gmx trjconv -f ../md.xtc -s ../md.tpr -o md_pbcWhole.xtc -pbc whole
echo -e "Protein\nSystem" | gmx trjconv -f md_pbcWhole.xtc -s ../md.tpr -o md_pbcCluster.xtc -pbc cluster
echo -e "Protein\nSystem" | gmx trjconv -f md_pbcCluster.xtc -s ../md.tpr -o md_pbcMOLcenter.xtc -pbc mol -ur compact -center
echo -e "Protein\nProtein" | gmx trjconv -f md_pbcMOLcenter.xtc -s ../md.tpr -o md_fit.xtc -fit rot+trans

echo "===== STEP 2: RMSD Calculations ====="
echo -e "Backbone\nBackbone" | gmx rms \
    -s ../md.tpr \
    -f md_fit.xtc \
    -n index.ndx \
    -o rmsd_protein.xvg \
    -tu ns

echo "===== STEP 3: RMSF Calculations ====="
echo "Protein" | gmx rmsf \
    -s ../md.tpr \
    -f md_fit.xtc \
    -n index.ndx \
    -o rmsf_protein.xvg \
    -res

echo "===== STEP 4: Radius of Gyration ====="
echo "Protein" | gmx gyrate \
    -s ../md.tpr \
    -f md_fit.xtc \
    -n index.ndx \
    -o rg_protein.xvg \
    -tu ns

echo "===== STEP 5: H-Bond Occupancy Calculations ====="
echo -e "Protein\nProtein" | gmx hbond \
	    -s ../md.tpr \
	        -f md_fit.xtc \
		    -n index.ndx \
    -num hbond_count.xvg \
    -tu ns

echo "===== STEP 6: SASA Analysis ====="
echo -e "Protein" | gmx sasa \
	    -s ../md.tpr \
	        -f md_fit.xtc \
		    -n index.ndx \
    -o sasa_protein.xvg \
    -tu ns

echo "===== STEP 7: Cluster Analysis ====="
echo "Backbone" | gmx cluster -f md_fit.xtc -s ../md.tpr -method gromos -cutoff 0.3 -g cluster.log -b 10000

echo "✅ ALL DONE ✅"

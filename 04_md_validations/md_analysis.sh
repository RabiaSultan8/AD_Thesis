echo "===== Post Simulation Analysis ====="

mkdir Analysis

cd Analysis

echo "===== STEP 1: Create index file ====="
gmx make_ndx -f ../md.tpr -o index.ndx << EOF
1 | 12
name 20 Protein_RNA
q
EOF

echo "===== STEP 2: PBC Correction ====="
# Step 1: Make molecules whole
echo "System" | gmx trjconv \
	    -f ../md.xtc -s ../md.tpr -n index.ndx \
    -o md_whole.xtc \
    -pbc whole

# Step 2: Cluster complex — keeps aptamer with Tau
echo -e "Protein_RNA\nSystem" | gmx trjconv \
	    -f md_whole.xtc -s ../md.tpr -n index.ndx \
	        -o md_cluster.xtc \
    -pbc cluster

# Step 3: Center and compact
echo -e "Protein_RNA\nSystem" | gmx trjconv \
    -f md_cluster.xtc -s ../md.tpr -n index.ndx \
    -o md_center.xtc \
    -pbc mol -ur compact -center

# Step 4: Fit — removes global tumbling
echo -e "Protein_RNA\nProtein_RNA" | gmx trjconv \
    -f md_center.xtc -s ../md.tpr -n index.ndx \
    -o md_fit.xtc \
    -fit rot+trans
	
echo "===== STEP 3: RMSD Calculations ====="
echo -e "Backbone\nBackbone" | gmx rms \
    -s ../md.tpr \
    -f md_fit.xtc \
    -n index.ndx \
    -o rmsd_protein.xvg \
    -tu ns
	
echo -e "Backbone\nRNA" | gmx rms \
    -s ../md.tpr \
    -f md_fit.xtc \
    -n index.ndx \
    -o rmsd_rna.xvg \
    -tu ns
	
echo "===== STEP 4: RMSF Calculations ====="
echo "Protein" | gmx rmsf \
    -s ../md.tpr \
    -f md_fit.xtc \
    -n index.ndx \
    -o rmsf_protein.xvg \
    -res
	
echo "RNA" | gmx rmsf \
    -s ../md.tpr \
    -f md_fit.xtc \
    -n index.ndx \
    -o rmsf_rna.xvg \
    -res

echo "===== STEP 5: H-Bond Occupancy Calculations ====="
echo -e "RNA\nProtein" | gmx hbond \
	    -s ../md.tpr \
	        -f md_fit.xtc \
		    -n index.ndx \
    -num hbond_count.xvg \
    -tu ns
	
echo "===== STEP 6: Radius of Gyration ====="
echo "Protein_RNA" | gmx gyrate \
    -s ../md.tpr \
    -f md_fit.xtc \
    -n index.ndx \
    -o rg_complex.xvg \
    -tu ns

echo "===== STEP 7: MM/PBSA ====="
mpirun -np 8 gmx_MMPBSA -O \
  -i ../mmpbsa.in \
  -cs ../md.tpr \
  -ct md_fit.xtc \
  -ci index.ndx \
  -cg 1 12 \
  -cp ../complex.top \
  -o MMPBSA_results.dat \
  -eo MMPBSA_energy.csv \
  -nogui
  
echo "===== Analysis Completed ====="

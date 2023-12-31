TINKER Files:

amoebabio18_add.prm - force field with PCA and oCNF
PYP_F28oCNF_min.xyz - starting file for trajectory
PYP_F62oCNF_min.xyz - starting file for trajectory
PYP_F92oCNF_min.xyz - starting file for trajectory
PYP_F96oCNF_min.xyz - starting file for trajectory
TinkerMD_EF_xyz.py  - Python script to analyze electric fields from Tinker MD simulations
TinkerMD_ProteinFields.py - Alternative Python script to analyze electric fields from Tinker MD simulations
TinkerMD_HBondAnalysis.py - Python script to analyze H-bonding from Tinker MD simulations (NOTE: PBC_postproc.py was used in the end)
TinkerMD_xyz2pdb.py - Translates ARC file to PDB (required due to non-natural amino acids)
submit_TINKER_analysis.sh - example submission script

GROMACS Files:
amber99sb-ildn.ff.zip - force field with PCA and oCNF; also contains residuetypes.dat for GMX pdb2gmx
PYP_F28oCNF_md.gro - starting file for trajectory
PYP_F62oCNF_md.gro - starting file for trajectory
PYP_F92oCNF_md.gro - starting file for trajectory
PYP_F96oCNF_md.gro - starting file for trajectory
GromacsMD_calc_EF.py - Python script to analyze electric fields from Gromacs MD simulationshh
GromacsMD_hbond_extended_v5.py - Python script to analyze H-bonding from Gromacs MD simulations (NOTE: PBC_postproc.py was used in the end)
submit_Gromacs_EFs.sh - Example script for submission to cluster
submit_hbond_extended.sh - Example script for submission to cluster

combine_PDB_trajectories.py - script for combining PDB trajectory files of repetition runs
PDB_postproc.py - script for analysis of H-bonding from PDB trajectory files; works for GROMACS and TINKER trajectories

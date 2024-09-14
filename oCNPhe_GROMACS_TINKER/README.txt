Files for publication 
"Critical Evaluation of Polarizable and Nonpolarizable Force Fields for Proteins Using Experimentally Derived Nitrile Electric Fields"
Jacob M. Kirsh, Jared Bryce Weaver, Steven G. Boxer, and Jacek Kozuch
Journal of the American Chemical Society 2024 146 (10), 6983-6991
DOI: 10.1021/jacs.3c14775

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
GROMACS folders.zip - folders with augmented force field file etc.

combine_PDB_trajectories.py - script for combining PDB trajectory files of repetition runs
PDB_postproc.py - script for analysis of H-bonding from PDB trajectory files; works for GROMACS and TINKER trajectories

Files for publication
"Hydrogen bond blueshifts in nitrile vibrational spectra are dictated by hydrogen bond geometry and dynamics"
Jacob M. Kirsh, and Jacek Kozuch
ChemRxiv, 2024
DOI: 10.26434/chemrxiv-2024-rjm5p

Data_for_HBblueshifts.txt - data of AMOEBA-based fields on oTN CN group, and DFT-based CN wavenumbers and TDMs

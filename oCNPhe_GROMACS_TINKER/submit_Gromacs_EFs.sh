#!/bin/bash

#SBATCH --job-name=PYP_F96oCNF_0q                  # replace name
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mail-user=jkozuch@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=all
#SBATCH --nodes=1
#SBATCH --ntasks=1                              # replace with value for your job
#SBATCH --mem=1G                                # replace with value for your job
#SBATCH --time=24:00:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module reset
module load fosscuda/2019b
module load GROMACS/2020-fosscuda-2019b

gmx traj -s PYP_md.tpr -f PYP_md.trr -n probe.ndx -ox co_coords.xvg
gmx traj -s PYP_md.tpr -f PYP_md.trr -n probe.ndx -of co_forces.xvg
gmx traj -s PYP_md_0q.tpr -f PYP_md_0q.trr -n probe.ndx -of co_forces_0q.xvg

python calc_EF_GROMACS.py -a 902 903 -at C N -ch 0.2278 -0.3538 -t PYP_md

gmx trjconv -s PYP_md.tpr -f PYP_md.trr -o PYP_md_center.xtc -center -pbc mol -ur compact << EOF
1
0
EOF
                          
gmx trjconv -s PYP_md.tpr -f PYP_md_center.xtc -o PYP_md_fit.xtc -fit rot+trans << EOF
1
0
EOF



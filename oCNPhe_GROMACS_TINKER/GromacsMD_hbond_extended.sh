#!/bin/bash
#SBATCH --job-name=PYP_eval                     # replace name
#SBATCH --mail-user=jkozuch@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=all
#SBATCH --nodes=1
#SBATCH --ntasks=1                              # replace with value for your job
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G                        # replace with value for your job
#SBATCH --time=24:00:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job
#SBATCH --export=NONE


module purge
module load GROMACS/2020-fosscuda-2019b

gmx hbond -s PYP_md.tpr -f PYP_md_fit.xtc -num zz_Hbond_num.xvg -dist zz_Hbond_dist.xvg -hbn zz_Hbond_ndx.ndx -n z_index_hbond.ndx -r 0.35 -a 30 -shell 0.35 << EOF
19
20
19
EOF


python z_gmx_hbond_extended_v5.py -p PYP_md -d 0.35 -a 30 -accatm 903 -accref 902






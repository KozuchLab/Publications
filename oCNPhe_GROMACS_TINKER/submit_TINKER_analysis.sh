#!/bin/bash

#SBATCH --job-name=Eval_F96_1                  # replace name
#SBATCH --partition=main
#SBATCH --mail-user=jkozuch@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=all
#SBATCH --nodes=1
#SBATCH --ntasks=1                              # replace with value for your job
#SBATCH --mem=12G                                # replace with value for your job
#SBATCH --time=24:00:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module add fosscuda/2020a          

python TinkerMD_ProteinFields.py -t abs -arc PYP_F96oCNF_md -a 1431 1432 -at C N -p 1.4150 1.1580 -f 2500

python TinkerMD_HBondAnalysis.py -arc PYP_F96oCNF_md -a 1431 1432

python TinkerMD_xyz2pdb.py -arc PYP_F96oCNF_md -a 1432 -d 10 




#!/bin/bash

#SBATCH --job-name=ReVal_penam                    # replace name
#SBATCH --mail-user=jkozuch@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=all
#SBATCH --nodes=1
#SBATCH --ntasks=1                              # replace with value for your job
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G                        # replace with value for your job
#SBATCH --time=7-00:00:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job
#SBATCH --export=NONE

source ~/.bash_profile

module purge
module load fosscuda/2020a gaussian Anaconda3/2021.05         # replace with value for your job

python main.py

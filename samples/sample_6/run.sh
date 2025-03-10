#!/bin/bash

#SBATCH -t 200:00:00                    # Walltime
#SBATCH -N 1                           # Number of nodes
#SBATCH -n 42                          # Number of processor cores (i.e. tasks)
#SBATCH -J "parcs_lpo"               # Job name
#SBATCH -p oldq                      # Partition name
#SBATCH -o output.txt                  # Output file name it will be printed in the $HOME directory by default
#SBATCH -e error.txt                   # Error file name it will be printed in the $HOME directory by default

# Move into the submission directory -------------------------------------------------------------------------------
cd $SLURM_SUBMIT_DIR

# Main -------------------------------------------------------------------------------------------------------------

rm -rf initial_*
rm -rf child_*
python ../../mofMain.py --input cycle1.yaml --cpus 4


##stress-ng --cpu $SLURM_CPUS_ON_NODE --timeout 60s --metrics-brief


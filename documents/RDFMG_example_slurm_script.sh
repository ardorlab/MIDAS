#!/bin/bash

#SBATCH -t 30:00:00                    # Walltime
#SBATCH -N 1                           # Number of nodes
#SBATCH -n 8                           # Number of processor cores (i.e. tasks)
#SBATCH -J MOF                         # Job name
#SBATCH -p defq                        # Partition name
##SBATCH --mail-user=UNITYID@ncsu.edu  # Email address
##SBATCH --mail-type=BEGIN             # Receive email when job start
##SBATCH --mail-type=END               # Receive email when job end
##SBATCH --mail-type=FAIL              # Receive email when job fail
#SBATCH -o output.txt                  # Output file name it will be printed in the $HOME directory by default
#SBATCH -e error.txt                   # Error file name it will be printed in the $HOME directory by default

# Move into the submission directory -------------------------------------------------------------------------------
cd $SLURM_SUBMIT_DIR

# Main -------------------------------------------------------------------------------------------------------------
rm -r initial_*
rm -r child_*

python mofMain.py --input input_file.yaml --cpus 8

##stress-ng --cpu $SLURM_CPUS_ON_NODE --timeout 60s --metrics-brief

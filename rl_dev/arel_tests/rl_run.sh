#!/bin/bash

#SBATCH -t 24:00:00                    # Walltime
#SBATCH -N 1                           # Number of nodes
#SBATCH -n 2                          # Number of processor cores (i.e. tasks)
#SBATCH -J "RL_testing"               # Job name
#SBATCH -p defq                      # Partition name
##SBATCH --mail-user=UNITYID@ncsu.edu  # Email address
##SBATCH --mail-type=BEGIN             # Receive email when job start
##SBATCH --mail-type=END               # Receive email when job end
##SBATCH --mail-type=FAIL              # Receive email when job fail
#SBATCH -o output.txt                  # Output file name it will be printed in the $HOME directory by default
#SBATCH -e error.txt                   # Error file name it will be printed in the $HOME directory by default

# Move into the submission directory -------------------------------------------------------------------------------
cd $SLURM_SUBMIT_DIR

# Main -------------------------------------------------------------------------------------------------------------

export CMSBIN=/cm/shared/apps/ncsu/CasmoSimulate/bin
rm -rf initial_*
rm -rf child*
rm -rf solution_?_*
rm -rf solution_??_*
python /home/rbdagost/optimization/MOF_exercises2/mofMain.py --input /home/rbdagost/optimization/MOF_exercises2/rl_rtsT_egreedy.yaml --cpus 1






##stress-ng --cpu $SLURM_CPUS_ON_NODE --timeout 60s --metrics-brief


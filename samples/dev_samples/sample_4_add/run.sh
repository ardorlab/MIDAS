#!/bin/bash

#SBATCH -t 200:00:00                    # Walltime
#SBATCH -N 1                           # Number of nodes
#SBATCH -n 2                          # Number of processor cores (i.e. tasks)
#SBATCH -J "blurple"               # Job name
#SBATCH -p gradq                      # Partition name
##SBATCH --mail-user=UNITYID@ncsu.edu  # Email address
##SBATCH --mail-type=BEGIN             # Receive email when job start
##SBATCH --mail-type=END               # Receive email when job end
##SBATCH --mail-type=FAIL              # Receive email when job fail
#SBATCH -o blurpoutput.txt                  # Output file name it will be printed in the $HOME directory by default
#SBATCH -e blurperror.txt                   # Error file name it will be printed in the $HOME directory by default

# Move into the submission directory -------------------------------------------------------------------------------
cd $SLURM_SUBMIT_DIR
# Main -------------------------------------------------------------------------------------------------------------

export CMSBIN=/cm/shared/apps/ncsu/CasmoSimulate/bin

python run_solutions.py

#casmo4e -v u1.00.02 lat21.inp 
#casmo4e -v u1.00.02 lat26.inp
#casmo4e -v u1.00.02 lat31.inp
#casmo4e -v u1.00.02 lat26p16.inp
#casmo4e -v u1.00.02 lat26p20.inp
#casmo4e -v u1.00.02 lat26p24.inp
#casmo4e -v u1.00.02 lat31p8.inp
#casmo4e -v u1.00.02 lat31p12.inp
#casmo4e -v u1.00.02 lat31p16.inp
#casmo4e -v u1.00.02 lat31p20.inp
#casmo4e -v u1.00.02 lat31p24.inp
#casmo4e -v u1.00.02 gad12.inp
#casmo4e -v u1.00.02 gad24.inp
#casmo4e -v u1.00.02 ifba80.inp
#casmo4e -v u1.00.02 ifba120.inp
#simulate3 pyrex_case.inp
##stress-ng --cpu $SLURM_CPUS_ON_NODE --timeout 60s --metrics-brief


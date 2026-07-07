#!/bin/bash 
# 
#SBATCH -J "r193-sur"
#SBATCH -p newq
#SBATCH -t 300:00:00
#SBATCH -N 4
#SBATCH --ntasks-per-node=1  # More explicit
#SBATCH --cpus-per-task=50
#SBATCH --mem-per-cpu=4G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --nodelist=node043,node044,node045,node046
#conda activate midas-surrogate

conda run -n midas-surrogate python ../../midasmain.py --input runcase.yaml --cpus 50 > temp.log
##conda run -n midas-surrogate python ../midasmain.py --input testoldcase.yaml --cpus 25 > temp.log
 

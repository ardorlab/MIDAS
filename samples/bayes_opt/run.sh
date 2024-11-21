#!/bin/bash

#SBATCH -J "NuScale"
#SBATCH -p newq
#SBATCH -t 100:00:00
#SBATCH -N 1
#SBATCH -n 10
conda activate ptest

python3 /home/cahowar9/cahowar9/MIDAS/midasmain.py --input /home/cahowar9/cahowar9/MIDAS/samples/bayes_opt/bayes_test.yaml --cpus 10

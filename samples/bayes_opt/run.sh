#!/bin/bash

#SBATCH -J "NuScale"
#SBATCH -p newq
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH -n 5
python3 midasmain.py --input bayes_test.yaml --cpus 5
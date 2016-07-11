#!/bin/bash
#SBATCH --job-name="SOLIDS"
#SBATCH --partition=slims
#SBATCH --mail-user=roberto.leon@gmail.com
#SBATCH --mail-type=ALL

module load cplex

srun ./main 3000

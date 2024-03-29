#!/bin/bash

#$ -l nc=4
#$ -p -50
#$ -r yes
#$ -q large.q

#SBATCH -n 4
#SBATCH --nice=50
#SBATCH --requeue
#SBATCH -p node03-06
SLURM_RESTART_COUNT=2

# echo $@
# echo $CONDA_PREFIX
# $CONDA_PREFIX/bin/Rscript src/WTS4_DimReduc.R $@
echo $@
Rscript src/WTS4_DimReduc.R $@
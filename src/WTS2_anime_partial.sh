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

echo $@
echo $CONDA_PREFIX
# $CONDA_PREFIX/bin/convert -delay 50 output/WTS2/heatmap/normalize_1/all/$1/*1.png $2
$CONDA_PREFIX/bin/convert -delay 50 `ls -v output/WTS2/heatmap/normalize_1/all/$1/new/*1.png` $2
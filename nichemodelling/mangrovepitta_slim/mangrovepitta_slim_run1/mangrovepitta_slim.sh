#!/bin/bash
#
#SBATCH --job-name=mangrovepitta_slim
#SBATCH --output=mangrovepitta_slim.txt
#SBATCH --account=2016005
#SBATCH --mail-user=davidtan@unm.edu
#SBATCH --mail-type=ALL
#
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=normal
#
############################################

module load miniconda3

## start SLiM environment

source activate slim4_env

cd /carc/scratch/projects/andersen2016005/pitta/nichemodelling/mangrovepitta_slim/mangrovepitta_slim_run1/

slim /carc/scratch/projects/andersen2016005/pitta/nichemodelling/mangrovepitta_slim/mangrovepitta_slim_run1/mp_slim.txt

#!/bin/bash
#
#SBATCH --job-name=gstacks
#SBATCH --output=gstacks.txt
#SBATCH --account=2016005
#SBATCH --mail-user=davidtan@unm.edu
#SBATCH --mail-type=ALL
#
#SBATCH --ntasks=8
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=normal
#
############################################

cd /carc/scratch/projects/andersen2016005/pitta/ddRAD/3_Ref_Map

module load miniconda3

## start stacks2 environment

source activate stacks2_env

gstacks -I /carc/scratch/projects/andersen2016005/pitta/ddRAD/3_Ref_Map -M /carc/scratch/projects/andersen2016005/pitta/ddRAD/3_Ref_Map/pitta_popmap.txt -O /carc/scratch/projects/andersen2016005/pitta/ddRAD/3_Ref_Map -t 8 -S .BWA.sorted.bam

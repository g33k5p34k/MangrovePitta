#!/bin/bash
#
#SBATCH --job-name=populations_bwp_mp_r0_noMAF
#SBATCH --output=populations_bwp_mp_r0_noMAF.txt
#SBATCH --account=2016005
#SBATCH --mail-user=davidtan@unm.edu
#SBATCH --mail-type=ALL
#
#SBATCH --ntasks=8
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=normal
#
############################################

cd /carc/scratch/projects/andersen2016005/pitta/ddRAD/3_Ref_Map

module load miniconda3

## start stacks2 environment

source activate stacks2_env

populations -P /carc/scratch/projects/andersen2016005/pitta/ddRAD/3_Ref_Map/ -O /carc/scratch/projects/andersen2016005/pitta/ddRAD/4_Stacksout/BWP_MP_R0_noMAF -M /carc/scratch/projects/andersen2016005/pitta/ddRAD/3_Ref_Map/bwp_mp_popmap.txt -t 8 --write-random-snp --fstats --smooth-fstats --vcf --genepop --structure --plink --phylip-var-all

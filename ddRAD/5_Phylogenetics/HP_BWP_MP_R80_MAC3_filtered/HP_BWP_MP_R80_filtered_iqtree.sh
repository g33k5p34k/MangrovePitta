#!/bin/bash
#
#SBATCH --job-name=iqtree_pitta
#SBATCH --output=iqtree_pitta.txt
#SBATCH --account=2016005
#SBATCH --mail-user=davidtan@unm.edu
#SBATCH --mail-type=ALL
#
#SBATCH --ntasks=8
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=normal
#
############################################

cd /carc/scratch/projects/andersen2016005/pitta/ddRAD/5_Phylogenetics/

module load miniconda3

## start stacks2 environment

source activate stacks2_env

iqtree -s /carc/scratch/projects/andersen2016005/pitta/ddRAD/5_Phylogenetics/HP_BWP_MP_R80_MAC3_filtered/HP_BWP_MP_R80_filtered_bpp.phy \
-p /carc/scratch/projects/andersen2016005/pitta/ddRAD/5_Phylogenetics/HP_BWP_MP_R80_MAC3_filtered/populations.all.partitions.phylip \
-m MFP \
-B 1000 \
-keep-ident \
-bnni
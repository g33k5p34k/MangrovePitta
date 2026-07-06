#!/bin/bash
#
#SBATCH --job-name=bpp_R80_filtered_subset1
#SBATCH --output=bpp_R80_filtered_subset1.txt
#SBATCH --account=2016005
#SBATCH --mail-user=davidtan@unm.edu
#SBATCH --mail-type=ALL
#
#SBATCH --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=96000
#SBATCH --partition=bigmem-1TB
#
############################################

module load miniconda3

## start bpp environment

source activate bpp_env

cd /carc/scratch/projects/andersen2016005/pitta/ddRAD/5_Phylogenetics/HP_BWP_MP_R80_MAC3_noadmix_clean_subset1000_1/

~/bpp/src/bpp --cfile /carc/scratch/projects/andersen2016005/pitta/ddRAD/5_Phylogenetics/HP_BWP_MP_R80_MAC3_noadmix_clean_subset1000_1/HP_BWP_MP_R80_MAC3_phylip_bpp_filtered.ctl

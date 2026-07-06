#!/bin/bash
#
#SBATCH --job-name=samtools_cat_sortreads
#SBATCH --output=samtools_cat_sortreads.txt
#SBATCH --account=2016005
#SBATCH --mail-user=davidtan@unm.edu
#SBATCH --mail-type=ALL
#
#SBATCH --ntasks=8
#SBATCH --time=16:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=normal
#
############################################

cd /carc/scratch/projects/andersen2016005/pitta/ddRAD/2_Reference_Mapping

module load miniconda3

## start stacks2 environment

source activate stacks2_env

while read name; do 
    samtools cat $name.BWA.bam $name.rem1.BWA.bam $name.rem2.BWA.bam | samtools sort -o /carc/scratch/projects/andersen2016005/pitta/ddRAD/3_Ref_Map/$name.BWA.sorted.bam
done < specimenlist.txt

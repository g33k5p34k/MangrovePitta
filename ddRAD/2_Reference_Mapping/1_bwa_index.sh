#!/bin/bash
#
#SBATCH --job-name=bwa_index
#SBATCH --output=bwa_index_out.txt
#SBATCH --account=2016005
#SBATCH --mail-user=davidtan@unm.edu
#SBATCH --mail-type=ALL
#
#SBATCH --ntasks=8
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=normal
#
############################################

cd /carc/scratch/projects/andersen2016005/pitta/ddRAD/2_Reference_Mapping

module load miniconda3

## start stacks2 environment

source activate stacks2_env

bwa-mem2 index /carc/scratch/projects/andersen2016005/pitta/genomes/HP_cucullata_Genome/ragtag_output/ragtag.scaffold.fasta

#!/bin/bash
#
#SBATCH --job-name=ragtag_cucullata
#SBATCH --output=ragtag_cucullata_out.txt
#SBATCH --account=2016005
#SBATCH --mail-user=davidtan@unm.edu
#SBATCH --mail-type=ALL
#
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=normal
#
############################################

cd /carc/scratch/projects/andersen2016005/pitta/genomes/HP_cucullata_Genome/

module load miniconda3

## start stacks2 environment

source activate stacks2_env

ragtag.py scaffold /carc/scratch/projects/andersen2016005/pitta/genomes/LTManakin_Genome/GCF_009829145.1_bChiLan1.pri_genomic.fna /carc/scratch/projects/andersen2016005/pitta/genomes/HP_cucullata_Genome/M.Irestedt_pitta_sordida.scf.fasta -o /carc/scratch/projects/andersen2016005/pitta/genomes/HP_cucullata_Genome/ragtag_output

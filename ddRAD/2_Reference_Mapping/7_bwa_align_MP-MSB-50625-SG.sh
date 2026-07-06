#!/bin/bash
#
#SBATCH --job-name=bwa_ref_align_MP-MSB-50625-SG
#SBATCH --output=bwa_ref_align_MP-MSB-50625-SG.txt
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

bwa-mem2 mem /carc/scratch/projects/andersen2016005/pitta/genomes/HP_cucullata_Genome/ragtag_output/ragtag.scaffold.fasta /carc/scratch/projects/andersen2016005/pitta/ddRAD/1_Process_Radtags/MP-MSB-50625-SG.1.fq.gz /carc/scratch/projects/andersen2016005/pitta/ddRAD/1_Process_Radtags/MP-MSB-50625-SG.2.fq.gz -t 8 > MP-MSB-50625-SG.BWA.sam
samtools view -bS MP-MSB-50625-SG.BWA.sam > MP-MSB-50625-SG.BWA.bam

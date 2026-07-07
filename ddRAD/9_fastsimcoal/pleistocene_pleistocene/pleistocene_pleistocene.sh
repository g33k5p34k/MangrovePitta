#!/bin/bash
#
#SBATCH --job-name=pitta_fastsimcoal_pleisto-pleisto
#SBATCH --output=pitta_fastsimcoal_pleisto-pleisto.txt
#SBATCH --account=2016005
#
#SBATCH --ntasks=20
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=normal
#
############################################
## change directory to working directory

cd /carc/scratch/projects/andersen2016005/pitta/ddRAD/9_fastsimcoal/pleistocene_pleistocene/

module load miniconda3

## start demographic history reconstruction environment

source activate dadi_env

fsc27093 -t /carc/scratch/projects/andersen2016005/pitta/ddRAD/9_fastsimcoal/pleistocene_pleistocene/BWP_MP_R80_fastsimcoal.tpl -e /carc/scratch/projects/andersen2016005/pitta/ddRAD/9_fastsimcoal/pleistocene_pleistocene/BWP_MP_R80_fastsimcoal.est -m -0 -C 10 -n 1000000 -L 100 -s 0 -M -c 0 --multiSFS
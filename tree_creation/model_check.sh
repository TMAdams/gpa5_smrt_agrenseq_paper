#!/bin/bash
#SBATCH --job-name="model_check"
#SBATCH --partition=medium
#SBATCH --mem=4G
#SBATCH -o model_check_%j.out
#SBATCH -e model_check_%j.err

Path_to_R_script=$1
input_file=$2

source activate r_phylogeny
Rscript $Path_to_R_script --inp $input_file

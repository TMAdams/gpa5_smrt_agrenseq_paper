#!/bin/bash
#SBATCH --job-name="bootstrap"
#SBATCH --partition=long
#SBATCH --cpus-per-task=32
#SBATCH --mem=30G
#SBATCH -o bootstrap_tree_%j.out
#SBATCH -e bootstrap_tree_%j.err

path_to_R_script=$1

source activate r_phylogeny
Rscript $path_to_R_script --cores $SLURM_CPUS_PER_TASK

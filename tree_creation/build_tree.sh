#!/bin/bash
#SBATCH --job-name="build_tree"
#SBATCH --partition=medium
#SBATCH --mem=4G
#SBATCH -o build_tree_%j.out
#SBATCH -e build_tree_%j.err

path_to_R_script=$1
inp=$2
inv=$3
gamma=$4
outp=$5

source activate r_phylogeny
Rscript $path_to_R_script --inp $inp --inv $inv --gamma $gamma --outp $outp

#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mail-type=end
#SBATCH --mail-user=shady.kuster@colostate.edu
#SBATCH --output=nonnmtCHAF_pcor_2025-07-25

#activate R
source /home/skuster/anaconda3/etc/profile.d/conda.sh
conda activate rforme

#run R script
Rscript nonnmt_pcor_CHAF.R

#deactivate envs
conda deactivate


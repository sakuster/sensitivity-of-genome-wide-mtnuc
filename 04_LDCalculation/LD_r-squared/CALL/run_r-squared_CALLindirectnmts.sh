#!/bin/bash

#SBATCH --time=02:00:00
#SBATACH --cpus-per-task=24
#SBATCH --mail-type=end
#SBATCH --mail-user=shady.kuster@colostate.edu
#SBATCH --output=indirectnmt_rsquared_CALL_08-29-24

#activate R
source /home/skuster/anaconda3/etc/profile.d/conda.sh
conda activate rforme

#run R script
Rscript indirectnmt_rsquared_calculation_CALL.R

#deactivate envs
conda deactivate


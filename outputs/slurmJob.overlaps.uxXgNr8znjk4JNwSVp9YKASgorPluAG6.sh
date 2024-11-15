#!/bin/bash
#SBATCH --job-name=overlaps
#SBATCH --partition=componc_cpu
#SBATCH --ntasks=1       
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
##SBATCH --mem=24G       
#SBATCH --time=6:00:00       
#SBATCH --error=%J_%I.err    
#SBATCH --output=%J_%I.out
source ~/miniforge3/etc/profile.d/conda.sh
conda activate r-env
Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R

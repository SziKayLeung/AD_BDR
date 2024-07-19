#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=100:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

## ---------------------------
## Purpose: Download Mayo RNA-Seq dataset for alignment to Iso-Seq annotation
## 
## 12/07/2022: Download 
## ---------------------------


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
export SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_BDR


##-------------------------------------------------------------------------

# Download, trim and align Mayo dataset 
bash $SC_ROOT/2_Differential_Analysis/3b_Mayo/2_download_isca.sh; wait
bash $SC_ROOT/2_Differential_Analysis/3b_Mayo/3_trim_kallisto.sh; wait 

# Download, trim and align ROSMAP dataset 
bash $SC_ROOT/2_Differential_Analysis/3b_Rosmap/2_download_isca.sh; wait
bash $SC_ROOT/2_Differential_Analysis/3b_Rosmap/3_trim_kallisto.sh; wait 


##-------------------------------------------------------------------------

# Prepare files for TappAS 
# generate expression and phenotype file from aggregated RNA-Seq datasets 
# RNA-Seq dataset: BDR, Mayo, ROSMAP
Rscript $SC_ROOT/2_Differential_Analysis/3c_public_hybrid_tappas_input.R
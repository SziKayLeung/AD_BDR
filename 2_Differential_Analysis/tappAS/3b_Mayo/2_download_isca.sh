#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=300:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=2_download_mayo_isca.o
#SBATCH --error=2_download_mayo_isca.e

## ---------------------------
## Purpose: Download Mayo RNA-Seq dataset for alignment to Iso-Seq annotation
## 
## 12/07/2022: Download 
## ---------------------------


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
export SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR
source $SC_ROOT/2_Differential_Analysis/bdr_differential.config
source activate sqanti2_py3

##-------------------------------------------------------------------------
echo "#************************************* DOWNLOAD DATA"
python $SC_ROOT/2_Differential_Analysis/3b_Mayo/1_download_dataset.py $MAYO_RAW_DIR
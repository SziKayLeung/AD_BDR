#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=50:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-52%15
#SBATCH --output=4b_talon_label-%A_%a.o
#SBATCH --error=4b_talon_label-%A_%a.e


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_BDR
source $SC_ROOT/1_ONT_Pipeline/bdr_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh

sample=${ALL_SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}

##-------------------------------------------------------------------------
# run_talon_label <input_sam> <output_dir>
run_talon_label ${WKD_ROOT}/4_tclean/${sample}/${sample}_clean.sam ${WKD_ROOT}/5_talon


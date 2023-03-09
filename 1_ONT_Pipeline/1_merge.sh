#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=12 # specify number of processors per node
#SBATCH --mem=200G # specify bytes of memory to reserve
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=1_demux.o
#SBATCH --error=1_demux.e


# 10/11/2022: ADBDR targeted datasets, run porechop

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR
source $SC_ROOT/1_ONT_Pipeline/bdr_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh


##-------------------------------------------------------------------------

# 1) run_merge <raw_directory> <sample_output_name>
run_merge ${RAW_FASTQ_1} $RAW_ROOT_DIR/P0063_20221026_100859_pass_merged.fastq.gz
run_merge ${RAW_FASTQ_2} $RAW_ROOT_DIR/P0063_20221031_10859_pass_merged.fastq.gz
merged_raw_files=$(ls $RAW_ROOT_DIR/*fastq.gz)
echo ${merged_raw_files}
cat ${merged_raw_files} > $RAW_ROOT_DIR/BDR_all_pass_merged.fastq.gz


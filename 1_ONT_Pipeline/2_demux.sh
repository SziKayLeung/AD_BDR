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
#SBATCH --array=0-8418%50
#SBATCH --output=2_demux-%A_%a.o
#SBATCH --error=2_demux-%A_%a.e


# 10/11/2022: ADBDR targeted datasets, run porechop 
# 6502 fastq files in RAW_FASTQ_1 dir
# 1917 fastq files in RAW_FASTq_2_dir
# total: 8419 files

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR
source $SC_ROOT/1_ONT_Pipeline/bdr_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh


##-------------------------------------------------------------------------

raw_fastq1_files=($(ls ${RAW_FASTQ_1}/*fastq.gz))
raw_fastq2_files=($(ls ${RAW_FASTQ_2}/*fastq.gz))
raw_merged_fastq_files=($(echo ${raw_fastq1_files[@]} ${raw_fastq2_files[@]}))

SamplePath=${raw_merged_fastq_files[${SLURM_ARRAY_TASK_ID}]}
Sample=$(basename ${SamplePath} .fastq.gz)

echo "Processing ${Sample}"

# 3) run_porechop <raw.fastq.gz> <output_dir>
run_porechop ${SamplePath} ${WKD_ROOT}/1_demultiplex/${Sample} > ${WKD_ROOT}/1b_demultiplex_merged/log/${Sample}.log

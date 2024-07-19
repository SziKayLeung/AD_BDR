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
#SBATCH --output=../OutputADBDR/ONTBatch2/2log/2_demux_batch2.o
#SBATCH --error=../OutputADBDR/ONTBatch2/2log/2_demux_batch2.e


# 09/01/2023: ADBDR targeted datasets Batch 2
# 11054 fastq files in RAW_FASTQ_3 dir

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_BDR
source $SC_ROOT/1_ONT_Pipeline/bdr_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh


##-------------------------------------------------------------------------

raw_fastq3_files=($(ls ${RAW_FASTQ_3}/*fastq.gz))
#echo "${#raw_fastq3_files[@]}"

for SLURM_ARRAY_TASK_ID in {0..11053}; do 

echo $SLURM_ARRAY_TASK_ID
SamplePath=${raw_fastq3_files[${SLURM_ARRAY_TASK_ID}]}
Sample=$(basename ${SamplePath} .fastq.gz)

echo "Processing ${Sample}"

# 3) run_porechop <raw.fastq.gz> <output_dir>
run_porechop ${SamplePath} ${WKD_ROOT}/1_demultiplex/Batch2/${Sample} > ${WKD_ROOT}/1_demultiplex/Batch2/log/${Sample}.log

done
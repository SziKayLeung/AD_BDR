#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-52%15
#SBATCH --output=../Output/ONTBatch2/3log/3_cutadapt_minimap2_tclean-%A_%a.o
#SBATCH --error=../Output/ONTBatch2/3log/3_cutadapt_minimap2_tclean-%A_%a.e


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR
source $SC_ROOT/1_ONT_Pipeline/bdr_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh

sample=${ALL_SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}

##-------------------------------------------------------------------------

# merge each sample into one fastq file 
merge_fastq_across_samples ${sample} ${WKD_ROOT}/1_demultiplex/Batch2 ${WKD_ROOT}/1b_demultiplex_merged/Batch2

# delinate polyA and polyT sequences, reverse complement polyT sequences, remove polyA from all sequences
post_porechop_run_cutadapt ${WKD_ROOT}/1b_demultiplex_merged/Batch2/${sample}_merged.fastq ${WKD_ROOT}/2_cutadapt_merge/Batch2

# map combined fasta to reference genome
run_minimap2 ${WKD_ROOT}/2_cutadapt_merge/Batch2/${sample}_merged_combined.fasta ${WKD_ROOT}/3_minimap/Batch2

# run transcript clean on aligned reads
run_transcriptclean ${WKD_ROOT}/3_minimap/Batch2/${sample}_merged_combined_sorted.sam ${WKD_ROOT}/4_tclean/Batch2

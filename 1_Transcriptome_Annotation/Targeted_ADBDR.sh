#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --reservation=research_project-mrc148213_5
#SBATCH --array=0 # 1 sample
#SBATCH --output=Targeted_ADBDR_Part1.o
#SBATCH --error=Targeted_ADBDR_Part1.e

# 18/04/2021: Batch 1

#************************************* DEFINE GLOBAL VARIABLES
# File directories
cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR; mkdir IsoSeq Post_IsoSeq RNASeq
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/IsoSeq
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/RNASeq

cd $Isoseq3_WKD; mkdir CCS LIMA REFINE CLUSTER
cd $Isoseq3_WKD/LIMA; mkdir BATCHES
cd $Isoseq3_WKD/REFINE; mkdir BATCHES
cd $Isoseq3_WKD/CLUSTER; mkdir BATCHES
cd $PostIsoseq3_WKD; mkdir mkdir MAP TOFU SQANTI2 KALLISTO TAMA SQANTI_TAMA_FILTER
cd $RNASeq_WKD; mkdir MAPPED

# For Pooled Targeted
### Important order of BAM files is the same as the sample names
BATCH_NAMES=(Targeted_Seq_1)
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/ad_bdr/1_Transcriptome_Annotation
cat $FUNCTIONS/Isoseq_Targeted_ADBDR_Raw.txt
BAM_FILES=(`cat $FUNCTIONS/Isoseq_Targeted_ADBDR_Raw.txt | egrep -v "^\s*(#|$)"`)

# Other input files and directory
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
# sourcing functions script
source $FUNCTIONS/Targeted_ADBDR_Functions.sh

module load Miniconda2/4.3.21
################################################################################################
echo "#*************************************  Isoseq3 [Function 1, 2]"
BATCH=${BATCH_NAMES[${SLURM_ARRAY_TASK_ID}]}
BAM_FILE=${BAM_FILES[${SLURM_ARRAY_TASK_ID}]}
# Isoseq3.4.0
    # run_CCS <input_ccs_bam> <prefix_output_name> <Output_directory>
    # run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
    # run_REFINE $Sample $Input_LIMA_directory $Output_directory
    # run_targeted_REFINE $Input_Pooled_Sample $Input_config_file $Input_Lima_sample $Input_LIMA_directory $Output_directory
    # run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
run_CCS ${BAM_FILE} ${BATCH} $Isoseq3_WKD/CCS
run_LIMA ${BATCH} $Isoseq3_WKD/CCS $Isoseq3_WKD/LIMA "multiplex"
run_LIMA ${BATCH} $Isoseq3_WKD/CCS $Isoseq3_WKD/LIMA/BATCHES "no_multiplex"
run_REFINE ${BATCH} $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE/BATCHES
run_CLUSTER ${BATCH} $Isoseq3_WKD/REFINE/BATCHES $Isoseq3_WKD/CLUSTER/BATCHES

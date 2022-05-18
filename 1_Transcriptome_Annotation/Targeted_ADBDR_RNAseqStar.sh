#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=3:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=ADBDR_Rnaseq_Star-%A_%a.o
#SBATCH --error=ADBDR_Rnaseq_Star-%A_%a.e
#SBATCH --array=0-29%15 #30samples

#************************************* DEFINE LOCAL VARIABLES
RNASeq_Filtered=/lustre/projects/Research_Project-193356/Project_10202/11_fastp_trimmed
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/RNASeq

samples=(BBN006_30024 BBN_24548 BBN002_28350 BBN_20194 BBN_24260 BBN003_30195 BBN002_29416 BBN_24287 BBN002_30035 BBN_15250 BBN002_29087 BBN003_26927 BBN_24289 BBN006_29902 BBN_25890 BBN002_29471 BBN_15247 BBN_18399 BBN002_26311 BBN_18405 BBN_24938 BBN006_29162 BBN_24253 BBN_19616 BBN_4240 BBN002_30920 BBN_15237 BBN_19632 BBN004_26227 BBN006_26347)
sample=${samples[${SLURM_ARRAY_TASK_ID}]}

# sourcing functions script
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/AD_BDR/1_Transcriptome_Annotation
source $FUNCTIONS/Targeted_ADBDR_Functions.sh

module load Miniconda2/4.3.21

################################################################################################
echo "#************************************* RNAseq [Function 9, 10]"
## 9) run_star <list_of_samples> <rnaseq_input_directory> <output_dir>
run_star ${sample} $RNASeq_Filtered $RNASeq_WKD/MAPPED




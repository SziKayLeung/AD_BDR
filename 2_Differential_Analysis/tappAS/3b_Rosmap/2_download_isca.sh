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
#SBATCH --output=rosmap_download.o
#SBATCH --error=rosmap_download.e

# 22/06/2022: Download ROSMAP RNA-Seq dataset (4Tb)

#************************************* DEFINE GLOBAL VARIABLES
DiffAnalysis_WKD=/lustre/projects/Research_Project-MRC148213/lsl693/IsoSeq/Targeted_Transcriptome/ADBDR/Differential
#RNASeq_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/ROSMAP/RNASeq
FUNCTIONS=/lustre/projects/Research_Project-MRC148213/lsl693/Scripts/AD_BDR/2_Differential_Analysis/ROSMAP
RNASeq_DIR=/lustre/projects/Research_Project-MRC190311/ROSMAP/RNASeq

#************************************* DEFINE LOCAL VARIABLES
mkdir -p $DiffAnalysis_WKD/ROSMAP
mkdir -p $RNASeq_DIR/Synapse $RNASeq_DIR/Trimmed

module load Miniconda2
source activate sqanti2_py3

echo "#************************************* DOWNLOAD DATA"
raw_downloads=(`find $RNASeq_DIR/Synapse -maxdepth 1 -name "*fastq*"`)
python $FUNCTIONS/1_download_dataset.py $RNASeq_DIR/Synapse


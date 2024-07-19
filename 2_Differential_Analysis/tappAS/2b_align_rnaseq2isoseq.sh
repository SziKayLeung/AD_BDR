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

## ---------------------------
## Purpose: Align RNA-Seq (core 96 samples) to Iso-Seq annotation for differential expression analysis
## 
## 04/11/2021: Run Kallisto of the RNA-Seq on the Iso-Seq scaffold (using the --rf-stranded) 
## ---------------------------

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
export SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_BDR
source $SC_ROOT/1_IsoSeq_Pipeline/bdr_isoseq.config
source $SC_ROOT/2_Differential_Analysis/01_source_functions.sh

##-------------------------------------------------------------------------

# extract sample names from directory
cd $RNASEQ_FILTERED_DIR; SAMPLES_ALL=($(ls *BBN*))

# for each sample, take only the first two parts as the "name" (for downstream function)
# remove duplicated entries (given R1 and R2)
declare -a arr
for i in ${SAMPLES_ALL[@]}; do arr+=("$(echo $i |cut -d'_' -f2,3 )"); done
SAMPLE_NAMES=$(echo "${arr[@]}" | xargs -n1 | sort -u | xargs)


##-------------------------------------------------------------------------
echo "#************************************* RNA-Seq Expression Matrix on Iso-Seq scaffold"
# index Iso-Seq fasta file
cd $HYBRID_DIFF_DIR
kallisto index -i $NAME.idx $ISO_WKD_ROOT/$NAME.collapsed_classification.filtered_lite.fasta 2> $NAME.index.log

# run_kallisto_1sample <input_RNASEQ_rawdir> <sample> <input_ref_name_idx> <output_dir>
counter=1
for i in ${SAMPLE_NAMES[@]}; do
  echo $counter
  echo $i
  run_kallisto_1sample $RNASEQ_FILTERED_DIR $i $NAME.idx $HYBRID_DIFF_DIR
  counter=$((counter+1))
done

# generate_rnaseq_counts <input_dir>
# create one big expression matrix file from the kallisto output directory
generate_rnaseq_counts $HYBRID_DIFF_DIR


## Missing 5 samples from the 96
#MissingSamples=(BBN002_30035 BBN006_26347 BBN002_30920 BBN002.28350 BBN002_29471)
##Columns `BBN006_31492`, `BBN_18811`, `BBN004_28505`, and `BBN006_29228` don't exist. BBN006.31492
#MissingSamplesdonotexist=(BBN006.29228 BBN_18811 BBN006_31492 BBN006_29228)
#for i in ${MissingSamples[@]}; do echo $i; run_kallisto_1sample $RNASeq_Filtered $i AllBDRTargeted.idx $DiffAnalysis_WKD/RNASeq_SQANTI3; done
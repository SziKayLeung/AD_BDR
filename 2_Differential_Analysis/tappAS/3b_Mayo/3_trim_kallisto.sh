#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-92%10 # 93 samples
#SBATCH --output=3_trim_kallisto-%A_%a.o
#SBATCH --error=3_trim_kallisto-%A_%a.e

# 13/07/2022: Alignment of Mayo extreme braak stage samples (0,1,5,6)

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
module load Trim_Galore
export SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_BDR
source $SC_ROOT/2_Differential_Analysis/bdr_differential.config
source $SC_ROOT/2_Differential_Analysis/01_source_functions.sh
source activate sqanti2_py3


##-------------------------------------------------------------------------
echo "#************************************* PROCESS RAW DATA"

mkdir -p $MAYO_TRIMMED_DIR $MAYO_DIFF_DIR

# Create Kallist index file if not already created
kallisto_index_file=$HYBRID_DIFF_DIR/$NAME.idx
kallisto_fasta_file=$ISO_WKD_ROOT/$NAME.collapsed_classification.filtered_lite.fasta

if [ -f $kallisto_index_file ]; then
  echo "Kallisto index file created"
else
  kallisto index -i $kallisto_index_file $kallisto_fasta_file 2> AllBDRTargeted.index.log
fi

# Align RNA-Seq files individually
ALL_SAMPLES_NAMES=$(grep "^[^#;]" $MAYO_SAMPLE_NAMES_FILES | awk '{print $1}')
Sample=${ALL_SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}
echo ${Sample}
trim_and_run_kallisto $MAYO_RAW_DIR ${Sample} $kallisto_index_file $MAYO_TRIMMED_DIR $MAYO_DIFF_DIR
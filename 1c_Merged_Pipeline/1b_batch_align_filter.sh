#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=1:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-6 #7 samples (6 single nuclei sorted + 1 isoseq)
#SBATCH --output=1b_batch_align_filter-%A_%a.o
#SBATCH --error=1b_batch_align_filter-%A_%a.e

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_BDR
LOGEN_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen
source $SC_ROOT/1c_Merged_Pipeline/bdr_iso_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous


##-------------------------------------------------------------------------

export dir=${MERGED_ROOT}
mkdir -p ${dir}/1_align

cp $ISO_ROOT/5_merged_cluster/AllBDRTargeted.clustered.hq.fasta ${dir}/1_align/AllIsoTargeted.clustered.hq.fasta 
iso_fa=${dir}/1_align/AllIsoTargeted.clustered.hq.fasta

all_scn_fa=(${dir}/1_align/*SCN*)
all_iso_scn_fa=(${all_scn_fa[@]} ${iso_fa})

fasta=${all_iso_scn_fa[${SLURM_ARRAY_TASK_ID}]}
sample=$(basename ${fasta} | cut -d "." -f 1)

##-------------------------------------------------------------------------

# align
echo "Aligning ${sample}: ${fasta} ..."
echo "Output: ${dir}/1_align/${sample}_mapped.bam"
source activate isoseq3
cd ${dir}/1_align
pbmm2 align --preset ISOSEQ --sort ${GENOME_FASTA} ${fasta} ${sample}_mapped.bam --log-level TRACE --log-file ${sample}_mapped.log

# filter_alignment <input_name> <input_mapped_dir>
# output = ${sample}_mapped.filtered.bam, ${sample}_mapped.filtered.sorted.bam
filter_alignment ${sample}_mapped ${dir}/1_align
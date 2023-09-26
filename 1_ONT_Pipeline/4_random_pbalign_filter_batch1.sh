#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=4:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address


# Batch2: realign with pbmm2 align and filter alignment 

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR
source $SC_ROOT/1_ONT_Pipeline/bdr_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh

fasta=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/D_ONT/4_tclean/Batch1/BC01/BC01_clean.fa
sample=$(basename ${fasta} | cut -d "_" -f 1)
dir=${WKD_ROOT}/5_cupcake/5_align/Batch1/


##-------------------------------------------------------------------------

# align
echo "Aligning ${sample}: ${fasta} ..."
echo "Output: ${dir}/${sample}_mapped.bam"
source activate isoseq3
cd ${dir}
pbmm2 align --preset ISOSEQ --sort ${GENOME_FASTA} ${fasta} ${sample}_mapped.bam --log-level TRACE --log-file ${sample}_mapped.log

# filter_alignment <input_name> <input_mapped_dir>
# output = ${sample}_mapped.filtered.bam, ${sample}_mapped.filtered.sorted.bam
filter_alignment ${sample}_mapped ${dir}

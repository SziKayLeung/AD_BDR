#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=0:45:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=1a_prepare_samples.o
#SBATCH --error=1a_prepare_samples.e

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
source $SC_ROOT/1c_Merged_Pipeline/bdr_iso_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous

##-------------------------------------------------------------------------

source activate nanopore
#dir=$ONT_ROOT/5_cupcake
# copy and replace filename in ONT transcript clean directory
#echo "Replacing filenames in ONT transcript clean directory"
#replace_filenames_with_csv.py --copy --ext=filtered.bam -i=$ONT_ROOT/5_cupcake/5_align/Batch1 -f=$META_ROOT/ADBDR_Batch1_rename.csv -d=${dir}/5_align/combined &> ${dir}/5_align/combined/B1_copy.log
#replace_filenames_with_csv.py --copy --ext=filtered.bam -i=$ONT_ROOT/5_cupcake/5_align/Batch2 -f=$META_ROOT/ADBDR_Batch2_rename.csv -d=${dir}/5_align/combined &> ${dir}/5_align/combined/B2_copy.log
#replace_filenames_with_csv.py --copy --ext=filtered.fa -i=$ONT_ROOT/5_cupcake/5_align/Batch1 -f=$META_ROOT/ADBDR_Batch1_rename.csv -d=${dir}/5_align/combined_fasta &> ${dir}/5_align/combined_fasta/B1_copy.log
#replace_filenames_with_csv.py --copy --ext=filtered.fa -i=$ONT_ROOT/5_cupcake/5_align/Batch2 -f=$META_ROOT/ADBDR_Batch2_rename.csv -d=${dir}/5_align/combined_fasta &> ${dir}/5_align/combined_fasta/B2_copy.log

# sorted nuclei, convert transcriptclean bam to fasta
export dir=${MERGED_ROOT}
mkdir -p ${dir}/1_align

snuclei_bam=($(ls ${SNUCLEI_ROOT}/*Transcript_Clean_aligned_clean_merged.bam)) 
ls ${snuclei_bam}

subName=(SCN04 SCN05 SCN06 SCN07 SCN08 SCN09)
for i in {0..5}; do 
 
  echo "Processing: ${snuclei_bam[i]}"
  echo "Output: ${subName[i]}.fa"
  
  samtools bam2fq ${snuclei_bam[i]} | seqtk seq -A > ${dir}/1_align/${subName[i]}.fa
done
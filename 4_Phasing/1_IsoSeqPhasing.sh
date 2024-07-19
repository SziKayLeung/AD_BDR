#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=1_IsoSeqPhasing.o
#SBATCH --error=1_IsoSeqPhasing.e

# 19/06/2024: SNP phasing in AD-BDR Iso-Seq using Iso-Phase

module load Miniconda2
source activate sqanti2_py3 

# convert bam to fastq 
#source activate nanopore
#for i in *bam; do 
# 	echo $i
#	sample=$(basename "$i" | cut -d "." -f 1 )
#	echo $sample
#	samtools bam2fq $i > $sample.flnc.fastq
#done
#cat *flnc.fastq > merged.flnc.fastq

softwareDir=/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake/phasing/utils
humanReferenceFasta=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
cupcakeDir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/7_tofu
refineDir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/3_refine

#awk '{print $1 "\t" $5}' ${cupcakeDir}/AllBDRTargeted.collapsed.read_stat.txt > ${cupcakeDir}/AllBDRTargeted.collapsed.read_stat_renamed.txt
#awk -F " " '{ if ( $5 != "NA" ) print $0 }' ${cupcakeDir}/AllBDRTargeted.collapsed.read_stat.txt > ${cupcakeDir}/AllBDRTargeted.collapsed.read_stat_renamed.txt

cd /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing
python ${softwareDir}/select_loci_to_phase.py \
${humanReferenceFasta} ${refineDir}/merged.flnc.fastq \
${cupcakeDir}/AllBDRTargeted.collapsed.gff \
${cupcakeDir}/AllBDRTargeted.collapsed.read_stat_renamed.txt \
-c 40 --fq
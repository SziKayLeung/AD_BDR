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
#SBATCH --output=4_lorals.o
#SBATCH --error=4_lorals.e

module load Miniconda2 
source activate nanopore 

#cd /lustre/projects/Research_Project-MRC148213/lsl693/references/annotation
#gtf2bed --gtf gencode.v40.annotation.gtf --bed gencode.v40.annotation.bed

export PATH="/lustre/projects/Research_Project-MRC148213/lsl693/software/lorals:$PATH"

cd /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/6_minimap/
#samtools view -bS AllBDRTargeted.sam | samtools sort -o AllBDRTargeted.sorted.bam
#samtools index AllBDRTargeted.sorted.bam AllBDRTargeted.sorted.bam.bai

TRANSCRIPTOME_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/gencode.v40.transcripts.fa
BDR_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/5_merged_cluster/AllBDRTargeted.clustered.hq.fastq
minimap2 -t 32 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $TRANSCRIPTOME_FASTA $BDR_FASTA > AllBDRTargeted2Transcriptome.sam 2> AllBDRTargeted2Transcriptome.sam.log
samtools view -bS AllBDRTargeted2Transcriptome.sam | samtools sort -o AllBDRTargeted2Transcriptome.sorted.bam
samtools index AllBDRTargeted2Transcriptome.sorted.bam AllBDRTargeted2Transcriptome.sorted.bam.bai


# for lorals, need sorted bam file
cd /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/lorals
#calc_ase -b /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/6_minimap/AllBDRTargeted.sorted.bam -f /lustre/projects/Research_Project-MRC148213/lsl693#/AD_BDR/phasing/IsoSeq_IsoPhase.vcf

#mv lorals_out/ase.tsv .
#rmdir lorals_out
#annotate_ase -i ase.tsv -b /lustre/projects/Research_Project-MRC148213/lsl693/references/annotation/gencode.v40.annotation.bed -o ./ase_annotated.tsv

calc_asts -m quant -b /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/6_minimap/AllBDRTargeted.sorted.bam -x /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/6_minimap/AllBDRTargeted2Transcriptome.sorted.bam -i ase.tsv -o ./asts.tsv

# note error with process_asts
# https://github.com/LappalainenLab/lorals/issues/9
# had to setup install again
process_asts -i asts.tsv_asts_quant_mod.tsv -g /lustre/projects/Research_Project-MRC148213/lsl693/references/human/gencode.v40.annotation.gene_transcript_ids.txt -o .

GENCODE_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
chr17_76471122
m54082_210419_161354/66781811
samtools mpileup --output-QNAME -f ${GENCODE_FASTA} -r chr17:76471121-76471123 /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/6_minimap/AllBDRTargeted.sorted.bam > chr17_7647112_mpileup.txt
grep 76471122 chr17_7647112_mpileup.txt | awk '{ print $7 }' | tr ',' '\n'  > chr17_7647112_refAlleleTranscripts.txt

cupcakeCollpased=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/7_tofu/AllBDRTargeted.collapsed.group.txt
cupcakeCollpasedreadStat=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/7_tofu/AllBDRTargeted.collapsed.read_stat.txt
grep -w -f chr17_7647112_refAlleleTranscripts.txt ${cupcakeCollpased} | awk '{ print $1 }' > chr17_7647112_refAlleleTranscripts_PB.txt
grep -w -f chr17_7647112_refAlleleTranscripts_PB.txt ${cupcakeCollpasedreadStat} | awk '{ print $1 }' > chr17_7647112_refAlleleTranscripts_ccs.txt

sed -i 's/\/ccs//g' chr17_7647112_refAlleleTranscripts_ccs.txt
grep -f chr17_7647112_refAlleleTranscripts_ccs.txt merged.phased.partial.tagged.sam > chr17_7647112_refAlleleTranscripts_ccs.phased.partial.tagged.sam

samtools view -S -b chr17_7647112_refAlleleTranscripts_ccs.phased.partial.tagged.sam > chr17_7647112_refAlleleTranscripts_ccs.phased.partial.tagged.bam
samtools index chr17_7647112_refAlleleTranscripts_ccs.phased.partial.tagged.bam chr17_7647112_refAlleleTranscripts_ccs.phased.partial.tagged.bam.bai
grep ATGAGACATCGTCCCTACA chr17_7647112_refAlleleTranscripts_ccs.phased.partial.tagged.sam > chr17_7647112_refAlleleTranscripts_ccs_C.phased.partial.tagged.sam
grep ATGAGACATCGTCCGTACA chr17_7647112_refAlleleTranscripts_ccs.phased.partial.tagged.sam > chr17_7647112_refAlleleTranscripts_ccs_G.phased.partial.tagged.sam
awk '{print $1'} chr17_7647112_refAlleleTranscripts_ccs_G.phased.partial.tagged.sam > chr17_7647112_refAlleleTranscripts_ccs_G.txt
awk '{print $1'} chr17_7647112_refAlleleTranscripts_ccs_C.phased.partial.tagged.sam > chr17_7647112_refAlleleTranscripts_ccs_C.txt


export PATH="/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/miscellaneous:$PATH"
SQANTIGTF=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/9_sqanti3/basic/AllBDRTargeted.collapsed_corrected.gtf
subset_fasta_gtf.py ${SQANTIGTF} --gtf -i chr17_7647112_refAlleleTranscripts_PB.txt -o chr17_7647112_refAlleleTranscripts_AllBDRTargeted.collapsed_corrected -d . 
mv AllBDRTargeted.collapsed_corrected_chr17_7647112_refAlleleTranscripts_AllBDRTargeted.collapsed_corrected.gtf chr17_7647112_refAlleleTranscripts_AllBDRTargeted.collapsed_corrected.gtf

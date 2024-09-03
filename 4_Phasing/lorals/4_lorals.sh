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


##-------------------------------------------------------------------------

module load Miniconda2 
source activate nanopore 

export PATH="/lustre/projects/Research_Project-MRC148213/lsl693/software/lorals:$PATH"
export PATH="/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/miscellaneous:$PATH"

GENCODE_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
TRANSCRIPTOME_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/gencode.v40.transcripts.fa
BDR_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/5_merged_cluster/AllBDRTargeted.clustered.hq.fastq
cupcakeCollpased=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/7_tofu/AllBDRTargeted.collapsed.group.txt
cupcakeCollpasedreadStat=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/7_tofu/AllBDRTargeted.collapsed.read_stat.txt
SQANTIGTF=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/9_sqanti3/basic/AllBDRTargeted.collapsed_corrected.gtf
REFERENCE=/lustre/projects/Research_Project-MRC148213/lsl693/references

##-------------------------------------------------------------------------

# need annotation bed file for annotate_ase
gtf2bed --gtf ${REFERENCE}/annotation/gencode.v40.annotation.gtf --bed ${REFERENCE}/annotation/gencode.v40.annotation.bed

# generate sorted bam files 
cd /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/6_minimap/
samtools view -bS AllBDRTargeted.sam | samtools sort -o AllBDRTargeted.sorted.bam
samtools index AllBDRTargeted.sorted.bam AllBDRTargeted.sorted.bam.bai

# align to the reference transcriptome
minimap2 -t 32 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $TRANSCRIPTOME_FASTA $BDR_FASTA > AllBDRTargeted2Transcriptome.sam 2> AllBDRTargeted2Transcriptome.sam.log
samtools view -bS AllBDRTargeted2Transcriptome.sam | samtools sort -o AllBDRTargeted2Transcriptome.sorted.bam
samtools index AllBDRTargeted2Transcriptome.sorted.bam AllBDRTargeted2Transcriptome.sorted.bam.bai

##-------------------------------------------------------------------------

## 1. Calculates the allelic coverage of each variant; requires genome aligned sorted bam file and a phased VCF file (provided by isoPhase).
cd /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/lorals
calc_ase -b /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/6_minimap/AllBDRTargeted.sorted.bam \
-f /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/IsoSeq_IsoPhase.vcf
mv lorals_out/ase.tsv .
rmdir lorals_out

## 2. Annotates the output of calc_ase and assigns it a gene.
annotate_ase -i ase.tsv -b /lustre/projects/Research_Project-MRC148213/lsl693/references/annotation/gencode.v40.annotation.bed -o ./ase_annotated.tsv

## 3. Calculate transcript counts assigned to each haplotype
calc_asts -m quant -b /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/6_minimap/AllBDRTargeted.sorted.bam -x /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/6_minimap/AllBDRTargeted2Transcriptome.sorted.bam -i ase.tsv -o ./asts.tsv

## 4. Calculates the number of reads containing the REF or ALT allele assigned to each transcript.
process_asts -i asts.tsv_asts_quant_mod.tsv -g /lustre/projects/Research_Project-MRC148213/lsl693/references/human/gencode.v40.annotation.gene_transcript_ids.txt -o .


##-------------------------------------------------------------------------

## Rhbdf2 variant: chr17: 76471122; ref=C, alt=T

## 1. create mpileups tabulating the Rhbdf2 variant (take +1 and -1 nucleotide) 
samtools mpileup --output-QNAME -f ${GENCODE_FASTA} -r chr17:76471121-76471123 /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/6_minimap/AllBDRTargeted.sorted.bam > chr17_7647112_mpileup.txt

# 2. extract the transcripts (before collapse) with reference allele ("C")
grep 76471122 chr17_7647112_mpileup.txt | awk '{ print $7 }' | tr ',' '\n'  > chr17_7647112_refAlleleTranscripts.txt

# 3. extract the collapsed PBID and ccs reads referring to the reference allele
grep -w -f chr17_7647112_refAlleleTranscripts.txt ${cupcakeCollpased} | awk '{ print $1 }' > chr17_7647112_refAlleleTranscripts_PB.txt
grep -w -f chr17_7647112_refAlleleTranscripts_PB.txt ${cupcakeCollpasedreadStat} | awk '{ print $1 }' > chr17_7647112_refAlleleTranscripts_ccs.txt
sed -i 's/\/ccs//g' chr17_7647112_refAlleleTranscripts_ccs.txt

grep -w PB.4702 ${cupcakeCollpased} > chr17_7647112_all_PB.txt

# Read the refAllele file into a sorted array
mapfile -t refAllele < <(awk '{print $1}' /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/lorals/chr17_7647112_refAlleleTranscripts.txt | sort)

# Read the allAllele file, split column 2 by commas, and store it in an array
allAlleleTranscripts=($(awk -F'\t' '{print $2}' /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/lorals/chr17_7647112_all_PB.txt | tr ',' '\n' | sort))

# Find elements in allAlleleTranscripts that are not in refAllele
for allele in "${allAlleleTranscripts[@]}"; do
    if ! printf "%s\n" "${refAllele[@]}" | grep -q -w "$allele"; then
        echo "$allele" >> chr17_7647112_altAlleleTranscripts.txt
    fi
done
grep -w -f chr17_7647112_altAlleleTranscripts.txt ${cupcakeCollpased} | awk '{ print $1 }' > chr17_7647112_altAlleleTranscripts_PB.txt
grep -w -f chr17_7647112_altAlleleTranscripts_PB.txt ${cupcakeCollpasedreadStat} | awk '{ print $1 }' > chr17_7647112_altAlleleTranscripts_ccs.txt
sed -i 's/\/ccs//g' chr17_7647112_altAlleleTranscripts_ccs.txt


# 4. subset the tagged merged.sam with the reference transcript, convert to bam and index for igv
grep -f chr17_7647112_refAlleleTranscripts_ccs.txt merged.phased.partial.tagged.sam > chr17_7647112_refAlleleTranscripts_ccs.phased.partial.tagged.sam
samtools view -S -b chr17_7647112_refAlleleTranscripts_ccs.phased.partial.tagged.sam > chr17_7647112_refAlleleTranscripts_ccs.phased.partial.tagged.bam
samtools index chr17_7647112_refAlleleTranscripts_ccs.phased.partial.tagged.bam chr17_7647112_refAlleleTranscripts_ccs.phased.partial.tagged.bam.bai

grep -f chr17_7647112_altAlleleTranscripts_ccs.txt merged.phased.partial.tagged.sam > chr17_7647112_altAlleleTranscripts_ccs.phased.partial.tagged.sam
# "TTTCCC" last c
grep GACCATGGCCAAAGCCTTTCCCACAATGTCCCATCTGAGAGCCTTATGGATGGGCTC merged.phased.partial.tagged.sam > chr17_7647112_refAlleleTranscripts_ccs_C.phased.partial.tagged.sam
grep GACCATGGCCAAAGCCTTTCCGACAATGTCCCATCTGAGAGCCTTATGGATGGGCTC merged.phased.partial.tagged.sam > chr17_7647112_refAlleleTranscripts_ccs_G.phased.partial.tagged.sam
grep GACCATGGCCAAAGCCTTTCCTACAATGTCCCATCTGAGAGCCTTATGGATGGGCTC merged.phased.partial.tagged.sam > chr17_7647112_altAlleleTranscripts_ccs_T.phased.partial.tagged.sam
grep GACCATGGCCAAAGCCTTTCCAACAAAGTCCCATCTGAGAGCCTTATGGATGGGCTC merged.phased.partial.tagged.sam > chr17_7647112_altAlleleTranscripts_ccs_A.phased.partial.tagged.sam

awk '{print $1'} chr17_7647112_refAlleleTranscripts_ccs_G.phased.partial.tagged.sam > chr17_7647112_refAlleleTranscripts_ccs_G.txt
awk '{print $1'} chr17_7647112_refAlleleTranscripts_ccs_C.phased.partial.tagged.sam > chr17_7647112_refAlleleTranscripts_ccs_C.txt
awk '{print $1'} chr17_7647112_altAlleleTranscripts_ccs_T.phased.partial.tagged.sam > chr17_7647112_altAlleleTranscripts_ccs_T.txt
awk '{print $1'} chr17_7647112_altAlleleTranscripts_ccs_A.phased.partial.tagged.sam > chr17_7647112_altAlleleTranscripts_ccs_A.txt

# subset into ref and alt gtf
subset_fasta_gtf.py ${SQANTIGTF} --gtf -i chr17_7647112_refAlleleTranscripts_PB.txt -o chr17_7647112_refAlleleTranscripts_AllBDRTargeted.collapsed_corrected -d . 
mv AllBDRTargeted.collapsed_corrected_chr17_7647112_refAlleleTranscripts_AllBDRTargeted.collapsed_corrected.gtf chr17_7647112_refAlleleTranscripts_AllBDRTargeted.collapsed_corrected.gtf

subset_fasta_gtf.py ${SQANTIGTF} --gtf -i chr17_7647112_altAlleleTranscripts_PB.txt -o chr17_7647112_altAlleleTranscripts_AllBDRTargeted.collapsed_corrected -d . 
mv AllBDRTargeted.collapsed_corrected_chr17_7647112_altAlleleTranscripts_AllBDRTargeted.collapsed_corrected.gtf chr17_7647112_altAlleleTranscripts_AllBDRTargeted.collapsed_corrected.gtf

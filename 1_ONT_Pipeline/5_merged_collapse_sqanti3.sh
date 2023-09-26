#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=5_merged_collapse_sqanti3.o
#SBATCH --error=5_merged_collapse_sqanti3.e


# Batch2: realign with pbmm2 align and filter alignment 

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
source $SC_ROOT/1_ONT_Pipeline/bdr_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing

export dir=$WKD_ROOT/5_cupcake
#export samplename=ontBDRSortedNuclei
export samplename=ontBDR
mkdir -p ${dir}/5_align/combined ${dir}/6_collapse ${dir}/5_align/combined_fasta ${dir}/7_sqanti3


##-------------------------------------------------------------------------

#source activate nanopore
#replace_filenames_with_csv.py --copy --ext=filtered.bam -i=$WKD_ROOT/5_cupcake/5_align/Batch1 -f=$META_ROOT/ADBDR_Batch1_rename.csv -d=${dir}/5_align/combined &> ${dir}/5_align/combined/B1_copy.log
#replace_filenames_with_csv.py --copy --ext=filtered.bam -i=$WKD_ROOT/5_cupcake/5_align/Batch2 -f=$META_ROOT/ADBDR_Batch2_rename.csv -d=${dir}/5_align/combined &> ${dir}/5_align/combined/B2_copy.log
#replace_filenames_with_csv.py --copy --ext=filtered.fa -i=$WKD_ROOT/5_cupcake/5_align/Batch1 -f=$META_ROOT/ADBDR_Batch1_rename.csv -d=${dir}/5_align/combined_fasta &> ${dir}/5_align/combined_fasta/B1_copy.log
#replace_filenames_with_csv.py --copy --ext=filtered.fa -i=$WKD_ROOT/5_cupcake/5_align/Batch2 -f=$META_ROOT/ADBDR_Batch2_rename.csv -d=${dir}/5_align/combined_fasta &> ${dir}/5_align/combined_fasta/B2_copy.log


##-------------------------------------------------------------------------

# merge alignment
#echo "Collapsing..."
#allfilteredmapped=($(ls ${dir}/5_align/combined/*filtered.bam)) 
#ls ${allfilteredmapped[@]}
#source activate nanopore
#samtools merge -f ${dir}/6_collapse/ontBDR_mapped.filtered.sorted.bam ${allfilteredmapped[@]}

#cp ${SortedBAM} ${dir}/6_collapse
#mv ${dir}/6_collapse/merged.bam ${dir}/6_collapse/sortedNuclei.bam
#filter_alignment sortedNuclei ${dir}/6_collapse
#export ontBDRBAMfiltered=${dir}/6_collapse/ontBDR_mapped.filtered.sorted.bam
#export sortedBAMfiltered=${dir}/6_collapse/sortedNuclei_mapped.filtered.sorted.bam
#samtools merge -f ${dir}/6_collapse/${samplename}_mapped.filtered.sorted.bam ${ontBDRBAMfiltered} ${sortedBAMfiltered}

# collapse
#echo "Collapsing..."
#echo "Output: ${dir}/2_collapse/${samplename}_collapsed.gff"
#cd ${dir}/6_collapse
#source activate isoseq3
#isoseq3 collapse ${dir}/6_collapse/${samplename}_mapped.filtered.sorted.bam ${samplename}_collapsed.gff \
#  --min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons \
#  --log-level TRACE --log-file ${samplename}_collapsed.log
  
# demultiplex 
#source activate nanopore
#adapt_cupcake_to_ont.py ${dir}/5_align/combined_fasta -o ${samplename}

#demux_cupcake_collapse.py \
#  ${dir}/6_collapse/${samplename}_collapsed.read_stat.txt \
#  ${dir}/5_align/combined_fasta/${samplename}_sample_id.csv\
#  --dataset=ont
  
# sqanti3
echo "Running SQANTI3..."
source activate sqanti2_py3
cd ${dir}/7_sqanti3
python $SQANTI3_DIR/sqanti3_qc.py ${dir}/6_collapse/${samplename}_collapsed.gff \
  $GENOME_GTF $GENOME_FASTA -t 30 --CAGE_peak $CAGE_PEAK --polyA_motif_list $POLYA \
  --genename --isoAnnotLite --gff3 $GFF3 --skipORF --report skip &> ${samplename}_sqanti_qc.log

# sqanti3 filter 
filteringJson=$SQANTI3_DIR/utilities/filter/filter_default_reducecoverage.json
python $SQANTI3_DIR/sqanti3_filter.py rules ${samplename}_collapsed_classification.txt \
  --faa=${samplename}_collapsed_corrected.fasta \
  --gtf=${samplename}_collapsed_corrected.gtf \
  -j=${filteringJson} --skip_report &> ${samplename}_sqanti_filter.log
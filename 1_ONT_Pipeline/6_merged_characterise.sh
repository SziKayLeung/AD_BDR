#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --error=6_merged_characterise.e
#SBATCH --output=6_merged_characterise.o

# 22/03/2023: characterisation of merged ONT transcripts (BDR from Batch 1 and Batch 2)

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21

SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR
source $SC_ROOT/1_ONT_Pipeline/bdr_ont.config

FICLE_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/FICLE/
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
export PATH=$PATH:${LOGEN_ROOT}/compare_datasets
export PATH=$PATH:${LOGEN_ROOT}/target_gene_annotation
export PATH=$PATH:${LOGEN_ROOT}/merge_characterise_dataset
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous 
export PATH=$PATH:${FICLE_ROOT}
export PATH=$PATH:${FICLE_ROOT}/reference
source ${LOGEN_ROOT}/bash_pipelines/post_sqanti.sh

export dir=$WKD_ROOT/5_cupcake
export samplename=ontBDR
mkdir -p ${dir}/8_characterise

##-------------------------------------------------------------------------
## Characterisation with CPAT and colour by abundance

# LOGEN: subset cupcake classification file by target genes
# merge cupcake classification file with abundance
source activate nanopore
subset_quantify_filter_tgenes.R \
--classfile ${dir}/7_sqanti3/${samplename}_collapsed_RulesFilter_result_classification.txt \
--expression ${dir}/6_collapse/demux_fl_count.csv \
--target_genes ${TGENES_TXT} 

# filter cupcake classification file with minimum number of reads and counts
subset_quantify_filter_tgenes.R \
--classfile ${dir}/7_sqanti3/${samplename}_collapsed_RulesFilter_result_classification.txt \
--expression ${dir}/6_collapse/demux_fl_count.csv \
--target_genes ${TGENES_TXT} \
--filter --nsample=5 --nreads=10

# working variables
finalanno=${dir}/7_sqanti3/${samplename}_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt 
finaliso=${dir}/7_sqanti3/${samplename}_collapsed_RulesFilter_result_classification.targetgenes_filtered_isoforms.txt

# LOGEN: subset fasta and gtf using the finalised list of target gene isoforms
source activate sqanti2_py3
subset_fasta_gtf.py --gtf ${dir}/7_sqanti3/${samplename}_collapsed.filtered.gtf -i ${finaliso} -o counts_filtered
subset_fasta_gtf.py --fa ${dir}/7_sqanti3/${samplename}_collapsed.filtered.faa -i ${finaliso} -o counts_filtered

# run_cpat <input_fasta> <output_name> <output_dir>
run_cpat ${dir}/7_sqanti3/${samplename}_collapsed.filtered_counts_filtered.fa ${samplename} ${dir}/8_characterise

# extract_best_orf <sample> <root_dir>
extract_best_orf ${samplename} ${dir}/8_characterise

# classify reads by dataset 
source activate nanopore
classify_reads_bydataset.py -a=${dir}/6_collapse/demux_fl_count.csv -p=${GRPPHENOTYPE}

# colour_by_abundance <sample> <input_gtf> <abundance_file> <root_dir> <species>
colour_by_abundance ${samplename} \
${dir}/7_sqanti3/${samplename}"_collapsed.filtered_counts_filtered.gtf" \
${dir}/6_collapse/demux_fl_count_classified.csv ${dir}/8_characterise human

# colour_by_abundance <sample> <input_gtf> <abundance_file> <root_dir> <species>
# convert gtf to bed12
mkdir -p $4/bed12Files
convert_gtf_bed12 ${dir}/7_sqanti3/${samplename}"_collapsed.filtered_counts_filtered.gtf"

colour_transcripts_by_countandpotential.py \
  --bed ${dir}/7_sqanti3/${samplename}"_collapsed.filtered_counts_filtered_sorted.bed12" \
  --cpat ${dir}/8_characterise/CPAT/${samplename}".ORF_prob.best.tsv" \
  --noORF ${dir}/8_characterise/CPAT/${samplename}".no_ORF.txt" \
  --a ${dir}/6_collapse/demux_fl_count_classified.csv \
  --o ${samplename} \
  --dir ${dir}/8_characterise/bed12Files/ \
  --species human


##-------------------------------------------------------------------------
## Characterisation with FICLE

# FICLE: subset_gene_reference 
subset_gene_reference ${dir}/8_characterise 

# FICLE: full characterisation
# coloured by ONT and Iso-Seq
source activate sqanti2_py3
for g in ${TGENES[@]}; do

  echo $g
  mkdir -p ${dir}/8_characterise/TargetGenes
  mkdir -p ${dir}/8_characterise/TargetGenes/Log
  output_dir=${dir}/8_characterise/TargetGenes
  
  ficle.py --gene=$g \
  --ref=${dir}/8_characterise/TargetGenesRef/ \
  --i_bed=${dir}/8_characterise/bed12Files/${samplename}_concat_counts_coloured.bed12 \
  --i_gtf=${dir}/7_sqanti3/${samplename}_collapsed.filtered_counts_filtered.gtf \
  --i_class=${dir}/7_sqanti3/${samplename}_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt \
  --orf=${dir}/8_characterise/CPAT/${samplename}.ORF_prob.best.tsv   \
  --o_dir=$output_dir &> ${dir}/8_characterise/TargetGenes/Log/$g"_characterise.log"

done

# FICLE: extract reference information
extract_reference_info.py --r ${GENOME_GTF} --glist ${TGENES[@]} --short_read ${RNASEQ_COUNTS} --split ${dir}/8_characterise/TargetGenesRef 
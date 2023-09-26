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
#SBATCH --output=1c_merged_collapse_sqanti3.o2
#SBATCH --error=1c_merged_collapse_sqanti3.e2


# Batch2: realign with pbmm2 align and filter alignment 

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
source $SC_ROOT/1c_Merged_Pipeline/bdr_iso_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing
export PATH=$PATH:${LOGEN_ROOT}/assist_isoseq_processing
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous

export dir=${MERGED_ROOT}
export samplename=all_iso_ont_scn
mkdir -p ${dir}/2_collapse ${dir}/3_sqanti3

##-------------------------------------------------------------------------

# merge alignment
echo "Merging..."
all_iso_scn_mapped=($(ls ${dir}/1_align/*filtered.bam)) 
all_ont_mapped=($(ls ${ONT_ROOT}/5_cupcake/5_align/combined/*filtered.bam)) 
all_iso_ont_scn_mapped=(${all_iso_scn_mapped[@]} ${all_ont_mapped[@]})
ls ${all_iso_ont_scn_mapped[@]}

source activate nanopore
#samtools merge -f ${dir}/2_collapse/${samplename}.filtered.sorted.bam ${all_iso_ont_scn_mapped[@]}

# collapse
echo "Collapsing..."
echo "Output: ${dir}/2_collapse/${samplename}_collapsed.gff"
#cd ${dir}/2_collapse
#source activate isoseq3
#isoseq3 collapse ${dir}/2_collapse/${samplename}.filtered.sorted.bam ${samplename}_collapsed.gff \
#  --min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons \
#  --log-level TRACE --log-file ${samplename}_collapsed.log

# demultiplex 
source activate nanopore
# ont BDR targeted dataset 
adapt_cupcake_to_ont.py ${ONT_ROOT}/5_cupcake/5_align -o ontBDR
# ont sorted nuclei dataset
adapt_cupcake_to_ont.py ${dir}/1_align -o scn -i mapped.fa 
# combine header files from ONT datasets
ontBDRSampleID=${ONT_ROOT}/5_cupcake/5_align/combined_fasta/ontBDR_sample_id.csv
scnSampleID=${dir}/1_align/scn_sample_id.csv
cat ${ontBDRSampleID} <(tail -n+2 ${scnSampleID}) > ${dir}/1_align/ont_scn_sample_id.csv

# isoseq
demux_merged_cluster.py ${ISO_ROOT}/3_refine ${ISO_ROOT}/5_merged_cluster/AllBDRTargeted.clustered.cluster_report.csv

source activate sqanti2_py3
demux_ont_isoseq_cupcake_collapse.py \
  ${dir}/2_collapse/${samplename}_collapsed.read_stat.txt \
  ${dir}/1_align/ont_scn_sample_id.csv \
  ${ISO_ROOT}/5_merged_cluster/AllBDRTargeted.clustered.demuxed.cluster_report.csv \
  -d=${dir}/2_collapse

# sqanti3
echo "Running SQANTI3..."
source activate sqanti2_py3
cd ${dir}/3_sqanti3
python $SQANTI3_DIR/sqanti3_qc.py ${dir}/2_collapse/${samplename}_collapsed.gff \
$GENOME_GTF $GENOME_FASTA -t 30 --CAGE_peak $CAGE_PEAK --polyA_motif_list $POLYA \
--genename --isoAnnotLite --gff3 $GFF3 --skipORF --report skip &> ${samplename}_sqanti_qc.log

# sqanti3 filter 
filteringJson=$SQANTI3_DIR/utilities/filter/filter_default_reducecoverage.json
python $SQANTI3_DIR/sqanti3_filter.py rules ${samplename}_collapsed_classification.txt \
--faa=${samplename}_collapsed_corrected.fasta \
--gtf=${samplename}_collapsed_corrected.gtf \
-j=${filteringJson} --skip_report &> ${samplename}_sqanti_filter.log
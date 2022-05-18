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
#SBATCH --output=Targeted_ADBDR_Part2c.o
#SBATCH --error=Targeted_ADBDR_Part2c.e

# 11/05/2021: Run pipeline from refine to sqanti and tama (output: Targeted_ADBDR_Part2.o)
# 03/11/2021: Run kallisto with ADBDR RNA-Seq (30 samples), SQANTI3 with RNA-Seq support and TAMA filter (output: Targeted_ADBDR_Part2b.o)

#************************************* DEFINE GLOBAL VARIABLES
# File directories
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/IsoSeq
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/RNASeq
RNASeq_Filtered=/lustre/projects/Research_Project-193356/Project_10202/11_fastp_trimmed
DiffAnalysis_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Differential

# cd $PostIsoseq3_WKD/MAP; mkdir Individual
# cd $PostIsoseq3_WKD; mkdir SQANTI3 SQANTI3_TAMA_FILTER
# cd $DiffAnalysis_WKD; mkdir RNASeq_SQANTI3 TAPPAS_INPUT

# For Pooled Targeted
### Important order of BAM files is the same as the sample names
BATCH_NAMES=(Targeted_Seq_1 Targeted_Seq_2 Targeted_Seq_3)
RAWDIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/AD_BDR/Raw_Data

# Other input files and directory
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
PROBES=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/AD_BDR/Raw_Data/Probes

# For Demultiplexing Samples
Pooled_Samples_Targeted1=(A01 A02 A03 A04 A05 A06 A07 A08 A09 A10)
Pooled_Samples_Targeted2=(B01 B02 B03 B04 B05 B06 B07 B08 B09 B10)
Pooled_Samples_Targeted3=(C01 C02 C03 C04 C05 C06 C07 C08 C09 C11) # note C11 refers to C10 but not able to grep efficiently due to same name as barcode
Barcoded_Targeted1_config_file=$RAWDIR/Barcode_Configs/Isoseq_ADBDR_Targeted1_barcode.config
Barcoded_Targeted2_config_file=$RAWDIR/Barcode_Configs/Isoseq_ADBDR_Targeted2_barcode.config
Barcoded_Targeted3_config_file=$RAWDIR/Barcode_Configs/Isoseq_ADBDR_Targeted3_barcode.config
ALL_SAMPLES_NAMES=(A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 B01 B02 B03 B04 B05 B06 B07 B08 B09 B10 C01 C02 C03 C04 C05 C06 C07 C08 C09 C11)

# sourcing functions script
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/AD_BDR/1_Transcriptome_Annotation
source $FUNCTIONS/Targeted_ADBDR_Functions.sh

module load Miniconda2/4.3.21
################################################################################################
echo "#*************************************  Isoseq3 [Function 1, 2, 3, 6]"
echo "Already processed in batch (Mouse/Targeted_Mouse_Part1.sh) Functions 1 - 3, 6 for all batched runs"

## 4) run_targeted_REFINE ${BATCH} $Input_config_file $Input_Lima_sample $Input_LIMA_directory $Output_directory
#for i in ${Pooled_Samples_Targeted1[@]}; do echo "Processing $i*********"; run_targeted_REFINE $i $Barcoded_Targeted1_config_file Targeted_Seq_1 $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE; done
#for i in ${Pooled_Samples_Targeted2[@]}; do echo "Processing $i*********"; run_targeted_REFINE $i $Barcoded_Targeted2_config_file Targeted_Seq_2 $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE; done
#for i in ${Pooled_Samples_Targeted3[@]}; do echo "Processing $i*********"; run_targeted_REFINE $i $Barcoded_Targeted3_config_file Targeted_Seq_3 $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE; done

# 6) run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
#for i in ${ALL_SAMPLES_NAMES[@]}; do run_CLUSTER $i $Isoseq3_WKD/REFINE $Isoseq3_WKD/CLUSTER; done

################################################################################################
echo "#************************************* Isoseq3 and Post_Isoseq3 [Function 5, 7, 8]"
## 5) merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
#merging_at_refine $Isoseq3_WKD/REFINE $Isoseq3_WKD/MERGED_CLUSTER AllBDRTargeted A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 B01 B02 B03 B04 B05 B06 B07 B08 B09 B10 C01 C02 C03 C04 C05 C06 C07 C08 C09 C11

## 7) run_map_cupcakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
#run_map_cupcakecollapse AllBDRTargeted $Isoseq3_WKD/MERGED_CLUSTER $PostIsoseq3_WKD/MAP $PostIsoseq3_WKD/TOFU

## 8) demux_targeted <refine_dir> <input_cluster_report> <input_tofu_readstat> <output_path_file>
#demux_targeted $Isoseq3_WKD/REFINE $Isoseq3_WKD/MERGED_CLUSTER/AllBDRTargeted.clustered.cluster_report.csv $PostIsoseq3_WKD/TOFU/AllBDRTargeted.collapsed.read_stat.txt $PostIsoseq3_WKD/TOFU/AllBDRTargeted.Demultiplexed_Abundance.txt

################################################################################################
echo "#************************************* RNAseq [Function 9, 10]"
## 9) run_star <list_of_samples> <input_directory> <output_dir>
# Running separately on 30 samples in batch in Targeted_ADBDR_RNAseqStar.sh
# .e files were empty; cat all output files > ADBDR_Rnaseq_Star.o

## 10) ADBDR_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
#SAMPLES_NAMES=(BBN006_30024 BBN_24548 BBN002_28350 BBN_20194 BBN_24260 BBN003_30195 BBN002_29416 BBN_24287 BBN002_30035 BBN_15250 BBN002_29087 BBN003_26927 BBN_24289 BBN006_29902 BBN_25890 BBN002_29471 BBN_15247 BBN_18399 BBN002_26311 BBN_18405 BBN_24938 BBN006_29162 BBN_24253 BBN_19616 BBN_4240 BBN002_30920 BBN_15237 BBN_19632 BBN004_26227 BBN006_26347)
#ADBDR_merge_fastq $RNASeq_Filtered $PostIsoseq3_WKD/KALLISTO AllBDRTargeted

################################################################################################
echo "#************************************* RNASeq & IsoSeq [Function 11]"
## 11) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
run_kallisto AllBDRTargeted $PostIsoseq3_WKD/TOFU/AllBDRTargeted.collapsed.rep.fa $PostIsoseq3_WKD/KALLISTO $PostIsoseq3_WKD/KALLISTO

################################################################################################
echo "#************************************* SQANTI2 [Function 12]"
## 12) run_sqanti2 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/noexp/lncrna>
#run_sqanti2 AllBDRTargeted.collapsed AllBDRTargeted.collapsed.gff $PostIsoseq3_WKD/TOFU $RNASeq_WKD/MAPPED $PostIsoseq3_WKD/KALLISTO/AllBDRTargeted.mod.abundance.tsv $PostIsoseq3_WKD/TOFU/AllBDRTargeted.Demultiplexed_Abundance.txt $PostIsoseq3_WKD/SQANTI2 basic

## run_sqanti3 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/rnaseq/lncrna>
run_sqanti3 AllBDRTargeted.collapsed AllBDRTargeted.collapsed.gff $PostIsoseq3_WKD/TOFU $RNASeq_WKD/MAPPED $PostIsoseq3_WKD/KALLISTO/AllBDRTargeted.mod.abundance.tsv $PostIsoseq3_WKD/TOFU/AllBDRTargeted.Demultiplexed_Abundance.txt $PostIsoseq3_WKD/SQANTI3 genome

################################################################################################
echo "#************************************* TAMA filter [Function 13,14]"
## 13) TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
#TAMA_remove_fragments $PostIsoseq3_WKD/SQANTI2/AllBDRTargeted.collapsed_classification.filtered_lite.gtf AllBDRTargeted $PostIsoseq3_WKD/TAMA
TAMA_remove_fragments $PostIsoseq3_WKD/SQANTI3/AllBDRTargeted.collapsed_classification.filtered_lite.gtf AllBDRTargeted $PostIsoseq3_WKD/SQANTI3_TAMA_FILTER

## 14) TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
sqname=AllBDRTargeted.collapsed_classification.filtered_lite
#TAMA_sqanti_filter $PostIsoseq3_WKD/TAMA/AllBDRTargeted.bed $PostIsoseq3_WKD/SQANTI2 $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" AllBDRTargeted $PostIsoseq3_WKD/SQANTI_TAMA_FILTER
TAMA_sqanti_filter $PostIsoseq3_WKD/SQANTI3_TAMA_FILTER/AllBDRTargeted.bed $PostIsoseq3_WKD/SQANTI3 $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" AllBDRTargeted $PostIsoseq3_WKD/SQANTI3_TAMA_FILTER

################################################################################################
echo "#************************************* QC [Function 15]"
# parse_stats_per_sample <input_ccs.bam_dir> <Input_LIMA_directory> <output_prefix_name>
parse_stats_per_sample $Isoseq3_WKD/CCS $Isoseq3_WKD/LIMA AllBDRTargeted

#1. run_minimap2 <prefix_sample> <input_dir> <ERCC/mm10> <output_dir>
#2. run_minimap2 <input_fasta> <reference_fasta> <general> <output_dir> <output_prefix>
for i in ${ALL_SAMPLES_NAMES[@]}; do run_minimap2 $i $Isoseq3_WKD/CLUSTER mm10 $PostIsoseq3_WKD/MAP/Individual; done

# on_target_rate <input_probe_bed.file> <input_polished.hq.fasta> <input_mappped.fastq.sam> <output_file_name>
cd $Isoseq3_WKD/CLUSTER; gunzip *clustered.hq.fasta.gz*
for i in ${ALL_SAMPLES_NAMES[@]}; do on_target_rate $PROBES/FINAL_HUMAN.bed $Isoseq3_WKD/CLUSTER/$i".clustered.hq.fasta" $PostIsoseq3_WKD/MAP/Individual/$i".clustered.hq.fastq.sam" $PostIsoseq3_WKD/MAP/Individual/$i".fasta.sam.probe_hit.txt"

##################################################################################################
# isoannolite_generate <input_dir> <input_gtf> <input_class> <input_junc> <species> <output_dir> <output_name>
# generate IsoAnnotLite output required for TAPPAS after sqanti (final files )
isoannolite_generate $PostIsoseq3_WKD/SQANTI_TAMA_FILTER AllBDRTargeted_sqantitamafiltered.classification.final.gtf AllBDRTargeted_sqantitamafiltered.classification.txt AllBDRTargeted_sqantitamafiltered.junction.txt Human $PostIsoseq3_WKD/TAPPAS AllBDRTargeted_tappasannot_from_SQANTI2.gff3
# counts_subset_4tappas <input_class> <output_class>
counts_subset_4tappas $PostIsoseq3_WKD/SQANTI_TAMA_FILTER/AllBDRTargeted_sqantitamafiltered.classification.txt $PostIsoseq3_WKD/TAPPAS/AllBDRTargeted_sqantitamafiltered.expression.txt

#### preSQANTI filter
isoannolite_generate $PostIsoseq3_WKD/SQANTI2 AllBDRTargeted.collapsed_corrected.gtf AllBDRTargeted.collapsed_classification.txt AllBDRTargeted.collapsed_junctions.txt Human $PostIsoseq3_WKD/TAPPAS AllBDRTargeted_preSQF_tappasannot_from_SQANTI2.gff3
counts_subset_4tappas $PostIsoseq3_WKD/SQANTI2/AllBDRTargeted.collapsed_classification.txt $PostIsoseq3_WKD/TAPPAS/AllBDRTargeted.collapsed_expression.txt

# tappas on knight
cd /mnt/data1/Szi/TAPPAS_ADBDR
scp -r sl693@login.isca.ex.ac.uk:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq/TAPPAS/* .
scp -r sl693@login.isca.ex.ac.uk:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/AD_BDR/Raw_Data/ADBDR_PhenotypeTAPPAS.txt .
/mnt/data1/Aaron/sw/jre1.8.0_181/bin/java -jar /mnt/data1/Szi/tappAS.jar

# after running tappas on knight
cd $PostIsoseq3_WKD/TAPPAS; mkdir Results TAPPAS_output
scp -r sLeung@knight.ex.ac.uk:/mnt/data1/Szi/TAPPAS_ADBDR/*tsv* $PostIsoseq3_WKD/TAPPAS/Results

scp -r sLeung@knight.ex.ac.uk:/home/sLeung/tappasWorkspace/Projects/Project.094519432.tappas/InputData/input_normalized_matrix.tsv $PostIsoseq3_WKD/TAPPAS/Results

scp -r sLeung@knight.ex.ac.uk:/home/sLeung/tappasWorkspace/Projects/Project.094519432.tappas/Data/transcript_matrix.tsv $PostIsoseq3_WKD/TAPPAS/TAPPAS_output
scp -r sLeung@knight.ex.ac.uk:/home/sLeung/tappasWorkspace/Projects/Project.094519432.tappas/Data/transcript_matrix_raw.tsv $PostIsoseq3_WKD/TAPPAS/TAPPAS_output
scp -r sLeung@knight.ex.ac.uk:/home/sLeung/tappasWorkspace/Projects/Project.094519432.tappas/Data/gene_matrix_raw.tsv $PostIsoseq3_WKD/TAPPAS/TAPPAS_output

scp -r sLeung@knight.ex.ac.uk:/home/sLeung/tappasWorkspace/Projects/Project.094519432.tappas/Data/gene_transcripts.tsv $PostIsoseq3_WKD/TAPPAS/TAPPAS_output

scp -r sLeung@knight.ex.ac.uk:/home/sLeung/tappasWorkspace/Projects/Project.094519432.tappas/Data/result_gene_trans.tsv $PostIsoseq3_WKD/TAPPAS/TAPPAS_output

scp -r sLeung@knight.ex.ac.uk:/home/sLeung/tappasWorkspace/Projects/Project.094519432.tappas/Data/gene_matrix.tsv $PostIsoseq3_WKD/TAPPAS/TAPPAS_output
scp -r sLeung@knight.ex.ac.uk:/home/sLeung/tappasWorkspace/Projects/Project.094519432.tappas/Data/exp_factors.txt $PostIsoseq3_WKD/TAPPAS/TAPPAS_output

# note remove the comment from the header for each tsv

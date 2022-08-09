## ---------- Script -----------------
##
## Script name: 
##
## Purpose of script: Functions script for ADBDR dataset 
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
##
## 
##   
##
##


## ---------- Directory and input files -----------------

WKD_ROOT="/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/A_IsoSeq"

METADATA_DIR <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/0_metadata/A_IsoSeq"
TargetGene <- read.table(paste0(METADATA_DIR, "/TargetGenes.tsv"))[["V1"]]
targetedpheno <- read.csv(paste0(METADATA_DIR, "/Targeted_Sample_Demographics.csv")) 
probes <- read.table(paste0(METADATA_DIR,"/Probes/FINAL_HUMAN.bed")) %>% mutate(gene = word(V4,c(3),sep=fixed("_")))

CCS_input_dir <- paste0(WKD_ROOT,"/1_ccs")
LiMA_input_dir <- paste0(WKD_ROOT,"/2_lima")
REFINE_input_dir <- paste0(WKD_ROOT,"3_refine")
CLUSTER_input_dir <- paste0(WKD_ROOT,"4_cluster")
#CLUSTER_Merge <- read.csv(paste0(WKD_ROOT,"5_merged_cluster/AllBDRTargeted.clustered.cluster_report.csv"))


## ---------- SQANTI classification files -----------------

# Classification file
class.names.files <- paste0(WKD_ROOT,"/9_sqanti3/basic/AllBDRTargeted.collapsed_classification.filtered_lite_classification.txt")
gtf.names.files <- paste0(WKD_ROOT, "/9_sqanti3/basic/AllBDRTargeted.collapsed_classification.filtered_lite.gtf")
class.files <- SQANTI_class_preparation(class.names.files,"nstandard")

# target genes and remove 3'ISM
targeted.class.files <- targeted_remove_3ISM(TargetGene, class.files) 

# Subset class files by phenotype
targeted.class.files = group_descriptions(targeted.class.files) 
AD_targeted = subset_class_phenotype(targeted.class.files,targetedpheno,"AD") %>% dplyr::rename(AD_TotalFL = TotalFL)
Control_targeted = subset_class_phenotype(targeted.class.files,targetedpheno,"Control") %>% dplyr::rename(Control_TotalFL = TotalFL)


## ---------- RNA-Seq -----------------

### RNA-Seq gene level counts (STAR followed by RSEM - Darren Project 10202)
# Adapt file to only include target genes and in a format compatible for downstream python script 
Rnaseq_gene_counts = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/B_RNASeq/genes_TPM_matrix.txt"


## ---------- Target Rate -----------------

targetrate_input_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/A_IsoSeq/6b_target_rate"
Probes_input <- list.files(path = paste0(WKD_ROOT, "/6b_target_rate"), pattern = "fasta.sam.probe_hit.txt", full.names = T)
Probes_files <- lapply(Probes_input, function(x) read.table(x, header=T, as.is=T, sep="\t"))
names(Probes_files) <- list.files(path = paste0(WKD_ROOT, "/6b_target_rate"), pattern = "fasta.sam.probe_hit.txt")


##-------------------------------------------------------------------------

anno_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq/ISO_CHAR/TargetGenes"


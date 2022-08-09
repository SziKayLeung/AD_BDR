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


output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/01_figures_tables"

## ---------- TappAS input files -----------------

# input directory
WKD_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/Differential/TAPPAS_OUTPUT/"
TAPPAS_PHENOTYPE_DIR = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/0_metadata/A_IsoSeq/Tappas/"
TargetGene <- read.table(paste0(METADATA_DIR, "/TargetGenes.tsv"))[["V1"]]

# output files from running tappAS
TAPPAS_INPUT_DIR = list(
  iso = paste0(WKD_ROOT,"IsoSeq_Expression"),
  iso_nf = paste0(WKD_ROOT,"IsoSeq_Expression_nonfiltered"),
  rna = paste0(WKD_ROOT,"RNASeq_Expression"),
  rna_public = paste0(WKD_ROOT,"Large_RNASeq_Expression"),
  rna_public_nf = paste0(WKD_ROOT,"Large_RNASeq_Expression_nonfiltered")
)

# phenotype 
TAPPAS_PHENOTYPE = list(
  iso = paste0(TAPPAS_PHENOTYPE_DIR,"ADBDR_PhenotypeTAPPAS.txt"),
  rna = paste0(TAPPAS_PHENOTYPE_DIR,"ADBDR_MinusIntermediate_RNASeqPhenotypeTAPPAS.txt"),
  # core 96 samples including intermediate samples; note 4 missing samples removed 
  rna_all = paste0(TAPPAS_PHENOTYPE_DIR,"ADBDR_AllRNASeqPhenotypeTAPPAS.txt"),
  rna_public = paste0(TAPPAS_PHENOTYPE_DIR,"LargeRNASeq_PhenotypeTAPPAS_actual.txt")
)
phenotype <- lapply(TAPPAS_PHENOTYPE, function(x) read.table(x, header = T))


# differential analysis
difftrans <- list(isoseq = read_noiseq_files(paste0(WKD_ROOT,"/noiseq_difftrans.xlsx"),"noiseq_isoseq"),
                  rnaseq = read_noiseq_files(paste0(WKD_ROOT, "/noiseq_difftrans.xlsx"),"noiseq_rnaseq"))


## ---------- SQANTI classification files -----------------

# Classification file
ISOSEQ_WKD_ROOT="/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/A_IsoSeq/"
class.names.files <- paste0(ISOSEQ_WKD_ROOT,"9_sqanti3/basic/AllBDRTargeted.collapsed_classification.filtered_lite_classification.txt")
gtf.names.files <- paste0(ISOSEQ_WKD_ROOT, "9_sqanti3/basic/AllBDRTargeted.collapsed_classification.filtered_lite.gtf")
class.files <- SQANTI_class_preparation(class.names.files,"nstandard")

# target genes and remove 3'ISM
targeted.class.files <- targeted_remove_3ISM(TargetGene, class.files) 

## ---------- Iso-Seq Abundance -----------------

FL_reads <- read.csv(paste0(ISOSEQ_WKD_ROOT, "7_tofu/AllBDRTargeted.Demultiplexed_Abundance.txt")) 
colnames(FL_reads)[1] <- "isoform"
colnames(FL_reads)[-1] <- paste0("FL.",colnames(FL_reads)[-1])



# all TPM, no filtering of isoforms
TPM <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Differential/TAPPAS_OUTPUT/IsoSeq_Expression_nonfiltered/InputData/input_normalized_matrix.tsv") %>% rownames_to_column(var = "isoform")

metapheno <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/Samples/96_samples_all_v2_singletab.csv")
metapheno$BBN.ID <- gsub(".", "_", metapheno$BBN.ID, fixed=TRUE)
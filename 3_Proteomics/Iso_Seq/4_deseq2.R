## ---------- Script -----------------
##
## Purpose: perform differential analysis on AD-BDR Iso-Seq targeted datasets using linear regression
## Transcript level separate analysis
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
##
## ---------- Notes ------------------
## two ONT batches: Batch 1 - Nov 2022, Batch 2 - March 2023 


## ---------- packages -----------------

suppressMessages(library("dplyr"))
suppressMessages(library("DESeq2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
suppressMessages(library("ggrepel"))
suppressMessages(library("wesanderson"))
suppressMessages(library("cowplot"))
suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))


## ---------- source functions -----------------

LOGEN_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/run_DESeq2.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/plot_transcript_level.R"))
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/pthemes.R"))

label_group <- function(genotype){
  if(genotype %in% c("Case","CASE")){group = "AD"}else{
    if(genotype %in% c("Control","CONTROL")){group = "Control"}}
  return(group)
}

## ---------- input -----------------

# directory names
dirnames <- list(
  bdr = "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/",
  output = "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/01_figures_tables/"
)

# target gene
TargetGene <- toupper(c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7",
                        "Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf"))

# read input files
input_files <- list(
  isoPhenotype = paste0(dirnames$bdr, "/0_metadata/A_IsoSeq/Targeted_Sample_Demographics.csv"),
  ontPhenotype = paste0(dirnames$bdr, "0_metadata/B_ONT/Selected_ONTTargeted_BDR.csv"), 
  # merged data (unfiltered but contain only target genes)
  classfiles = paste0(dirnames$bdr, "A_IsoSeq/9_sqanti3/basic/AllBDRTargeted.collapsed_classification.filtered_lite_classification.txt")
)


# proteomics input 
load(file = paste0(SC_ROOT,"/3_Proteomics/proteinInput.RData"))

input <- list()
input$ontPhenotype <- read.table(input_files$ontPhenotype, sep = ",", header = T) %>% mutate(BBN.ID = str_replace(BBN.ID,"_","")) 
input$isoPhenotype <- read.table(input_files$isoPhenotype, sep = ",", header = T)
input$classfiles <- SQANTI_class_preparation(input_files$classfiles,"ns") %>% filter(associated_gene %in% TargetGene)
input$p.classfiles <- proteinInput$pFiltered.class.files
row.names(input$classfiles) <- input$classfiles$isoform


# phenotype
phenoFiles <- list(
  isoGrp = input$isoPhenotype %>% mutate(sample = Column.ID, group = Phenotype),
  isoBraak = input$isoPhenotype %>% mutate(sample = Column.ID, group = as.numeric(Braak))
)

### protein 
# aggregate sum by same peptide sequence
pFL <- proteinInput$t.class.files %>% filter(associated_gene %in% TargetGene) %>% dplyr::select(base_acc, contains("FL.")) 
pFLsum <- aggregate(. ~ base_acc, pFL, sum)

# datawrangle for input to run_DESeq2()
expressionFiles <- list(
  iso = input$classfiles %>% dplyr::select(starts_with("FL.")),
  isoProtein = pFLsum %>% tibble::column_to_rownames(., var = "base_acc")
)
row.names(expressionFiles$iso) <- input$classfiles$isoform

## ---------- ONT: Creating DESeq2 object and analysis -----------------

ResTran <- list(
  tWald = run_DESeq2(test="Wald",expressionFiles$iso,phenoFiles$isoGrp,threshold=10,controlname="Control",design="case_control",groupvar="factor"),
  tWaldBraak = run_DESeq2(test="Wald",expressionFiles$iso,phenoFiles$isoBraak,threshold=10,design="case_control",groupvar="numeric"),
  pWald = run_DESeq2(test="Wald",expressionFiles$isoProtein,phenoFiles$isoGrp,threshold=10,controlname="Control",design="case_control",groupvar="factor"),
  pWaldBraak = run_DESeq2(test="Wald",expressionFiles$isoProtein,phenoFiles$isoBraak,threshold=10,design="case_control",groupvar="numeric")
)

annoResTran <- list(
  tWald = anno_DESeq2(ResTran$tWald,input$classfiles,phenoFiles$isoGrp,controlname="Control",level="transcript",sig=0.1),
  tWaldBraak = anno_DESeq2(ResTran$tWaldBraak,input$classfiles,phenoFiles$isoBraak,controlname="0",level="transcript",sig=0.1),
  pWald = anno_DESeq2(ResTran$pWald,input$classfiles,phenoFiles$isoGrp,controlname="Control",level="transcript",sig=0.1),
  pWaldBraak = anno_DESeq2(ResTran$pWaldBraak,input$classfiles,phenoFiles$isoBraak,controlname="0",level="transcript",sig=0.1)
)

 

## ---------- Output -----------------

saveRDS(annoResTran, file = paste0(dirnames$output, "IsoProtein_DESeq2TranscriptLevel.RDS"))
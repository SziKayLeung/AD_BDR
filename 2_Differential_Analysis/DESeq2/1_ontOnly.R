## ---------- Script -----------------
##
## Purpose: perform differential analysis on AD-BDR ONT targeted datasets using linear regression
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

# read input files
input_files <- list(
  ontPhenotype = paste0(dirnames$bdr, "0_metadata/B_ONT/Selected_ONTTargeted_BDR.csv"), 
  # merged data (unfiltered but contain only target genes)
  classfiles = paste0(dirnames$bdr, "D_ONT/5_cupcake/7_sqanti3/ontBDR_collapsed_RulesFilter_result_classification.targetgenes_counts.txt")
)


input <- list()
input$ontPhenotype <- read.table(input_files$ontPhenotype, sep = ",", header = T) %>% mutate(BBN.ID = str_replace(BBN.ID,"_","")) 
input$classfiles <- SQANTI_class_preparation(input_files$classfiles,"ns")

# datawrangle for input to run_DESeq2()
expressionFiles <- list(
  ontB1 = input$classfiles %>% dplyr::select(starts_with("B1.")),
  ontB2 = input$classfiles %>% dplyr::select(starts_with("B2.")),
  ontMerged = input$classfiles %>% dplyr::select(starts_with("B1."),starts_with("B2."))
)

# sum the expression across ONT Batch 1 and Batch 2
expressionFiles$ontMergedSum <- rbind(expressionFiles$ontB1 %>% tibble::rownames_to_column(var = "isoform") %>% 
         reshape2::melt(id="isoform",value.name="counts",variable.name="sample") %>% 
         mutate(sample_id = str_remove(sample,"B1.")),
      expressionFiles$ontB2 %>% tibble::rownames_to_column(var = "isoform") %>% 
        reshape2::melt(id="isoform",value.name="counts",variable.name="sample") %>% 
        mutate(sample_id = str_remove(sample,"B2."))) %>% 
  group_by(isoform,sample_id) %>% summarise(counts = sum(counts), .groups = 'drop')


expressionFiles$ontMergedSum  <- tidyr::spread(expressionFiles$ontMergedSum,"sample_id","counts") %>% 
  tibble::column_to_rownames(., var = "isoform")

# phenotype
phenoFiles <- list(
  ontB1Grp = input$ontPhenotype %>% mutate(sample = paste0("B1.",sample), group = phenotype),
  ontB2Grp =input$ontPhenotype %>% mutate(sample = paste0("B2.",sample), group = phenotype),
  ontB1Braak = input$ontPhenotype %>% mutate(sample = paste0("B1.",sample), group = as.numeric(BraakTangle_numeric)),
  ontB2Braak = input$ontPhenotype %>% mutate(sample = paste0("B2.",sample), group = as.numeric(BraakTangle_numeric)),
  ontMerged = rbind(input$ontPhenotype %>% mutate(sample = paste0("B1.",sample), group = phenotype),
                         input$ontPhenotype %>% mutate(sample = paste0("B2.",sample), group = phenotype)),
  ontMergedBraak = rbind(input$ontPhenotype %>% mutate(sample = paste0("B1.",sample), group = as.numeric(BraakTangle_numeric)),
                   input$ontPhenotype %>% mutate(sample = paste0("B2.",sample), group = as.numeric(BraakTangle_numeric))),
  ontGrpBraak = input$ontPhenotype %>% mutate(sample = sample, group = as.numeric(BraakTangle_numeric))
)

# remove outliers
phenoFiles <- lapply(phenoFiles, function(x) x %>% filter(!BBN.ID %in% c("BBN00226311","BBN00326927","BBN9944")) %>% dplyr::select(sample, group))


## ---------- ONT: Creating DESeq2 object and analysis -----------------

ResTran <- list(
  # ONT Batch 1
  B1Wald = run_DESeq2(test="Wald",expressionFiles$ontB1,phenoFiles$ontB1Grp,threshold=10,controlname="Control",design="case_control",groupvar="factor"),
  B1WaldBraak = run_DESeq2(test="Wald",expressionFiles$ontB1,phenoFiles$ontB1Braak,threshold=10,design="case_control",groupvar="numeric"),
  # ONT Batch 2
  B2Wald = run_DESeq2(test="Wald",expressionFiles$ontB2,phenoFiles$ontB2Grp,threshold=10,controlname="Control",design="case_control",groupvar="factor"),
  B2WaldBraak = run_DESeq2(test="Wald",expressionFiles$ontB2,phenoFiles$ontB2Braak,threshold=10,design="case_control",groupvar="numeric"),
  # Merged
  ontMerged = run_DESeq2(test="Wald",expressionFiles$ontB2,phenoFiles$ontMerged,threshold=10,controlname="Control",design="case_control",groupvar="numeric"),
  ontMergedBraak =  run_DESeq2(test="Wald",expressionFiles$ontMerged,phenoFiles$ontMergedBraak,threshold=10,design="case_control",groupvar="numeric"),
  ontMergedBraakSum = run_DESeq2(test="Wald",expressionFiles$ontMergedSum,phenoFiles$ontGrpBraak,threshold=10,design="case_control",groupvar="numeric")
  
)

annoResTran <- list(
  B1Wald = anno_DESeq2(ResTran$B1Wald,input$classfiles,phenoFiles$ontB1Grp,controlname="Control",level="transcript",sig=0.05),
  B1WaldBraak = anno_DESeq2(ResTran$B1WaldBraak,input$classfiles,phenoFiles$ontB1Braak,controlname="0",level="transcript",sig=0.05),
  B2Wald = anno_DESeq2(ResTran$B2Wald,input$classfiles,phenoFiles$ontB2Grp,controlname="Control",level="transcript",sig=0.05),
  B2WaldBraak = anno_DESeq2(ResTran$B2WaldBraak,input$classfiles,phenoFiles$ontB2Braak,controlname="0",level="transcript",sig=0.05),
  ontMerged =  anno_DESeq2(ResTran$ontMerged,input$classfiles,phenoFiles$ontMerged,controlname="Control",level="transcript",sig=0.05),
  ontMergedBraak =  anno_DESeq2(ResTran$ontMergedBraak,input$classfiles,phenoFiles$ontMergedBraak,controlname="0",level="transcript",sig=0.05),
  ontMergedBraakSum =  anno_DESeq2(ResTran$ontMergedBraakSum,input$classfiles,phenoFiles$ontGrpBraak,controlname="0",level="transcript",sig=0.05)
)

## ---------- Output -----------------

saveRDS(annoResTran, file = paste0(dirnames$output, "Ont_DESeq2TranscriptLevel.RDS"))
write.table(annoResTran$ontMerged$anno_res %>% filter(padj < 0.05) %>% dplyr::select(isoform),
            paste0(dirnames$output, "OntMerged_DESeq2_sigiso.txt"),quote = F, row.names = F, col.names = F)
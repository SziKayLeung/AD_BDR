## ---------- Script -----------------
##
## Purpose: perform differential analysis on AD-BDR ONT targeted datasets using linear regression
## Transcript level separate analysis
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
# https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html


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

LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "/transcriptome_stats/read_sq_classification.R"))
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
  bdr = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/",
  output = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/01_figures_tables/"
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
  ontMerged <- input$classfiles %>% dplyr::select(starts_with("B1."),starts_with("B2."))
)

phenoFiles <- list(
  ontB1Grp = input$ontPhenotype %>% mutate(sample = paste0("B1.",sample), group = phenotype),
  ontB2Grp =input$ontPhenotype %>% mutate(sample = paste0("B2.",sample), group = phenotype),
  ontB1Braak = input$ontPhenotype %>% mutate(sample = paste0("B1.",sample), group = as.numeric(BraakTangle_numeric)),
  ontB2Braak = input$ontPhenotype %>% mutate(sample = paste0("B2.",sample), group = as.numeric(BraakTangle_numeric)),
  ontMergedBraak = rbind(input$ontPhenotype %>% mutate(sample = paste0("B1.",sample), group = as.numeric(BraakTangle_numeric)),
                   input$ontPhenotype %>% mutate(sample = paste0("B2.",sample), group = as.numeric(BraakTangle_numeric)))
)

# remove outliers
phenoFiles <- lapply(phenoFiles, function(x) x %>% filter(!BBN.ID %in% c("BBN00226311","BBN00326927","BBN9944")) %>% dplyr::select(sample, group))


## ---------- ONT: Creating DESeq2 object and analysis -----------------

ResTran <- list(
  # ONT Batch 1
  B1Wald = run_DESeq2(test="Wald",expressionFiles$ontB1,phenoFiles$ontB1Grp,threshold=10,controlname="Control",design="case_control",groupvar="factor"),
  B1WaldBraak = run_DESeq2(test="Wald",expressionFiles$ontB1,phenoFiles$ontB1Braak,threshold=10,design="case_control",groupvar="numeric"),
  # ONT Batch 2
  B2Wald = run_DESeq2(test="Wald",expressionFiles$ontB2,phenoFiles$ontB2Grp,threshold=10,controlname="Control",design="case_control"),groupvar="factor",
  B2WaldBraak = run_DESeq2(test="Wald",expressionFiles$ontB2,phenoFiles$ontB2Braak,threshold=10,design="case_control",groupvar="numeric"),
  # Merged
  ontMerged =  run_DESeq2(test="Wald",expressionFiles$ontMerged,phenoFiles$ontMergedBraak,threshold=10,design="case_control",groupvar="numeric")
)

annoResTran <- list(
  B1Wald = anno_DESeq2(ResTran$B1Wald,input$classfiles,phenoFiles$ontB1Grp,controlname="Control",level="transcript",sig=1.0),
  B1WaldBraak = anno_DESeq2(ResTran$B1WaldBraak,input$classfiles,phenoFiles$ontB1Braak,controlname="0",level="transcript",sig=1.0),
  B2Wald = anno_DESeq2(ResTran$B2Wald,input$classfiles,phenoFiles$ontB2Grp,controlname="Control",level="transcript",sig=1.0),
  B2WaldBraak = anno_DESeq2(ResTran$B2WaldBraak,input$classfiles,phenoFiles$ontB2Braak,controlname="0",level="transcript",sig=1.0),
  ontMerged =  anno_DESeq2(ResTran$ontMerged,input$classfiles,phenoFiles$ontMergedBraak,controlname="0",level="transcript",sig=1.0)
)


## ---------- Output -----------------

saveRDS(annoResTran, file = paste0(dirnames$output, "Ont_DESeq2TranscriptLevel.RDS"))
write.table(annoResTran$ontMerged$anno_res %>% filter(padj < 0.05) %>% dplyr::select(isoform),
            paste0(dirnames$output, "OntMerged_DESeq2_sigiso.txt"),quote = F, row.names = F, col.names = F)


## ---------- ONT: Creating DESeq2 object and analysis -----------------

setorderAD = c("Control","Intermediate_1","Intermediate_2","AD")
plot_transexp_overtime("SNCA",annoResTran$B2Wald$norm_counts,design="multiple_case_control",show="specific",isoSpecific=c("PB.70995.7653"), setorder=setorderAD)
plot_transexp_overtime("APP",annoResTran$B2WaldBraak$norm_counts,design="multiple_case_control",show="specific",isoSpecific=c("PB.57985.21"))
plot_transexp_overtime("BIN1",annoResTran$ontMerged$norm_counts,design="multiple_case_control",show="specific",isoSpecific=c("PB.50706.8"))
TargetedDESeq$ontResTranAnno$ontMerged$anno_res <- TargetedDESeq$ontResTranAnno$ontMerged$anno_res %>% filter(padj < 0.05)  

pdf("DifferentialAnalysis.pdf")
for(i in 1:nrow(TargetedDESeq$ontResTranAnno$ontMerged$anno_res)){
  gene = TargetedDESeq$ontResTranAnno$ontMerged$anno_res$associated_gene[i]
  iso = TargetedDESeq$ontResTranAnno$ontMerged$anno_res$isoform[i]
  print(plot_transexp_overtime(gene,annoResTran$ontMerged$norm_counts,design="multiple_case_control",show="specific",isoSpecific=c(iso)))
}
dev.off()

p <- annoResTran$ontMerged$norm_counts %>% filter(isoform == iso) %>% 
  filter(group %in% c(0,6)) %>%
  ggplot(., aes(x = reorder(sample, normalised_counts), y = normalised_counts, fill = group)) +
  geom_bar(stat = "identity") +
  facet_grid(~group, scales="free") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)


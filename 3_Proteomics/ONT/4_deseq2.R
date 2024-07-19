#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## proteogenonomics analysis following paper review
##    process merged (ONT + Iso-Seq) targeted data after collapsing by protein ORF following G.Shenkyman pipeline
##    re-determined reference isoform using most abundant transcript after collapsing by ORF
##    differential transcript analysis at protein level 
## --------------------------------

## ---------- packages -----------------

suppressMessages(library("utils"))
LOGEN_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/run_DESeq2.R"))

## ---------- config file -----------------

# directory names
dirnames <- list(
  bdr = "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/",
  output = "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/01_figures_tables/",
  protein = "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/C_Proteomics/ONT/"
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
classFiles <- input$classfiles

protein = list(
  cpat = read.table(paste0(dirnames$protein,"5_calledOrfs/ADBDR_best_orf.tsv"), sep ="\t", header = T),
  t2p.collapse = read.table(paste0(dirnames$protein,"6_refined_database/ADBDR_orf_refined.tsv"), sep = "\t", header = T)
)
TargetGene <- toupper(c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7",
                        "Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf"))

## ---------- process data from pipeline -----------------

## data-wrangle orf_refined.tsv 
# filter to target genes
# generate column of the number of transcripts collapsed by delimiting the pb_accs column
# output: a table of the pb_accs, base_acc (the representative collapsed isoform selected by G.Shenkyman pipeline) and numtxCollapsed
protein$t2p.collapse <- protein$t2p.collapse %>% 
  filter(gene%in%TargetGene) %>% 
  mutate(numtxCollapsed = count.fields(textConnection(as.character(pb_accs)), sep = "|"))
char <- strsplit(as.character(protein$t2p.collapse $pb_accs), '|', fixed = T)
t2p.collapse.dissected <- data.frame(pb_accs=unlist(char), base_acc=rep(protein$t2p.collapse$base_acc, sapply(char, FUN=length)))
protein$t2p.collapse <- merge(t2p.collapse.dissected,protein$t2p.collapse[,c("base_acc","numtxCollapsed")])

## re-determine representative colalsped isoform: using ONT abundance (sum across all samples) rather than arbitrary (G.Shenkyman pipeline)
# take the ONT_sum read counts from the classification file
# max = grouping by the base_acc (i.e. the previously selected isoform), select the rows with the maximum ONT FL reads
# create an index to remap and create a "corrected_acc" column with the corresponding isoform that has the highest number of ONT FL reads
protein$t2p.collapse <- merge(protein$t2p.collapse,classFiles[,c("isoform","nreads")],by.x = "pb_accs", by.y = "isoform", all.x = TRUE)
max = protein$t2p.collapse %>% group_by(base_acc) %>% filter(nreads == max(nreads))
idx <- match(protein$t2p.collapse$base_acc, max$base_acc)
protein$t2p.collapse = transform(protein$t2p.collapse, corrected_acc = ifelse(!is.na(idx), as.character(max$pb_accs[idx]), base_acc))

## include in the original classification file the collapsed PB.ID 
classFiles <- merge(classFiles, protein$t2p.collapse[,c("pb_accs","numtxCollapsed","base_acc","corrected_acc")], by.x = "isoform", by.y = "pb_accs", all.x = TRUE)

## Statistics
message("Total number of RNA transcripts: ", nrow(classFiles))
message("Number of protein-coding RNA transcripts: ", nrow(classFiles[!is.na(classFiles$base_acc),]))
message("Number of non-protein-coding RNA transcripts: ", nrow(classFiles[is.na(classFiles$base_acc),]))


## ---------- differential expression analysis: proteogenomics -----------------

## 1. expression file
# aggregate sum by same peptide sequence
# datawrangle for input to run_DESeq2()
pFL <- classFiles%>% filter(!is.na(corrected_acc) & associated_gene %in% TargetGene) 
ontB2pFL <- pFL %>% select(corrected_acc, contains("B2"), -nreads)
expressionFiles <- list(
  ontB2Protein =  aggregate(. ~ corrected_acc, ontB2pFL, sum) %>% tibble::column_to_rownames(., var = "corrected_acc")
)
colnames(expressionFiles$ontB2Protein) <- word(colnames(expressionFiles$ontB2Protein),c(2),sep=fixed("."))


## 2. phenotype file 
phenoFiles <- list(
  ontB2Braak = input$ontPhenotype  %>% mutate(group = as.numeric(BraakTangle_numeric), sample = word(sample,c(1), sep = fixed("."))),
  ontB2Grp = input$ontPhenotype %>% filter(phenotype != "Intermediate") %>% mutate(group = phenotype, sample = word(sample,c(1), sep = fixed(".")))
)

# remove outliers
phenoFiles <- lapply(phenoFiles, function(x) x %>% filter(!BBN.ID %in% c("BBN00226311","BBN00326927","BBN9944")) %>% dplyr::select(sample, group))

## ---------- Creating DESeq2 object and analysis -----------------

ResTran <- list(
  # ONT Batch 2
  B2WaldGroup = run_DESeq2(test="Wald",expressionFiles$ontB2,phenoFiles$ontB2Grp,threshold=10,controlname="Control",design="case_control",groupvar="factor"),
  B2WaldBraak = run_DESeq2(test="Wald",expressionFiles$ontB2Protein,phenoFiles$ontB2Braak,threshold=10,design="case_control",groupvar="numeric")
)

annoResTran <- list(
  B2WaldGroup = anno_DESeq2(ResTran$B2WaldGroup,input$classfiles,phenoFiles$ontB2Grp,controlname="Control",level="transcript",sig=0.1),
  B2WaldBraak = anno_DESeq2(ResTran$B2WaldBraak,input$classfiles,phenoFiles$ontB2Braak,controlname="0",level="transcript",sig=0.1)
)

## ---------- Output -----------------
write.table(classFiles, paste0(dirnames$bdr, "D_ONT/5_cupcake/7_sqanti3/ontBDR_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered_pCollapsed.txt"),sep="\t",quote = F)
write.table(protein$t2p.collapse, paste0(dirnames$protein,"6_refined_database/ADBDR_orf_refined.tsv"),sep="\t",quote = F)
saveRDS(annoResTran, file = paste0(dirnames$output, "/DESeq2ONTProteinLevel.RDS"))

#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: process output from long-read proteogenomics pipeline
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------
## 

## ------------ packages ------------

SC_ROOT = "/lustre/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR"
LOGEN <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN, "transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN, "transcriptome_stats/plot_basic_stats.R"))
source(paste0(LOGEN, "merge_characterise_dataset/run_ggtranscript.R"))


## ------------ directory paths ------------

root_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/"
dirnames <- list(
  root = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/C_Proteomics/",
  targ_iso_metadata = paste0(root_dir,"AD_BDR/0_metadata/A_IsoSeq/"),
  targ_ont_metadata = paste0(root_dir,"AD_BDR/0_metadata/B_ONT/"),
  targ_ont_root = paste0(root_dir, "AD_BDR/D_ONT/5_cupcake/")
)


## ------------ input ------------

inputAll <- list(
  # transcript SQANTI3 classification file
  t.class.files = SQANTI_class_preparation(paste0(dirnames$root,"2_longread_processed/AllBDRTargeted.collapsed_classification.filtered_lite_classification.txt"),"ns"),
  # protein unfiltered long-read proteogenomics classification file
  pUnfiltered.class.files = pSQANTI_class_preparation(paste0(dirnames$root,"7_classified_protein/ADBDR_unfiltered.protein_classification.tsv")),
  # protein filtered long-read proteogeonomics classification file
  pFiltered.class.files = pSQANTI_class_preparation(paste0(dirnames$root,"7_classified_protein/ADBDR.sqanti_protein_classification.tsv")),
  # list of transcripts collapsed by protein reading frame
  t2p.collapse = fread(paste0(dirnames$root,"6_refined_database/ADBDR_orf_refined.tsv")),
  # gtf of aligned peptide ORF
  peptide_orf = paste0(dirnames$root, "5_calledOrfs/ADBDR.gtf"),
  # cpat
  cpat_output = fread(paste0(dirnames$root,"5_calledOrfs/ADBDR.ORF_prob.best.tsv")),
  no_cpat_orf = fread(paste0(dirnames$root, "5_calledOrfs/ADBDR.no_ORF.txt"),header = F)
)

TargetGene <- toupper(c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7",
                        "Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf"))

# subset to only list of target genes
input <- list(
  t.class.files = inputAll$t.class.files[which(inputAll$t.class.files$associated_gene%in%TargetGene),],
  pUnfiltered.class.files = inputAll$pUnfiltered.class.files[which(inputAll$pUnfiltered.class.files$tx_gene%in%TargetGene),],
  pFiltered.class.files = inputAll$pFiltered.class.files[which(inputAll$pFiltered.class.files$tx_genename%in%TargetGene),],
  t2p.collapse = inputAll$t2p.collapse %>% filter(inputAll$t2p.collapse$gene%in%TargetGene),
  peptide_orf = as.data.frame(rtracklayer::import(inputAll$peptide_orf)),
  cpat_output = inputAll$cpat_output,
  no_coat_orf = inputAll$no_cpat_orf
) 

input$gtf <- as.data.frame(rtracklayer::import(paste0(dirnames$root,"2_longread_processed/AllBDRTargeted.collapsed_classification.filtered_lite.gtf")))
input$peptide_orf <- input$peptide_orf %>% mutate(isoform = word(transcript_id,c(1),sep=fixed("_")))
input$merged_peptide_gtf <- rbind(input$gtf[,c("seqnames","strand","start","end","type","transcript_id","gene_id")] ,
                         input$peptide_orf[,c("seqnames","strand","start","end","type","transcript_id","gene_id")])


# number of RNA transcripts collapsed by protein sequence
input$t2p.collapse <- input$t2p.collapse %>% mutate(numtxCollapsed = count.fields(textConnection(pb_accs), sep = "|"))
char <- strsplit(as.character(input$t2p.collapse$pb_accs), '|', fixed = T)
t2p.collapse.dissected <- data.frame(pb_accs=unlist(char), base_acc=rep(input$t2p.collapse$base_acc, sapply(char, FUN=length)))
input$t2p.collapse <- merge(t2p.collapse.dissected,input$t2p.collapse[,c("base_acc","numtxCollapsed")])
input$t.class.files <- merge(input$t.class.files, input$t2p.collapse, by.x = "isoform", by.y = "pb_accs")

# phenotype
phenotype <- list(
  targ_iso = read.csv(paste0(dirnames$targ_iso_metadata, "Targeted_Sample_Demographics.csv"), header = T),
  targ_ont = read.csv(paste0(dirnames$targ_ont_metadata, "Selected_ONTTargeted_BDR.csv"), header = T),
  scn =  read.csv(paste0(root_dir, "AD_BDR/0_metadata/RBSCNPhenotype.csv"), header = T)
)
phenotype$batch_ont <- rbind(phenotype$targ_ont %>% mutate(sample = paste0("B1.",sample)),phenotype$targ_ont %>% mutate(sample = paste0("B2.",sample)))


## ------------ output ------------ 

proteinInput <- input
save(proteinInput, file = paste0(SC_ROOT,"/3_Proteomics/proteinInput.RData"))

# aggregate sum by same peptide sequence
pFL <- input$t.class.files %>% dplyr::select(base_acc, contains("FL."))
pFLsum <- aggregate(. ~ base_acc, pFL, sum)

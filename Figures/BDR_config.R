# Szi Kay Leung: sl693@exeter.ac.uk

LOGEN <- "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen"
source(paste0(LOGEN,"/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN,"/target_gene_annotation/summarise_gene_stats.R"))
source(paste0(LOGEN,"/compare_datasets/dataset_identifer.R"))
source(paste0(LOGEN,"/merge_characterise_dataset/run_ggtranscript.R"))

## --------------------------- 
# variables
TargetGene <- toupper(c("Abca1","Sorl1","Mapt","Bin1","Tardbp","App","Abca7",
                "Ptk2b","Ank1","Fyn","Clu","Cd33","Fus","Picalm","Snca","Apoe","Trpa1","Rhbdf2","Trem2","Vgf"))

## --------------------------- 
# directory names
root_dir <- "/lustre/projects/Research_Project-MRC148213/lsl693/"
dirnames <- list(
  # targeted sequencing (Iso-Seq, ONT)
  targ_iso_metadata = paste0(root_dir,"AD_BDR/0_metadata/A_IsoSeq/"),
  targ_ont_metadata = paste0(root_dir,"AD_BDR/0_metadata/B_ONT/"),
  targ_ont_root = paste0(root_dir, "AD_BDR/D_ONT/5_cupcake/"),
  targ_iso_root = paste0(root_dir, "AD_BDR/A_IsoSeq/"),
  targ_anno = paste0(root_dir,"AD_BDR/D_ONT/5_cupcake/8_characterise/TargetGenes/"),
  targ_output = paste0(root_dir, "/AD_BDR/01_figures_tables/"),
  
  # targeted sequencing + single cell
  targ_merged = paste0(root_dir, "/AD_BDR/E_MergedTargeted/"),
  references = paste0(root_dir,"references/annotation")
)


## ------------- Phenotype files -------------------

phenotype <- list(
  targ_iso = read.csv(paste0(dirnames$targ_iso_metadata, "Targeted_Sample_Demographics.csv"), header = T),
  targ_ont = read.csv(paste0(dirnames$targ_ont_metadata, "Selected_ONTTargeted_BDR.csv"), header = T),
  scn =  read.csv(paste0(root_dir, "AD_BDR/0_metadata/RBSCNPhenotype.csv"), header = T)
)
phenotype$batch_ont <- rbind(phenotype$targ_ont %>% mutate(sample = paste0("B1.",sample)),phenotype$targ_ont %>% mutate(sample = paste0("B2.",sample)))


## --------------------------- 
# Final classification file
class.names.files <- list(
  targ_ont_all = paste0(dirnames$targ_ont_root, "7_sqanti3/ontBDR_collapsed_RulesFilter_result_classification.targetgenes_counts.txt"),
  targ_ont_filtered = paste0(dirnames$targ_ont_root, "7_sqanti3/ontBDR_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt"),
  targ_merged = paste0(dirnames$targ_merged, "3_sqanti3/all_iso_ont_scn_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt")
) 
class.files <- lapply(class.names.files, function(x) SQANTI_class_preparation(x,"nstandard"))

## adapt to human
### for downstream subsetting of the global transcriptome by WT and TG mice sample
#sub_class.files <- lapply(wholesamples, function(x) subset_class_by_sample(class.files$glob_iso,x))
#names(sub_class.files) <- wholesamples

#group_class.files <- bind_rows(
#  subset_class_phenotype(class.files$glob_iso, phenotype$whole_rTg4510_iso, "WT") %>% mutate(Dataset = "WT"),
#  subset_class_phenotype(class.files$glob_iso, phenotype$whole_rTg4510_iso, "TG") %>% mutate(Dataset = "TG")
#)


## ---------- DESeq2 results -----------------

TargetedDESeq <- list(
  ontResTranAnno = readRDS(file = paste0(dirnames$targ_output, "/Ont_DESeq2TranscriptLevel.RDS"))
) 


## ---------- merged targeted results overview -----------------

Targeted <- list(
  Genes = array(read.table(paste0(dirnames$targ_ont_metadata, "TargetGenes.tsv"))[["V1"]]),
  
  # final transcript classifications of all merged samples
  Gene_class = lapply(list.files(path = dirnames$targ_anno, pattern = "Final_Transcript_Classifications.csv", recursive = TRUE, full = T), 
                      function(x) read.csv(x)),
  
  # CPAT 
  cpat = read.table(paste0(dirnames$targ_ont_root, "8_characterise/CPAT/ontBDR.ORF_prob.best.tsv"), header = T) %>% 
    mutate(coding_status = ifelse(Coding_prob >= 0.364, "Coding","Non_Coding")),
  
  # noORF file from CPAT
  noORF = read.table(paste0(dirnames$targ_ont_root, "8_characterise/CPAT/ontBDR.no_ORF.txt")) %>% 
    mutate(coding_status = "No_ORF") %>% `colnames<-`(c("seq_ID", "coding_status")),
  
  # reference 
  ref_gencode = read.csv(paste0(dirnames$targ_ont_root, "8_characterise/TargetGenesRef/TargetGene_Reference_LengthNum.csv")),
  ref_altcon = read.csv(paste0(dirnames$targ_ont_root, "8_characterise/TargetGenesRef/TargetGene_Reference_AltConExons.csv")),
  
  # ont abundance 
  ont_abundance = class.files$targ_ont_filtered %>% dplyr::select(isoform, contains("B2."), contains("B1.")) %>% mutate(annot_transcript_id = isoform)
  
)
names(Targeted$Gene_class) = lapply(list.files(path = dirnames$targ_anno, pattern = "Final_Transcript_Classifications.csv", recursive = TRUE), 
                                             function(x) word(x, c(1), sep = fixed("/")))
Merged_gene_class_df <- all_summarise_gene_stats(Targeted$Gene_class, class.files$targ_ont_filtered, Targeted$cpat, Targeted$noORF, Targeted$Genes)


## --------------------------
# gtf
gtf <- list(
  targ_merged = rtracklayer::import(paste0(dirnames$targ_merged, "3_sqanti3/all_iso_ont_scn_collapsed.filtered_counts_filtered.gtf")),
  ref_target = rtracklayer::import(paste0(dirnames$references,"/gencode.v40.annotation.20Targets.gtf"))
)
gtf <- lapply(gtf, function(x) as.data.frame(x))

gtf$targ_merged <- rbind(gtf$targ_merged[,c("seqnames","strand","start","end","type","transcript_id","gene_id")] ,
                         gtf$ref_target[,c("seqnames","strand","start","end","type","transcript_id","gene_id")])

gtf$targ_merged <- gtf$targ_merged %>% mutate(co = paste0(seqnames,":",start,"-",end))

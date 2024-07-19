library("dplyr")
library("stringr")
library("ggplot2")

massID = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/Proteomics/1_raw/30_samples_correspondance.tsv", header = T)
isoID = read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/AD_BDR/Raw_Data/Targeted_Sample_Demographics.csv")
complete_massID = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/Proteomics/1_raw/filenames.tsv")

massID$BBN.ID <- gsub(".", "_", massID$BBN.ID, fixed = TRUE)
isoID$BBN.ID <- gsub(".", "_", isoID$BBN.ID, fixed = TRUE)
complete_massID <- complete_massID %>% mutate("No." = word(word(V1,c(1),sep=fixed(".")),c(3),sep=fixed("_")))
AllID <- merge(merge(massID, isoID),complete_massID)

gencode = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/Proteomics/9_metamorpheus/ADBDR_gencode_search_results/Task1SearchTask/AllQuantifiedProteinGroups.gencode.tsv",fill = TRUE, header = T, sep = "\t",row.names = NULL)
colnames(gencode) <- colnames(gencode)[2:ncol(gencode)]    

filtered = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/Proteomics/9_metamorpheus/ADBDR_filtered_search_results/Task1SearchTask/AllQuantifiedProteinGroups.filtered.tsv",fill = TRUE, header = T, sep = "\t",row.names = NULL)
colnames(filtered) <- colnames(filtered)[2:ncol(filtered)]   

filtered_peptides = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/Proteomics/9_metamorpheus/ADBDR_filtered_search_results/Task1SearchTask/AllQuantifiedPeptides.tsv", header = T, sep = "\t",row.names = NULL)
colnames(filtered_peptides) <- colnames(filtered_peptides)[2:ncol(filtered_peptides)]   

nrow(gencode[grepl("SNCA", gencode$Gene),])
nrow(filtered[grepl("SNCA", filtered$Gene),])
intensity <- filtered[grepl("SNCA", filtered$Gene),] %>% select(Protein.Accession,Gene,starts_with("Intensity")) %>% reshape2::melt() %>% 
  mutate(MassID = as.numeric(word(variable, c(4), sep = fixed("_")))) %>% 
  full_join(., massID, by = c("MassID" = "No.")) %>% 
  mutate(BBN.ID = gsub(".", "_", BBN.ID, fixed = TRUE)) %>% 
  full_join(., isoID[,c("Phenotype","BBN.ID","Sample.ID")],by="BBN.ID") 

ggplot(intensity, aes(x=Protein.Accession,y = value, colour = Phenotype)) + 
  geom_boxplot() + 
  geom_point(aes(fill = Phenotype), size = 2, shape = 21, position = position_jitterdodge()) 


# SQANTI, TAMA filtered file
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")
PostIsoSeq_root_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq"
class.names.files <- paste0(PostIsoSeq_root_dir, "/SQANTI3_nojunc/AllBDRTargeted.collapsed_classification.filtered_lite_classification.txt")
class.files <- SQANTI_class_preparation(class.names.files,"nstandard")

class.files %>% filter(isoform %in% c("PB.7941.178","PB.7941.28")) %>% select(isoform,starts_with("FL.")) %>% reshape2::melt() %>% 
  full_join(.,isoID[,c("Column.ID","Phenotype")],by=c("variable"="Column.ID") ) %>% 
  ggplot(.,aes(x =isoform,y=value,colour=Phenotype)) + geom_boxplot() +
  geom_point(aes(fill = Phenotype), size = 2, shape = 21, position = position_jitterdodge()) 



####################################
# Peptide quantification 
intensity_peptides <- filtered_peptides %>% filter(Gene.Names == "SNCA") %>% select(Base.Sequence, starts_with("Intensity")) %>% 
  reshape::melt() %>%
  mutate(MassID = as.numeric(word(variable, c(4), sep = fixed("_")))) %>% 
  full_join(., massID, by = c("MassID" = "No.")) %>% 
  mutate(BBN.ID = gsub(".", "_", BBN.ID, fixed = TRUE)) %>% 
  full_join(., isoID[,c("Phenotype","BBN.ID","Column.ID")],by="BBN.ID") 

ggplot(intensity_peptides, aes(x=Base.Sequence,y = value, colour = Phenotype)) + 
  geom_boxplot() + 
  geom_point(aes(fill = Phenotype), size = 2, shape = 21, position = position_jitterdodge()) +
  coord_flip()

intensity_peptides %>% filter(Base.Sequence == "KLIRNSASR")
class.files %>% filter(isoform %in% c("PB.7941.178","PB.7941.28")) %>% select(isoform,starts_with("FL.")) %>% reshape2::melt() %>% 
  full_join(.,isoID[,c("Column.ID","Phenotype")],by=c("variable"="Column.ID")) %>% group_by(variable) %>% 
  tally(value) %>% full_join(., intensity_peptides %>% filter(Base.Sequence == "KLIRNSASR"),by=c("variable"="Column.ID")) %>%
  ggplot(.,aes(x=n,y=value)) + geom_point()

View(class.files %>% filter(isoform %in% c("PB.7941.178","PB.7941.28")) %>% select(isoform,starts_with("FL.")) %>% reshape2::melt())

###################
quant_peptides <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/Proteomics/9_metamorpheus/ADBDR_filtered_search_results/Task1SearchTask/QuantifiedPeptides.tsv", sep = "\t", header = T, as.is=T,fill=T)

quant_peptides[grepl("KLIRNSASR",quant_peptides$Sequence),] %>% select(starts_with("Intensity")) %>% 
  reshape2::melt() %>% mutate(MassID = as.numeric(word(variable, c(4), sep = fixed("_")))) %>% 
  full_join(., massID, by = c("MassID" = "No.")) %>% 
  mutate(BBN.ID = gsub(".", "_", BBN.ID, fixed = TRUE)) %>% 
  full_join(., isoID[,c("Phenotype","BBN.ID","Column.ID")],by="BBN.ID") 

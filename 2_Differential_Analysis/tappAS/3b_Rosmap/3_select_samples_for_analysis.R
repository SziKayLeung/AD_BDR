#!/usr/bin/env Rscript
# Szi Kay Leung
# 11/07/2022: Generate a file of the RNA-Seq samples from extreme braak scores (0,1,5,6) that are already downloaded for subsequent trimming and kallisto alignment

# library
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))

# load files 
metadir = "/lustre/projects/Research_Project-MRC148213/lsl693/Scripts/AD_BDR/2_Differential_Analysis/ROSMAP/Metadata/"
final_merged <- read.csv(final_merged,paste0(metadir, "/synapseid_phenotype.csv"), quote=F)

# Samples that are already downloaded
downloaded1 = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/ROSMAP/RNASeq/Synapse/downloaded_round1.txt") %>% 
  mutate(specimenID = word(V1, c(1), sep = fixed(".")))
downloaded2 = read.table("/lustre/projects/Research_Project-MRC190311/ROSMAP/RNASeq/Synapse/downloaded_round2.txt") %>% 
  mutate(specimenID = word(V1, c(1), sep = fixed(".")))

#******************** 1. Subset samples for trimming and kalisto alignment 
# Classify if samples are downloaded 
# Subset extreme AD and control samples (braak score 0, 1, 5, and 6)
final_merged_extreme <- final_merged %>% 
  mutate(Status = ifelse(specimenID %in% c(downloaded1$specimenID, downloaded2$specimenID),"Downloaded","Not Downloaded")) %>% 
  filter(braaksc %in% c("0","1","5","6"))

# Subset the samples that are already downloaded
final_merged_extreme_downloaded <- final_merged_extreme[final_merged_extreme$Status == "Downloaded",]
length(final_merged_extreme_downloaded$specimenID)

# Output
write.table(final_merged_extreme_downloaded$specimenID,paste0(metadir, "/synapseid_extreme_downloaded.txt"), quote=F,row.names = F, col.names = F)

## ---------------------------
##
## Script name: 1b_mass_spec_transfer_validate.R
##
## Purpose of script: 
##    1. Cross validating the BDR samples downloaded and sequenced using mass-spectrometry
##    2. Creating Experimental_Design.tsv file in the raw directory for downstream metamorpheus
##
## Author: Szi Kay Leung
##
## Date Created: 02-08-2022
##
## Email: sl693@exeter.ac.uk
##
## ---------------------------
##
## Notes: 
##   BDR_downloaded file = a catalogue of files in ISCA directory of the raw files downloaded using filesender
##   BDR meta file from Valentin/Ehsan Pishva
##   Total 96 samples 
##
## ---------------------------

## ---------------------------

# library
suppressMessages(library("stringr"))
suppressMessages(library("dplyr"))


## ---------------------------

# load BDR meta file
BDR_meta <- read.csv("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/Proteomics/0_metadata/BDR_96_samples_proteomics_meta.csv")

# load file containing names of files downloaded
BDR_downloaded <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/Proteomics/1_raw/ADBDR_all_files.txt")
BDR_downloaded <- BDR_downloaded %>% mutate(ID = as.integer(word(word(V1,c(3),sep=fixed("_")),c(1),sep=fixed("."))))


## ---------------------------

# check that all the files downloaded match the BDR metafile
# i.e. no files missing
setdiff(BDR_meta$No.,BDR_downloaded$ID)
setdiff(BDR_downloaded$ID,BDR_meta$No.)
length(BDR_downloaded$ID)

## --------------------------- 

# create a file for Experimental.csv (Metamorpheus)
Experimental_File <- merge(BDR_downloaded, BDR_meta, by.x = "ID", by.y = "No.") %>% 
  select(V1,Case.Control) %>% 
  mutate(FileName = V1,
         Condition = Case.Control) %>% 
  select(FileName,Condition)

Experimental_File$Sample <- c(seq(1:length(which(Experimental_File$Condition == "control"))),
  seq(1:length(which(Experimental_File$Condition == "intermediate"))),
  seq(1:length(which(Experimental_File$Condition == "case"))))


Experimental_File$Biorep <- "1"  
Experimental_File$Fraction <- "1"  
Experimental_File$Techrep <- "1"  

write.table(Experimental_File,"/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/1_raw/1_proteomics_raw/ExperimentalDesign.tsv",
            quote=F,col.names=T,row.names=F,sep="\t")

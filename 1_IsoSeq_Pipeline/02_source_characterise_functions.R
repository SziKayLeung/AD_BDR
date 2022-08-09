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

## ---------- Packages -----------------

suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(rjson)) # json files
suppressMessages(library(plyr)) # revalue
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(reshape))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(dplyr))
suppressMessages(library(stringr)) 
suppressMessages(library(viridis)) 
suppressMessages(library(wesanderson)) 
suppressMessages(library(extrafont))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))
suppressMessages(library(tibble))
suppressMessages(library(VennDiagram))
suppressMessages(library(directlabels))
suppressMessages(library(cowplot))
suppressMessages(library(readxl))
suppressMessages(library(ggdendro))
suppressMessages(library(pheatmap))
suppressMessages(library(extrafont))
suppressMessages(loadfonts())

## ----------Functions-----------------

# load all the functions
source_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential"
file.sources = list.files(path = source_dir, pattern="*.R", full = T)
sapply(file.sources,source,.GlobalEnv)

## ----------Plot colours-----------------

# plot label colour
label_colour <- function(genotype){
  if(genotype == "Braak0"){colour = wes_palette("Royal1")[1]}else{
    if(genotype == "Braak1"){colour = wes_palette("Royal2")[3]}else{
      if(genotype == "Braak6"){colour = wes_palette("Royal1")[2]}else{
        if(genotype == "novel"){colour = wes_palette("Darjeeling1")[4]}else{
          if(genotype == "known"){colour = wes_palette("Darjeeling1")[5]}else{
            if(genotype == "targeted"){colour = wes_palette("Darjeeling1")[2]}else{
              if(genotype == "whole"){colour = wes_palette("Darjeeling1")[1]}else{
                if(genotype == "whole+targeted"){colour = wes_palette("Darjeeling2")[1]}else{
                  if(genotype == "Control"){colour = wes_palette("Royal1")[2]}else{
                    if(genotype == "Case"){colour = wes_palette("Royal1")[1]}else{
                      if(genotype == c("Yes")){colour = alpha(wes_palette("Cavalcanti1")[4],0.8)}else{
                        if(genotype == c("No")){colour = alpha(wes_palette("Cavalcanti1")[5],0.5)}else{
                }}}}}}}}}}}}
  return(colour)
}


## ---------- Dataset defined functions -----------------

# number_of_reads
# Aim: huge wrapper function to read in multiple files from CCS, LIMA and REFINE directory, and output 5 plots
number_of_reads <- function(){
  
  Reads <- input_isoseq_files("AllBDRTargeted_CCS_output.csv","AllBDRTargeted_LIMA_summary.csv")
  
  # classify Braak in Reads
  Reads$Braak <- lapply(Reads$Sample.ID, function(x)
    if(x %in% targetedpheno[targetedpheno$Braak == "6","Sample.ID"]){"6"
    } else if (x %in% targetedpheno[targetedpheno$Braak == "0","Sample.ID"]){"0"
    } else if (x %in% targetedpheno[targetedpheno$Braak == "1","Sample.ID"]){"1"
    } else {"Batch"
    }
  )
  Reads$Braak <- unlist(Reads$Braak)
  
  # classify samples into Batches
  Reads$Batch <- lapply(Reads$Sample.ID, function(x)
    if(grepl("A", x)){"Batch 1"
    } else if (grepl("B", x)){"Batch 2"
    } else if (grepl("C", x)){"Batch 3"
    } else if (grepl("Targeted_Seq_1", x)){"Batch 1"
    } else if (grepl("Targeted_Seq_2", x)){"Batch 2"
    } else if (grepl("Targeted_Seq_3", x)){"Batch 3"
    } else {"Batch"
    }
  )
  Reads$Batch  <- unlist(Reads$Batch)
  
  # Genotype
  Reads = Reads %>% 
    mutate(Genotype = ifelse(Braak == "6","AD","Control")) %>%
    mutate(Genotype = ifelse(Braak == "Batch","NA",Genotype))
  
  # classify Genotype in CCS_values_mod for downstream plotting
  CCS_values_mod$Braak <- lapply(CCS_values_mod$sample, function(x)
    if(x %in% targetedpheno[targetedpheno$Braak == "6","Sample.ID"]){"6"
    } else if (x %in% targetedpheno[targetedpheno$Braak == "0","Sample.ID"]){"0"
    } else if (x %in% targetedpheno[targetedpheno$Braak == "1","Sample.ID"]){"1"
    } else {"Batch"
    }
  )
  CCS_values_mod <<- CCS_values_mod
  
  return(Reads)
}


QC_yield_plot <- function(){
  #Reads <- number_of_reads()
  #write.csv(Reads,paste0(OUTPUT_DIR,"/ADBDR_IsoSeqTargetedReadsStats.csv"), row.names = F)
  Reads <- read.csv(paste0(OUTPUT_DIR,"/ADBDR_IsoSeqTargetedReadsStats.csv"))
  Reads$Description <- factor(Reads$Description, 
                              levels = c("Polymerase Reads","CCS Reads","FL Reads","FLNC Reads","Poly-A FLNC Reads","Transcripts"))
  Reads$Genotype <- factor(Reads$Genotype, c("NA","Control","AD"))
  
  p1 <- Reads %>% filter(Description != "Transcripts") %>% 
    ggplot(., aes(x = Description, y = value, colour = Batch, group = Batch)) +
    geom_line(linetype = "dotted") + geom_point() +  mytheme + theme(legend.position = c(0.8,0.8)) + labs(x = "", y = "Number of Reads (Thousands)") +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_colour_discrete(name = "",labels = c("Batch 1 (n = 10)","Batch 2 (n = 10) ","Batch 3 (n = 10)"))
  
  # Number of Transcripts - Braak Staging
  p2 <- Reads %>% filter(Description == "Transcripts") %>% 
    ggplot(., aes(x = Braak, y = value, colour = Braak)) +
    geom_boxplot() + geom_point(size = 3) +  mytheme + 
    labs(x = "Braak staging", y = "Number of FL Transcripts (Thousands)") +
    scale_color_manual(values = c(label_colour("Braak0"),label_colour("Braak1"),label_colour("Braak6"))) + theme(legend.position = "none") +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3,accuracy = 1))
  
  # Number of Transcripts - Genotype
  p3 <- Reads %>% filter(Description == "Transcripts") %>% 
    ggplot(., aes(x = Genotype, y = value, colour = Genotype)) +
    geom_boxplot() + geom_point(size = 3) +  mytheme + 
    labs(x = "Phenotype", y = "Number of FL Transcripts (Thousands)") +
    scale_color_manual(values = c(label_colour("Case"),label_colour("Control"))) + theme(legend.position = "none") +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3,accuracy = 1)) 
  
  p4 <- p3 + facet_grid(~Batch) + theme(strip.background = element_blank())
  
  cat("Number of CCS Reads in Batch 1, 2 and 3:", 
      sum(Reads[Reads$Description == "CCS Reads","value"])/1000000,"M \n")
  cat("Mean Number of CCS Reads in Batch 1, 2 and 3:", 
      mean(Reads[Reads$Description == "CCS Reads","value"])/1000000, "M \n")
  cat("SD of CCS Reads in Batch 1, 2 and 3:", 
      sd(Reads[Reads$Description == "CCS Reads","value"])/1000, "K \n")
  cat("Min of CCS Reads in Batch 1, 2 and 3:", 
      min(Reads[Reads$Description == "CCS Reads","value"])/1000, "K \n")
  cat("Max of CCS Reads in Batch 1, 2 and 3:", 
      max(Reads[Reads$Description == "CCS Reads","value"])/1000, "K \n")

  cat("Number of Poly-A FLNC Reads in Batch 1, 2 and 3:", 
      sum(Reads[Reads$Description == "Poly-A FLNC Reads","value"])/1000000,"M \n")
  cat("Mean Number of Poly-A FLNC Reads in Batch 1, 2 and 3:", 
      mean(Reads[Reads$Description == "Poly-A FLNC Reads","value"])/1000, "K \n")
  cat("SD of Poly-A FLNC Reads in Batch 1, 2 and 3:", 
      sd(Reads[Reads$Description == "Poly-A FLNC Reads","value"])/1000, "K \n")
  cat("Min of Poly-A FLNC Reads in Batch 1, 2 and 3:", 
      min(Reads[Reads$Description == "Poly-A FLNC Reads","value"])/1000, "K \n")
  cat("Max of Poly-A FLNC Reads in Batch 1, 2 and 3:", 
      max(Reads[Reads$Description == "Poly-A FLNC Reads","value"])/1000, "K \n")
  
  # not normally distributed therefore wilcoxon rank sum test
  #with(Reads_plot %>% filter(Description == "Transcripts"), shapiro.test(value[Phenotype == "WT"]))
  #with(Reads_plot %>% filter(Description == "Transcripts"), shapiro.test(value[Phenotype == "TG"]))
  #var.test(value ~ Phenotype,Reads_plot %>% filter(Description == "Transcripts")) #cannot assume variance
  #wilcox.test(value ~ Genotype,Reads %>% filter(Description == "Transcripts")) 
  
  # correlation of FL transcripts and RIN
  #transcript_RIN <- merge(Reads_plot %>% filter(Description == "Transcripts"),tg4510_samples, by.x = "sample",by.y = "Sample.ID", all.x = T)
  #shapiro.test(transcript_RIN$value) # spearman's rank
  #shapiro.test(transcript_RIN$RIN) 
  #cor.test(transcript_RIN$value,transcript_RIN$RIN, method = "spearman", exact = FALSE)
  
  return(list(p1,p2,p3))
}

on_target_plot <- function(){
  Probes <- merge(ldply(Probes_files , function(x) nrow(x)),
                  ldply(Probes_files , function(x) length(which(x$num_base_overlap != "0"))),by = ".id") %>%
    `colnames<-`(c("file", "Total_mapped_reads", "reads_probe_hit")) %>% 
    mutate(perc = reads_probe_hit/Total_mapped_reads * 100) %>%
    mutate(sample = word(.$file, c(1), sep = fixed("."))) %>% 
    full_join(., targetedpheno, by = c("sample" = "Sample.ID"))

  Probes$Phenotype <- factor(Probes$Phenotype, c("NA","Control","AD"))
  
  p1<- ggplot(Probes, aes(x = Phenotype, y = perc, colour = Phenotype)) + geom_boxplot() +
    geom_point(size = 3, position = position_jitterdodge()) + 
    facet_grid(~Batch) + 
    scale_colour_manual(values = c(label_colour("Case"),label_colour("Control"))) +
    mytheme + labs(x = "Phenotype", y = "On-Target Rate (%)") + theme(legend.position = "none")
  
  for(i in 1:3){cat("Mean on target rate in Batch",i,":", 
                    mean(Probes[Probes$Batch == paste0("Batch ",i),"perc"]),"\n")}
  return(p1)
}

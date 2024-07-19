#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: basic functions for dataset comparisons 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------------------------------


## ---------- packages -----------------

suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("uniProt.ws"))
suppressMessages(library("httr"))


## ---------- load data -----------------

# import Ehsan's results from differential proteomic analysis

load("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/C_Proteomics/13_ehsan/BDR_Braak_DEPr_results_limma.Rdata")

# import txt file tabulating uniProtID accession ID and entry name
uniProtID <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/references/human/homo_sapiens_uniprot.txt", sep = "\t")
uniProtID$Description <- lapply(uniProtID$V3, function(x) strsplit(as.character(x),"OX")[[1]][1])
uniProtID$geneName <- stringr::str_remove(uniProtID$V4,"GN=")
colnames(uniProtID) <- c("header","uniProtID","description","V4","geneName")


## ---------- data-wrange data -----------------

# arrange by FDR
# annotate accession ID by entry name
top.table.bdr <- merge(top.table.bdr,uniProtID[,c("uniProtID","geneName")], by.x = 0, by.y = "uniProtID")
top.table.bdr <- top.table.bdr %>% arrange(adj.P.Val)
colnames(top.table.bdr)[1] <- "uniProtID_accession"
rownames(top.table.bdr) <- NULL

## ---------- plot -----------------

nrow(top.table.bdr[top.table.bdr$adj.P.Val < 0.05,])
# plot top-rank
# note Braak column in pheno file is wrong
y[row.names(top.table.bdr)[1], ] %>% reshape2::melt(value.name = "counts") %>% merge(pheno[,c("Gender","NFT_Braak")], by = 0) %>% 
  mutate(Braak = stringr::str_remove(NFT_Braak, "Braak tangle stage")) %>%
  ggplot(., aes(x = Braak, y = counts)) + geom_boxplot() +
  theme_classic() +
  labs(x = "Braak stage", y = "Counts")


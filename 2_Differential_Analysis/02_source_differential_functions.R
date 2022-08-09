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
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/plot_tappas_analysis.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/General/5_TappAS_Differential/sqanti_general.R")

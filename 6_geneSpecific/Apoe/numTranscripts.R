#!/usr/bin/env Rscript
## ----------Script-----------------
## 
## SKLeung: APOE focus; number of transcripts
## --------------------------------


## ------- packages ---------------- 

library("dplyr")
library("stringr")
library("ggplot2")

LOGEN_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "transcriptome_stats/read_sq_classification.R"))


## ---------- input -----------------

dirnames <- list(
  bdr = "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/",
  ficle = "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/D_ONT/5_cupcake/8_characterise/TargetGenes/TargetGenes"
)

# read input files
input_files <- list(
  ontPhenotype = paste0(dirnames$bdr, "0_metadata/B_ONT/Selected_ONTTargeted_BDR.csv"), 
  # merged data (unfiltered but contain only target genes)
  classfiles = paste0(dirnames$bdr, "D_ONT/5_cupcake/7_sqanti3/ontBDR_collapsed_RulesFilter_result_classification.txt"),
  # counts 
  counts = paste0(dirnames$bdr, "D_ONT/5_cupcake/6_collapse/demux_fl_count.csv")
)

input <- list()
input$ontPhenotype <- read.table(input_files$ontPhenotype, sep = ",", header = T) %>% mutate(BBN.ID = str_replace(BBN.ID,"_","")) 
input$classfiles <- SQANTI_class_preparation(input_files$classfiles,"ns")
input$counts <- fread(input_files$counts, data.table = FALSE)

# remove outlier samples with 0 reads
input$ontPhenotype <- input$ontPhenotype %>% filter(!ID %in% c("BBN002.26311;A234/15","BBN003.26927;20150365","BBN_9944;A329/10")) 

Apoe <- input$classfiles %>% filter(associated_gene == "APOE")
Apoe <- merge(Apoe, input$counts, by.x = "isoform", by.y = "id", all.x = T)
 
ApoeB1 <- Apoe %>% dplyr::select(contains("B1")) %>% 
  colSums(.>0) %>% 
  reshape2::melt(value.name = "num_Apoe_transcripts") %>% 
  mutate(sample = str_remove(row.names(.), "B1.")) %>%
  merge(., input$ontPhenotype, by = "sample")

ApoeB2 <- Apoe %>% dplyr::select(contains("B2")) %>% 
  colSums(.>0) %>% 
  reshape2::melt(value.name = "num_Apoe_transcripts") %>% 
  mutate(sample = str_remove(row.names(.), "B2.")) %>%
  merge(., input$ontPhenotype, by = "sample")

ggplot(ApoeB1, aes(x = RINe, y = num_Apoe_transcripts, colour = APOE_geno)) + 
  geom_point(size = 2) + 
  scale_y_log10() + 
  theme_classic() + 
  labs(x = "RIN", y = "Number of transcripts")

ggplot(ApoeB2, aes(x = RINe, y = num_Apoe_transcripts, colour = APOE_geno)) + 
  geom_point(size = 2) + 
  scale_y_log10() + 
  theme_classic() + 
  labs(x = "RIN", y = "Number of transcripts")

cor.test(ApoeB2$RINe, ApoeB2$num_Apoe_transcripts)
```
        Pearson's product-moment correlation

data:  ApoeB2$RINe and ApoeB2$num_Apoe_transcripts
t = 3.0716, df = 48, p-value = 0.003502
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.1430931 0.6143397
sample estimates:
     cor
0.405298
```





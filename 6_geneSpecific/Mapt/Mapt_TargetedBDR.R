#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## output a list of transcripts that are documented in the TALON gtf but not in the abundance file
## sanity check for TALON pipeline
## Input:
##  --gtf = TALON gtf generated from talon_create_GTF
##  --counts = TALON abundance file
##  --output = output path and name of file
## --------------------------------


## ------- packages ---------------- 

library("dplyr")
library("ggrepel")

LOGEN_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
sapply(list.files(path = paste0(LOGEN_ROOT,"target_gene_annotation"), pattern="*summarise*", full = T), source,.GlobalEnv)
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN_ROOT, "transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "merge_characterise_dataset/run_ggtranscript.R"))


## ---------- input -----------------

dirnames <- list(
  bdr = "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/",
  ficle = "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/D_ONT/5_cupcake/8_characterise/TargetGenes/TargetGenes"
)

# read input files
input_files <- list(
  ontPhenotype = paste0(dirnames$bdr, "0_metadata/B_ONT/Selected_ONTTargeted_BDR.csv"), 
  # merged data (unfiltered but contain only target genes)
  classfiles = paste0(dirnames$bdr, "D_ONT/5_cupcake/7_sqanti3/ontBDR_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt"),
  # gtf
  gtf = paste0(dirnames$bdr, "D_ONT/5_cupcake/7_sqanti3/ontBDR_collapsed.filtered_counts_filtered.gtf"),
  # ref gtf
  refgtf = paste0(dirnames$bdr, "E_MergedTargeted/4_characterise/TargetGenesRef/MAPT_gencode.gtf"),
  # differential expression analysis
  ONTDeseq = paste0(dirnames$bdr, "01_figures_tables/Ont_DESeq2TranscriptLevel.RDS"),
  # mapt classification 
  maptClassification = paste0(dirnames$ficle,"/MAPT/Stats/MAPT_further_classifications.csv")
)

input <- list()
input$ontPhenotype <- read.table(input_files$ontPhenotype, sep = ",", header = T) %>% mutate(BBN.ID = str_replace(BBN.ID,"_","")) 
input$classfiles <- SQANTI_class_preparation(input_files$classfiles,"ns")
input$gtf <- as.data.frame(rtracklayer::import(input_files$gtf))
input$refgtf <- as.data.frame(rtracklayer::import(input_files$refgtf))
input$mergedgtf <- rbind(input$gtf[,c("seqnames","strand","start","end","type","transcript_id","gene_id")] ,
                         input$refgtf[,c("seqnames","strand","start","end","type","transcript_id","gene_id")])
input$ONTDeseq <- readRDS(input_files$ONTDeseq)
input$maptClassification <- read.csv(input_files$maptClassification)

# remove outlier samples with 0 reads
input$ontPhenotype <- input$ontPhenotype %>% filter(!ID %in% c("BBN002.26311;A234/15","BBN003.26927;20150365","BBN_9944;A329/10")) 

# import mapt classifications 0N4R
maptClassificationList <- lapply(unique(input$maptClassification$MAPT_classification), function(x) 
  input$maptClassification[input$maptClassification$MAPT_classification == x, "isoform"])
names(maptClassificationList) <- unique(input$maptClassification$MAPT_classification)
maptClassificationList <- 

## ---------- Mapt Dendrogram  -----------------

# mapt exon 2 - updated exon 5
# mapt exon 3 - updated exon 6
# mapt exon 10 - updated exon 16

plot_dendro_Tgene(dirnames$ficle, "MAPT")


## ---------- Mapt classifications  -----------------

# 4R visualisations
IsoList <- data.frame(
  Isoform = unlist(IsoList <- list(
    Reference = c("ENST00000344290.10","ENST00000415613.6"),
    `0N4R` = as.character(maptClassificationList$`0N4R`),
    `1N4R` = as.character(maptClassificationList$`1N4R`),
    `2N4R` = as.character(maptClassificationList$`2N4R`)
    
  )),
  Category = rep(names(IsoList), lengths(IsoList))
)
IsoList$colour <- c(rep(NA,length(IsoList$Category[IsoList$Category != "DTE"])))
Mapt4R <- ggTranPlots(input$mergedgtf, input$classfiles,
            isoList = c(as.character(IsoList$Isoform)),
            selfDf = IsoList, gene = "MAPT", squish=FALSE)


# 3R visualisations
IsoList <- data.frame(
  Isoform = unlist(IsoList <- list(
    Reference = c("ENST00000344290.10","ENST00000415613.6"),
    `0N3R` = as.character(maptClassificationList$`0N3R`),
    `1N3R` = as.character(maptClassificationList$`1N3R`),
    `2N3R` = as.character(maptClassificationList$`2N3R`)
  )),
  Category = rep(names(IsoList), lengths(IsoList))
)
IsoList$colour <- c(rep(NA,length(IsoList$Category[IsoList$Category != "DTE"])))
Mapt3R <- ggTranPlots(input$mergedgtf, input$classfiles,
            isoList = c(as.character(IsoList$Isoform)),
            selfDf = IsoList, gene = "MAPT",squish=FALSE)

plot_grid(Mapt3R,Mapt4R)


## ---------- Mapt classifications  -----------------

# summarise all subsets of 4R and 3R isoforms
maptClassificationListSimplified <- list(
  R4 = as.character(unlist(list(maptClassificationList$`0N4R`,maptClassificationList$`1N4R`,maptClassificationList$`2N4R`))),
  R3 = as.character(unlist(list(maptClassificationList$`0N3R`,maptClassificationList$`1N3R`,maptClassificationList$`2N3R`)))
)

# sum the normalized counts from 4R and 3R isoforms
sum_mapt_classification <- function(lst, category){
  dat <- input$ONTDeseq$B2Wald$norm_counts %>% filter(isoform %in% lst[[category]]) %>% 
    group_by(sample) %>% tally(normalised_counts) 
  colnames(dat) <- c("sample", category)
  
  return(dat)
}


# plot ratio
# note the black dot with the higest value is a FLTD sample
do.call(cbind,lapply(names(maptClassificationListSimplified), function(x) sum_mapt_classification(maptClassificationListSimplified, x))) %>% 
  tibble::column_to_rownames(., var = "sample") %>%
  mutate(sample = str_remove(row.names(.), "B2.")) %>% 
  mutate(ratio = R4/R3) %>% 
  merge(ontPhenotype, ., by = "sample") %>% 
  filter(BraakTangle_numeric %in% c(0,1,2,5,6)) %>%
  mutate(phenotype = factor(ifelse(BraakTangle_numeric %in% c(0,1,2), "control", "AD"),levels = c("control","AD"))) %>%
  ggplot(., aes(x = phenotype, y = ratio, colour = phenotype)) + geom_boxplot(outliers = FALSE) +
  geom_jitter(size=2, alpha=0.9) +
  labs(x = "Phenotype", y = "Ratio 4R: 3R") +
  theme_classic() +
  scale_colour_manual(values = c("black","red")) +
  theme(legend.position='none')

# plot individual subsets
do.call(cbind,lapply(names(maptClassificationList), function(x) sum_mapt_classification(maptClassificationList, x))) %>%
  tibble::column_to_rownames(., var = "sample") %>%
  select(-sample) %>%
  mutate(sample = str_remove(row.names(.), "B2.")) %>% 
  reshape2::melt(., id = "sample", variable.name = "category", value.name = "counts") %>% 
  merge(ontPhenotype, ., by = "sample") %>% 
  filter(BraakTangle_numeric %in% c(0,1,2,5,6)) %>%
  mutate(phenotype = factor(ifelse(BraakTangle_numeric %in% c(0,1,2), "control", "AD"),levels = c("control","AD"))) %>%
  ggplot(., aes(x = phenotype, y = counts, colour = phenotype)) + geom_boxplot(outliers = FALSE) +
  facet_grid(~category) +
  geom_jitter(size=2, alpha=0.9) +
  labs(x = "Phenotype", y = "Normalised counts") +
  theme_classic() +
  scale_colour_manual(values = c("black","red")) +
  theme(legend.position='none')


# plot isoform fraction for each sample by 4R vs 3R
dat <- input$ONTDeseq$B2WaldBraak$norm_counts %>% filter(associated_gene == "MAPT") %>% 
  select(sample, isoform, normalised_counts) %>% 
  spread(., key = sample, value = normalised_counts) %>% 
  merge(., input$maptClassification, by = "isoform") %>%
  select(-isoform)
dat <- aggregate(. ~ MAPT_classification, dat, sum) %>% tibble::column_to_rownames(., var = "MAPT_classification")

dat <- as.data.frame(apply(dat,2,function(x){x/sum(x) * 100})) %>%
  tibble::rownames_to_column(., var = "MAPT_classification") %>%
  reshape2::melt(variable.name = "sample", value.name = "sum_counts") %>%
  mutate(sample = str_remove(sample, "B2.")) 


dat <- merge(ontPhenotype, dat, by = "sample") %>% 
  filter(BraakTangle_numeric %in% c(0,1,2,5,6)) %>%
  mutate(phenotype = factor(ifelse(BraakTangle_numeric %in% c(0,1,2), "control", "AD"),levels = c("control","AD")))

identify_overallMapt_class <- function(classification){
  
  if(classification %in% c("0N3R", "1N3R", "2N3R")){return("3R")
  }else if(classification %in% c("0N4R", "1N4R", "2N4R")){return("4R")
  }else{return("Truncated")}
  
}

dat$MaptClassificationSum  <- apply(dat, 1, function(x) identify_overallMapt_class(x[["MAPT_classification"]]))

# Plotting
dat %>% mutate(MAPT_classification = factor(MAPT_classification, levels = c("0N4R", "1N4R", "2N4R", "0N3R","1N3R","2N3R", "E10Truncated", "E2E3'3R", "E2E3'4R"))) %>% 
  ggplot(., aes(x = sample, y = sum_counts, fill = MAPT_classification)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("#FFC0CB", "#FF69B4", "#FF0000", "#D3D3D3", "#808080", "#000000", "#ADD8E6", "#0000FF", "#00008B")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(~ phenotype, scale = "free") +
  labs(x = "Sample", y = "Isoform Fraction")


## ---------- Mapt transcripts specific to AD  -----------------

ontPhenotype = read.csv(paste0("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/0_metadata/B_ONT/Selected_ONTTargeted_BDR.csv"))
phenotype96 = read.csv("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/0_metadata/96_samples_v2.csv")
phenotype96 <- phenotype96 %>% mutate(ID = paste0(BBN.ID,";",Brain.ID)) %>% 
  select(ID, BBN.ID, Brain.ID, Institute, Age.x, Sex.x, Braak_tangle, APOE_geno, Thal_amyloid, cerad_densitiy, AD, FAD, LB_stage, LB_type, LBD, FTLD, 
         severe_CAA, severe_SVD, Vascularity, parkinsons, TDP43, Dementia, S.C, AGONALST...link.to.death) 
ontPhenotype <- merge(phenotype96, ontPhenotype[,c("ID","sample")], by = "ID")


# mapt transcripts specific to AD
MaptCounts <- input$ONTDeseq$B2WaldBraak$norm_counts[input$ONTDeseq$B2WaldBraak$norm_counts$associated_gene == "MAPT",]
MaptCounts <- MaptCounts[grepl("B2", MaptCounts$sample),] %>% mutate(sample = word(sample,c(2),sep=fixed(".")))
MaptCounts <- MaptCounts %>% mutate(sample = word(sample,c(1),sep=fixed(".")))
ontPhenotype <- ontPhenotype %>% mutate(sample = word(sample,c(1),sep=fixed(".")))

# remove samples that have LB stage
MaptCounts <- merge(MaptCounts, ontPhenotype[,c("sample","APOE_geno","LBD","Braak_tangle")], by = "sample")
MaptCounts <- MaptCounts %>% filter(Braak_tangle %in% c(0,1,2,5,6)) %>% mutate(ADByBraak = ifelse(Braak_tangle %in% c(0,1,2), "control", "AD"))
MaptCountsGroup <- MaptCounts %>% group_by(isoform, ADByBraak) %>% tally(normalised_counts) 
notDetectedLowBraak <- MaptCountsGroup[MaptCountsGroup$ADByBraak == "control" & MaptCountsGroup$n == 0,]

dat <- input$classfiles[input$classfiles$isoform == "PB.39281.399",] %>% select(contains("B2")) %>% 
  reshape2::melt(variable.name = "sample", value.name = "FLReads") %>% 
  mutate(sample = word(sample,c(2),sep=fixed("."))) %>%
  merge(., ontPhenotype[,c("sample","APOE_geno","LBD","Braak_tangle")], by = "sample") %>%
  filter(Braak_tangle %in% c(0,1,2,5,6)) %>% 
  mutate(ADByBraak = ifelse(Braak_tangle %in% c(0,1,2), "control", "AD")) 

p1 <- ggplot(dat, aes(x = as.factor(Braak_tangle), y = FLReads, colour = as.factor(LBD), label = as.factor(sample))) + geom_point() +
  geom_label_repel(aes(label=ifelse(FLReads>0,as.character(sample),'')),hjust=0.5,vjust=0, show_guide = F) +
  theme_classic() +
  labs(x = "Braak stage", y = "FL Reads", colour = "LBD") 

PB.39281.399Sample <- ontPhenotype[ontPhenotype$sample %in% dat[dat$ADByBraak == "AD" & dat$FLReads > 0,"sample"],]

p2 <- ggTranPlots(input$mergedgtf, input$classfiles,
            isoList = c("ENST00000344290.10", "PB.39281.399","PB.39281.4292","PB.39281.4685"),
            gene = "MAPT",simple = TRUE,colours = c("blue","red","red","red"))


plot_grid(p2, p1, ncol = 1, rel_heights = c(0.3,0.7))
# protein
#cd /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/C_Proteomics/ONT/9_metamorpheus/filtered/ADBDR_filtered_search_results/Task1SearchTask/
#header=$(head -n 1 AllPeptides.filtered.psmtsv)
#grep "PB.39281.399" AllPeptides.filtered.psmtsv > PB.39281.399.filtered.psmtsv
#echo "$header" | cat - PB.39281.399.filtered.psmtsv > temp && mv temp PB.39281.399.filtered.psmtsv


proteinPeptides <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/C_Proteomics/ONT/9_metamorpheus/filtered/ADBDR_filtered_search_results/Task1SearchTask/PB.39281.399.filtered.psmtsv", sep = "\t", header = T)

# fetal dataset
sqantiDir <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/6_sqanti/sqanti/"
FetalclassFilesName <- paste0(sqantiDir, "WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification_noMonoIntergenic.txt")
classFiles <- SQANTI_class_preparation(FetalclassFilesName,"ns")
maptClassFilesmultiNIC <- classFiles[classFiles$associated_gene == "MAPT" & classFiles$structural_category == "NIC" & classFiles$exons != 1,]

tauMAPT = read.csv("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/devFicle/TargetGenes/MAPT/Stats/MAPT_general_exon_level.csv")
write.table(tauMAPT[tauMAPT$isoform %in% maptClassFilesmultiNIC$isoform,"isoform"],paste0(sqantiDir,"MAPTNIC.txt"),quote=F,row.names=F)
grep -w -f MAPTNIC.txt WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_Whole_2reads2samples_nomonointergenic.gtf > MaptNIC_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_Whole_2reads2samples_nomonointergenic.gtf

# PacBio dataset
pacbioPhenotype <- read.csv("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/0_metadata/A_IsoSeq/Selected_BDR_13012021.csv")
PacBioClassFile <- SQANTI_class_preparation("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/9_sqanti3/basic/AllBDRTargeted.collapsed_classification.filtered_lite_classification.txt", "ns")
#write.table(PacBioClassFile[PacBioClassFile$associated_gene == "MAPT" & PacBioClassFile$structural_category == "NIC","isoform"],
#            "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/9_sqanti3/basic/MAPTNIC.txt",quote=F,row.names=F)
#grep -w -f MAPTNIC.txt AllBDRTargeted.collapsed_corrected.gtf > MAPTNIC.gtf

# PB.4515.49, PB.4515.46
PacBioClassFile[PacBioClassFile$isoform == "PB.4515.46",] %>% select(isoform, contains("FL.")) %>% 
  reshape2::melt(id = "isoform") %>% 
  mutate(barcode = stringr::word(variable, c(2), sep = fixed("."))) %>%
  merge(., pacbioPhenotype, by.y = "Sample", by.x = "barcode") %>%
  ggplot(., aes(x = Braak, y = value)) + geom_point() +
  theme_classic() +
  labs(x = "Braak stage", y = "FL Reads") 


##---------------- novel exon

input$classfiles[input$classfiles$isoform == "PB.39281.2140",] %>% select(contains("B2")) %>% 
  reshape2::melt(variable.name = "sample", value.name = "FLReads") %>% 
  mutate(sample = word(sample,c(2),sep=fixed("."))) %>%
  merge(., ontPhenotype[,c("sample","APOE_geno","LBD","Braak_tangle")], by = "sample") %>%
  filter(Braak_tangle %in% c(0,1,2,5,6)) %>% 
  mutate(ADByBraak = ifelse(Braak_tangle %in% c(0,1,2), "control", "AD")) %>%
  ggplot(., aes(x = as.factor(Braak_tangle), y = FLReads, colour = as.factor(LBD), label = as.factor(sample))) + geom_point() +
  geom_label_repel(aes(label=ifelse(FLReads>0,as.character(sample),'')),hjust=0.5,vjust=0, show_guide = F) +
  theme_classic() +
  labs(x = "Braak stage", y = "FL Reads", colour = "LBD") 
            
            
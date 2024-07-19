LOGEN_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
sapply(list.files(path = paste0(LOGEN_ROOT,"target_gene_annotation"), pattern="*summarise*", full = T), source,.GlobalEnv)
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN_ROOT, "transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "merge_characterise_dataset/run_ggtranscript.R"))

classFilesName <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/6_sqanti/sqanti/WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification_noMonoIntergenic.txt"

classFiles <- SQANTI_class_preparation(classFilesName,"ns")
maptClassFilesmulti <- classFiles[classFiles$associated_gene == "MAPT" & classFiles$exons != 1,]
write.table(maptClassFilesmulti %>% select(-"V1"), "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/devFicle/Mapt_WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification_noMonoIntergenic.txt", row.names = F, quote = F, sep = "\t")

ontPhenotype = read.csv(paste0("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/0_metadata/WholeTargetedphenotype_fixedsex.csv"))

# gtf files
inputGtf = as.data.frame(rtracklayer::import("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/6_sqanti/sqanti/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_Whole_2reads2samples_nomonointergenic.gtf"))
refGtf = as.data.frame(rtracklayer::import("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/E_MergedTargeted/4_characterise/TargetGenesRef/MAPT_gencode.gtf"))
mergedGtf <- rbind(inputGtf[,c("seqnames","strand","start","end","type","transcript_id","gene_id")] ,
                   refGtf[,c("seqnames","strand","start","end","type","transcript_id","gene_id")])

ONTDeseq <- readRDS("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/01_figures_tables/Ont_DESeq2TranscriptLevel.RDS")
counts <- read.csv("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/10_deseq/2_DTE/MAPT_DESeq2_whole_development_normAll.csv", header = F)
colnames(counts) <- c("isoform","sample","normalised_counts","associated_gene","associated_transcript","structural_category","subcategory","exons")


# mapt exon 2 - updated exon 5
# mapt exon 3 - updated exon 6
# mapt exon 10 - updated exon 16

# dendrogram 
ficleDir = "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/devFicle/TargetGenes/"
plot_dendro_Tgene(ficleDir, "MAPT")

maptClassification <- read.csv(paste0(ficleDir,"/MAPT/Stats/MAPT_further_classifications.csv"))
maptClassificationList <- lapply(unique(maptClassification$MAPT_classification), function(x) 
  maptClassification[maptClassification$MAPT_classification == x, "isoform"])
names(maptClassificationList) <- unique(maptClassification$MAPT_classification)

IsoList <- data.frame(
  Isoform = unlist(IsoList <- list(
    #Reference = unique(refGtf[refGtf$gene_name == "MAPT" & !is.na(refGtf$transcript_id), "transcript_id"]),
    Reference = c("ENST00000344290.10","ENST00000415613.6"),
    `0N4R` = as.character(maptClassificationList$`0N4R`),
    `1N4R` = as.character(maptClassificationList$`1N4R`),
    `2N4R` = as.character(maptClassificationList$`2N4R`)
    
  )),
  Category = rep(names(IsoList), lengths(IsoList))
)
IsoList$colour <- c(rep(NA,length(IsoList$Category[IsoList$Category != "DTE"])))

Mapt4R <- ggTranPlots(mergedGtf, classFiles,
                      isoList = c(as.character(IsoList$Isoform)),
                      selfDf = IsoList, gene = "MAPT", squish=TRUE)


IsoList <- data.frame(
  Isoform = unlist(IsoList <- list(
    #Reference = unique(refGtf[refGtf$gene_name == "MAPT" & !is.na(refGtf$transcript_id), "transcript_id"]),
    Reference = c("ENST00000344290.10","ENST00000415613.6"),
    `0N3R` = as.character(maptClassificationList$`0N3R`),
    `1N3R` = as.character(maptClassificationList$`1N3R`),
    `2N3R` = as.character(maptClassificationList$`2N3R`)
  )),
  Category = rep(names(IsoList), lengths(IsoList))
)
IsoList$colour <- c(rep(NA,length(IsoList$Category[IsoList$Category != "DTE"])))
Mapt3R <- ggTranPlots(mergedGtf, classFiles,
                      isoList = c(as.character(IsoList$Isoform)),
                      selfDf = IsoList, gene = "MAPT",squish=TRUE)


plot_grid(Mapt3R,Mapt4R)

## 


maptClassificationListSimplified <- list(
  R4 = as.character(unlist(list(maptClassificationList$`0N4R`,maptClassificationList$`1N4R`,maptClassificationList$`2N4R`))),
  R3 = as.character(unlist(list(maptClassificationList$`0N3R`,maptClassificationList$`1N3R`,maptClassificationList$`2N3R`)))
)

sum_mapt_classification <- function(category, breakdown=FALSE){
  
  if(isTRUE(breakdown)){
    dat <- counts %>% filter(isoform %in% maptClassificationList[[category]]) %>% 
      group_by(sample) %>% tally(normalised_counts) 
    
  }else{
    dat <- counts %>% filter(isoform %in% maptClassificationListSimplified[[category]]) %>% 
      group_by(sample) %>% tally(normalised_counts)  
  }
  colnames(dat) <- c("sample", category)
  
  return(dat)
}


dat <- do.call(cbind,lapply(names(maptClassificationListSimplified), function(x) sum_mapt_classification(x))) %>% 
  tibble::column_to_rownames(., var = "sample") %>%
  mutate(sample = str_remove(row.names(.), "B2.")) %>% 
  mutate(ratio = R4/R3)
dat <- merge(ontPhenotype, dat, by = "sample") 
ggplot(dat, aes(x = group, y = ratio, colour = group)) + geom_boxplot(outliers = FALSE) +
  geom_jitter(size=2, alpha=0.9) +
  labs(x = "Phenotype", y = "Ratio 4R: 3R") +
  theme_classic() +
  scale_colour_manual(values = c("black","red")) +
  theme(legend.position='none')


lst <- lapply(names(maptClassificationList), function(x) sum_mapt_classification(x, breakdown = TRUE))
dat <- purrr::reduce(lst, full_join, by = "sample")
dat <- merge(ontPhenotype[,c("sample","group")], dat, by = "sample") 
dat <- reshape2::melt(dat, variable.name = "category", value.name = "counts") 
dat <- dat %>% mutate(broadCate = ifelse(grepl("3R", category), "3R","4R"))
ggplot(dat, aes(x = category, y = counts, colour = group)) + geom_boxplot(outliers = FALSE) +
  facet_grid(~broadCate, scales = "free") +
  theme_classic() +
  labs(x = "Category", y = "Normalised counts")

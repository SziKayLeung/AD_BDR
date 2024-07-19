library("dplyr")
library("stringr")

phenotype <- read.csv("/lustre/projects/Research_Project-MRC193462/BDR/Phenotypes/pathology/BDR_pathology_data_for_analysis.csv")

NeuroChipVariants <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/genotype/NeuroChipDiseaseVariants.txt", sep = "\t", 
                                header = T)

rareVariants <- read.csv("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/genotype/AnnotatedRareVariantsPresentinBDR.csv")

variant <- list.files(path = "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/genotype/imputed", pattern = "Variant_Phenotype.txt", 
           recursive = T, full.names = T)
variant <- lapply(variant, function(x) read.table(x, sep = "\t", header = T))
variantMerge <- bind_cols(variant)
variantMerge <- variantMerge %>% select("BBNId...1", "Brain_ID...2", "DNA_ID...3", contains("X"))
colnames(variantMerge) <- stringr::str_remove(colnames(variantMerge), "X")

setdiff(colnames(variantMerge),paste0("X", NeuroChipVariants$NeuroChip_variant_location_hg19))

variantsDetected <- word(colnames(variantMerge), c(1), sep = fixed("_"))
variantsDetected <- str_replace_all(variantsDetected, "\\.", ":")

NeuroChipVariants[NeuroChipVariants$NeuroChip_variant_location_hg19 %in% variantsDetected,]
NeuroChipVariants[NeuroChipVariants$HGMD_dbSNP_rsID %in% unique(word(rareVariants$paste.colnames.raw..i....1...sep...._.,c(1),sep=fixed("_"))),]

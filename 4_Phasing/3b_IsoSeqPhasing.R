library("dplyr")
library("stringr")

LOGEN_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "/transcriptome_stats/read_sq_classification.R"))


IsoPhaseVcf <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/IsoSeq_IsoPhase.vcf")
phasingSummary <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/summarized.isophase_results.txt", header = T)
class.names.files <- "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/9_sqanti3/basic/AllBDRTargeted.collapsed_classification.filtered_lite_classification.txt"
class.file <- SQANTI_class_preparation(class.names.files,"ns")
class.file <- class.file %>% mutate(associated_gene_id = paste0("PB.", word(isoform,c(2),sep=fixed("."))))

head(class.file)

phasingSummary <- phasingSummary %>% mutate(associated_gene_id = word(word(locus,c(2), sep = fixed("/")),c(1), sep = fixed("_")))
phasingSummary <- merge(phasingSummary, class.file[,c("associated_gene_id","associated_gene")], by = "associated_gene_id")
phasingSummary <- distinct(phasingSummary)

# final tagged files
phasingTagged <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/finalTagged.txt", header = F, col.names = c("locus"))
phasingTagged <- phasingTagged  %>% mutate(associated_gene_id = word(locus,c(1), sep = fixed("_")))
phasingTagged <- merge(phasingTagged, class.file[,c("chrom","associated_gene_id","associated_gene")], by = "associated_gene_id")
phasingTagged <- distinct(phasingTagged)

TargetGene <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/0_metadata/B_ONT/TargetGenes.tsv")[["V1"]]
phasingTagged[phasingTagged$associated_gene %in% TargetGene,]

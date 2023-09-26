library("dplyr")
library("tibble")
library(DRIMSeq)
library(DEXSeq)


input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR"
# Read in SQANTI2 classification file of all merged data
class <- read.table(paste0(input_dir,"/Post_IsoSeq/SQANTI_TAMA_FILTER/AllBDRTargeted_sqantitamafiltered.classification.txt"),header=T, as.is = T, sep = "\t")

TargetGene <- c("ABCA1","SORL1","MAPT","BIN1","TARDBP","APP","ABCA7","PTK2B","ANK1","FYN","CLU","CD33","FUS","PICALM","SNCA","APOE","TRPA1","RHBDF2","TREM2","VGF")

class <- class %>% filter(toupper(associated_gene) %in% TargetGene)


# phenotype dataset 
pheno <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/AD_BDR/Raw_Data/ADBDR_Phenotype.csv"
coldata  <- read.csv(pheno, header = T) %>% `colnames<-`(c("sample_id", "condition","Age","Sex"))
samps <- coldata
samps$condition <- factor(samps$condition)
table(samps$condition)


# datawrangle for input to DESeq2 
class[["transcript_name_id"]] <- paste0(class$associated_transcript,"_", class$isoform)
counts  <- class %>% dplyr::select(transcript_name_id,associated_gene,starts_with("FL.")) 
colnames(counts)[1] <- "feature_id"
colnames(counts)[2] <- "gene_id"

counts  <- class %>% dplyr::select(isoform,starts_with("FL.")) 
write.table(counts, "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq/TAPPAS/counts.txt", quote = F, row.names = F,sep = "\t")

d <- dmDSdata(counts=counts, samples=samps)
trs_cts <- counts(d)
n <- 30
n.small <-5
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=10)
d
table(table(counts(d)$gene_id))
design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
colnames(design_full)


sample.data <- DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
sample.data$condition <- relevel(sample.data$condition, ref = "Control")

dxd <- DEXSeqDataSet(countData=count.data, sampleData=sample.data, design=~sample + exon + condition:exon, featureID=trs_cts$feature_id, groupID=trs_cts$gene_id)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd, reducedModel=~sample + exon)
dxd <- estimateExonFoldChanges( dxd, fitExpToVar="condition")
dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
plotMA(dxr, cex=0.8, alpha=0.05) 
plotDispEsts(dxd)
qval <- perGeneQValue(dxr) 
dxr.g<-data.frame(gene=names(qval), qval)
dxr.g <- dxr.g[order(dxr.g$qval),]
dxr_out <- as.data.frame(dxr[,c("featureID", "groupID", "pvalue")])
dxr_out <- dxr_out[order(dxr$pvalue),]

output_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Differential"
write.table(dxr.g, file=paste0(output_dir,"/results_dtu_gene.tsv"), sep="\t")
write.table(dxr_out, file=paste0(output_dir,"/results_dtu_transcript.tsv"), sep="\t")

colnames(dxr)[grep("log2fold", colnames(dxr))] <- "log2fold"
MADTUdata <- data.frame(dxr)[order(dxr$padj),c("exonBaseMean", "log2fold", "pvalue", "padj")]
MADTUdata$exonBaseMean <- log2(MADTUdata$exonBaseMean)
colnames(MADTUdata)[which(colnames(MADTUdata)=="exonBaseMean")] <- "Log2MeanExon"
colnames(MADTUdata)[which(colnames(MADTUdata)=="log2fold")] <- "Log2FC"
write.table(MADTUdata, file=paste0(output_dir,"/results_dexseq.tsv"), sep="\t")


cat("stageR analysis\n")
library(stageR)

cat("Running stageR analysis on the differential transcript usage results.\n")
pConfirmation <- matrix(dxr$pvalue, ncol=1)
dimnames(pConfirmation) <- list(dxr$featureID, "transcript")
pScreen <- qval
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)
# note: the choice of 0.05 here means you can *only* threshold at 5% OFDR later
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.25)
suppressWarnings({dex.padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE)})

write.table(dex.padj, file=paste0(output_dir,"/results_dtu_stageR.tsv"), sep="\t")

dxr1 = DEXSeqResults( dxd )
dxr1
mcols(dxr1)$description
table ( dxr1$padj < 0.1 )
table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) )
plotMA( dxr1, cex=0.8 )
plotDEXSeq(dxr1, "FYN")


# drimseq
library(DRIMSeq)
d <- dmDSdata(counts=counts, samples=samps)
d
n <- 12
n.small <- 6
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=10)
d       
table(table(counts(d)$gene_id))
design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
colnames(design_full)
set.seed(1)
system.time({
  d <- dmPrecision(d, design=design_full)
  d <- dmFit(d, design=design_full)
  d <- dmTest(d, coef="conditionControl")
})
res <- DRIMSeq::results(d)
head(res)
res.txp <- DRIMSeq::results(d, level="feature")
head(res.txp)
no.na <- function(x) ifelse(is.na(x), 1, x)
res$pvalue <- no.na(res$pvalue)
res.txp$pvalue <- no.na(res.txp$pvalue)
idx <- which(res$adj_pvalue < 0.05)[1]
res[idx,]
plotProportions(d, res$gene_id[idx], "condition")


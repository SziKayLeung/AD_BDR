#### edgeR Differential Expression - 'data' contains raw counts
suppressMessages(library(edgeR))

indir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq/TAPPAS/TAPPAS_output"

phenotype = read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/AD_BDR/Raw_Data/ADBDR_Phenotype.csv", header = T) 

de_edgeR <- function(data, factors) {
  tmm_factors = calcNormFactors(round(data))
  myedgeR = DGEList(counts = round(data),
                    group = as.factor(factors[,1]),
                    norm.factors = tmm_factors, 
                    remove.zeros = FALSE)
  myedgeR = estimateDisp(myedgeR)
  res = exactTest(myedgeR, pair=sort(levels(as.factor(myfactors[,1])), decreasing = T))
  res.sort = topTags(res, n = nrow(data))[[1]]
  summary(de <- decideTestsDGE(res, p=0.05, adjust="BH"))
  return(res.sort)
}

de_edgeR_glm <- function(data, factors) {
  tmm_factors = calcNormFactors(round(data))
  myedgeR = DGEList(counts = round(data),
                    group = as.factor(factors[,1]),
                    norm.factors = tmm_factors, 
                    remove.zeros = FALSE)
  
  #design <- model.matrix(~0+group, data=myedgeR$samples)
  #colnames(design) <- c("Ctrl","AD")
  
  # Covariate Sex
  design <- model.matrix(~0+group+Sex, data=myedgeR$samples)
  colnames(design) <- c("Ctrl","AD","Sex")
  #lrt <- glmLRT(fit, coef=2:3)
  
  # data exploration by phenotype 
  #plotMDS(myedgeR, col = c(rep(label_colour("AD"),15), rep(label_colour("Control"),15)))
  
  # data exploration by sex
  #Sexcol = recode(Sex, female = wes_palette("Royal2")[3], male = wes_palette("Royal2")[5])
  #plotMDS(myedgeR, col = paste0('', Sexcol, ''))
  
  
  myedgeR = estimateDisp(myedgeR,design)
  fit <- glmFit(myedgeR, design)
  ctrlvsAD <- makeContrasts(Ctrl-AD, levels=design)
  lrt <- glmLRT(fit,contrast=ctrlvsAD)
  lrt$table$FDR <- p.adjust(lrt$table$PValue, method="BH")
  
  
  # qlf 
  #qlf<-glmQLFit(myedgeR, design)
  #qlt<-glmQLFTest(qlf, contrast = ctrlvsAD)
  #deg<-as.data.frame(topTags(qlt, n=Inf))
  
  return(lrt$table)
}



#### Read expression factors

cat("\nReading factors file data...")
myfactors=read.table(file.path(indir, "exp_factors.txt"), row.names=1, sep="\t", header=TRUE)
#myfactors = merge(myfactors,phenotype, by.x = 0,by.y = "sampleID")

run_de_edgeR <- function(type, method){
  cat("\nReading raw counts",type,"matrix file data...")
  expMatrix = read.table(file.path(paste0(indir,"/", type,"_matrix_raw.tsv")), row.names=1, sep="\t", header=TRUE)
  cat("\nRead ", nrow(expMatrix), " raw counts",type,"expression data rows")
  
  # Create Sex category 
  datapheno = expMatrix %>% reshape2::melt() %>% 
    full_join(., phenotype, by = c("variable"= "sampleID")) %>% .[,c("variable","condition","Age","Sex")] %>% distinct()
  Sex <- factor(datapheno$Sex)
  
  if(method == "exact_test"){
    result = de_edgeR(data=expMatrix, factors=myfactors)
  } else if(method == "glm"){
    result = de_edgeR_glm(data=expMatrix, factors=myfactors)
  } else {
    print("exact_test or glm required")
  }
  
  return(result)
}

result_trans = run_de_edgeR("transcript","exact_test")
result_gene = run_de_edgeR("gene","exact_test")

result2_gene = run_de_edgeR("gene","glm")
result2_trans = run_de_edgeR("transcript","glm")







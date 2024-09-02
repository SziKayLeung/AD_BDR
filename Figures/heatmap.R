suppressMessages(library(pheatmap))
suppressMessages(library(viridis)) 

normCounts = annoResTran$ontMerged$norm_counts %>% filter(isoform %in% TargetedDESeq$ontResTranAnno$ontMerged$anno_res$isoform)
cf = class.files$targ_ont_filtered %>% filter(isoform %in% normCounts$isoform)
pheno = phenotype$targ_ont

dat = normCounts[,c("isoform","normalised_counts","sample")] %>%
  mutate(log2normalised = log2(normalised_counts)) %>% 
  select(isoform, log2normalised, sample) %>%
  spread(., isoform, log2normalised) %>% tibble::column_to_rownames(var = "sample")
dat[dat == "-Inf"] <- 0
dat.t <- t(dat)

coldata = normCounts %>% 
  dplyr::select(sample, group) %>% distinct(.keep_all = TRUE) %>%
  mutate(sampleid = paste0(word(sample,c(2), sep = fixed(".")),".",word(sample,c(3), sep = fixed(".")))) %>% 
  merge(., pheno[,c("sample","Gender","APOE_geno")], by.x = "sampleid", by.y = "sample") %>%
  column_to_rownames(var = "sample") %>% select(-sampleid)  
colnames(coldata) = c("Phenotype","Gender","ApoeStatus")

rowdata = cf[cf$isoform %in% colnames(dat),c("isoform","structural_category")] %>% dplyr::select(-isoform)
colnames(rowdata) = c("Category")

annotation_colors = list(
  Phenotype = c("0" = "white","1"="white","2"="black","3"="black","4"="black","5"="red","6" = "red"), 
  Gender = c("male" = "blue","female" = "pink"),
  ApoeStatus = c("e2/e2" = "white","e2/e3" = "#3B9AB2","e2/e4" = "#78B7C5", "e3/e3" = "#EBCC2A", "e3/e4" = "#E1AF00", "e4/e4" = "#F21A00"),
  Category = c(FSM = alpha("#00BFC4",0.8),ISM = alpha("#00BFC4",0.3),NIC = alpha("#F8766D",0.8),NNC = alpha("#F8766D",0.3))
)

pheatmap(dat.t, 
         annotation_col = coldata, 
         annotation_row = rowdata, 
         annotation_legend = FALSE,
         show_colnames = FALSE, 
         show_rownames = TRUE, 
         color = viridis(10),
         annotation_colors = annotation_colors,
         fontsize_col = 20,
         labels_row = cf[match(rownames(dat.t), cf$isoform),"associated_gene"], 
         labels_col = FALSE,
         legend = FALSE,
         annotation_names_row = FALSE)

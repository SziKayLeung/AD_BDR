library("dplyr")
library("VennDiagram")
library("RColorBrewer")
class.files$SCN <- class.files$targ_merged %>% select(contains("SCN")) 
colnames(class.files$SCN) <- as.character(paste0(phenotype$scn$col[match(names(class.files$SCN), phenotype$scn$col)],
                                                 "_",
                                                 phenotype$scn$group[match(names(class.files$SCN), phenotype$scn$col)]))

class.files$targ_merged <- cbind(class.files$targ_merged %>% select(!contains("SCN")),class.files$SCN)
class.files$targ_merged <- cbind(class.files$targ_merged,
      cbind(
        class.files$targ_merged %>% select(contains("SCN")) %>% apply(., 1, sum),
        class.files$targ_merged %>% select(contains("Sox10")) %>% apply(., 1, sum),
        class.files$targ_merged %>% select(contains("NeuN")) %>% apply(., 1, sum)
      ) %>% `colnames<-`(c("SCNSum","Sox10Sum","NeuNSum"))
)
class.files$targ_merged$Dataset <- apply(class.files$targ_merged, 1, function(x) identify_dataset_by_counts (x[["Sox10Sum"]], x[["NeuNSum"]], "Sox10","NeuN"))


class.files$SCN <- class.files$targ_merged %>% select(associated_gene, associated_transcript, isoform, structural_category, contains("SCN")) %>%
  reshape2::melt(id=c("associated_gene","associated_transcript","isoform","structural_category")) %>%
  filter(variable != "sum") %>%
  mutate(celltype = word(variable, c(3),sep=fixed("_"))) %>% 
  filter(value >= 1) %>% as.data.frame()

# Number of transcripts in SCN dataset
p1 <- class.files$SCN %>% filter(!grepl("Sum", variable)) %>% 
  group_by(isoform,associated_gene, celltype) %>% tally() %>% 
  ggplot(.,aes(x = associated_gene, y = n, fill = celltype)) + geom_bar(stat = "identity") + 
  labs(x = "associated gene", y = "Number of transcripts") + mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


isoformByCellType <- list(
  Total = class.files$SCN %>% filter(celltype == "Total"),
  Sox10 = class.files$SCN %>% filter(celltype == "Sox10"),
  NeuN = class.files$SCN %>% filter(celltype == "NeuN")
)
p2 <- venn.diagram(
  x = list(isoformByCellType$Total$isoform, isoformByCellType$Sox10$isoform, isoformByCellType$NeuN$isoform),
  category.names = c("Total" , "Sox10" , "NeuN"),
  fill = brewer.pal(3, "Pastel2"),
  filename = NULL
)


dat <- class.files$targ_merged %>% select(associated_gene, associated_transcript, isoform, structural_category, Dataset, contains("SCN")) %>%
  reshape2::melt(id=c("associated_gene","associated_transcript","isoform","structural_category","Dataset")) %>% 
  filter(value >= 1, !grepl("Sum", variable)) %>%
  mutate(celltype = word(variable, c(3),sep=fixed("_"))) 

dat$Dataset<- replace(dat$Dataset, dat$Dataset == "NA", "Total")

p3 <- dat %>% group_by(isoform, associated_gene,Dataset,structural_category) %>% tally() %>%
  ggplot(., aes(x = associated_gene, y = n, fill = structural_category)) + geom_bar(stat = "identity") + 
  facet_grid(~Dataset,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = "Gene", y = "Number of transcripts") + mytheme +
  theme(legend.position = "top")


# visualization
dat[dat$Dataset %in% c("NeuN","Sox10") & dat$structural_category %in% c("NIC","NNC"),] %>% filter(associated_gene == "CLU")

CluIso <- data.frame(
  Isoform = unlist(CluIso <- list(
    Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == "CLU" & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
    NeuN = as.character(c("PB.293270.2073","PB.293270.29996","PB.293270.2266","PB.293270.2193","PB.293270.2077","PB.293270.2458")),
    Sox10 = as.character(c("PB.293270.2397","PB.293270.2993")))),
  Category = rep(names(CluIso), lengths(CluIso)))
CluIso$colour <- c(rep(NA,length(CluIso$Category[CluIso$Category != "DTE"])))


p4 <- ggTranPlots(gtf$targ_merged, class.files$targ_merged, isoList = c(as.character(CluIso$Isoform)),
            selfDf = CluIso, gene="CLU")


# plots
plot_grid(p1, plot_grid(p2), p3, p4, labels = c("A","B","C","D"))
plot_grid(p3, p4, labels = c("A","B"))

class.files$SCN[class.files$SCN$isoform == "PB.57985.21",]        

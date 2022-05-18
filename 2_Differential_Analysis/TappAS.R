# TAPPAS Results for ADBDR 

suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("wesanderson"))
suppressMessages(library("cowplot"))

library(extrafont)
#font_install('fontcm')
loadfonts()

mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=18,  family="CM Roman"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = c(.90, 0.95),
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6),
                 legend.text = element_text(size = 18,family="CM Roman"),
                 axis.text.x= element_text(size=16,family="CM Roman"),
                 axis.text.y= element_text(size=16,family="CM Roman"))

# plot label colour
label_colour <- function(genotype){
  if(genotype == "Braak0"){colour = wes_palette("Royal1")[1]}else{
    if(genotype == "Braak1"){colour = wes_palette("Royal2")[3]}else{
      if(genotype == "Braak6"){colour = wes_palette("Royal1")[2]}else{
        if(genotype == "novel"){colour = wes_palette("Darjeeling1")[4]}else{
          if(genotype == "known"){colour = wes_palette("Darjeeling1")[5]}else{
            if(genotype == "targeted"){colour = wes_palette("Darjeeling1")[2]}else{
              if(genotype == "whole"){colour = wes_palette("Darjeeling1")[1]}else{
                if(genotype == "whole+targeted"){colour = wes_palette("Darjeeling2")[1]}else{
                  if(genotype == "Control"){colour = wes_palette("Royal1")[2]}else{
                    if(genotype == "AD"){colour = wes_palette("Royal1")[1]}else{
                    }}}}}}}}}}
  return(colour)
}

output_plot_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/AD_BDR/Figures_Thesis"
input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq/TAPPAS/Results"

PostIsoSeq_root_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq"
targeted.class.files <- read.table(paste0(PostIsoSeq_root_dir, "/SQANTI_TAMA_FILTER/AllBDRTargeted_sqantitamafiltered.classification.txt"), sep = "\t", as.is = T, header = T)
phenotype <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/AD_BDR/Raw_Data/ADBDR_PhenotypeTAPPAS.txt", header = T) %>% mutate(col_names = paste0(group,".",sample))

# read in files generated from TAPPAS
files <- list.files(path = input_dir, pattern = ".tsv", full.names = T)
files <- lapply(files, function(x) read.table(x, sep = "\t", header = T))
names(files) <- list.files(path = input_dir, pattern = ".tsv")

# number of transcripts that are filtered for statistical purposes
files$tappAS_Transcripts_InputExpressionMatrix.tsv <- 
  merge(targeted.class.files[,c("isoform","associated_gene","associated_transcript","structural_category")], files$tappAS_Transcripts_InputExpressionMatrix.tsv, by.x = "isoform", by.y = "Id")

# tally of the number of transcripts filtered per gene 
files$tappAS_Transcripts_InputExpressionMatrix.tsv %>% group_by(associated_gene, structural_category, Filtered) %>% tally() %>% ggplot(., aes(x = associated_gene, y = n, fill = Filtered)) + geom_bar(stat = "identity") + labs(x = "Target Gene", y = "Number of Isoforms") + mytheme +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(labels = c("Retained","Removed due to low coverage"), values = c(wes_palette("Rushmore1")[3],wes_palette("Rushmore1")[5])) + theme(legend.position = c(0.8,0.8))

files$tappAS_Transcripts_InputExpressionMatrix.tsv %>% group_by(associated_gene, structural_category, Filtered) %>% tally() %>% filter(Filtered != "NO") %>% ggplot(., aes(x = associated_gene, y = n, fill = structural_category)) + geom_bar(stat = "identity") + labs(x = "Target Gene", y = "Number of Isoforms Removed") + mytheme +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position = c(0.8,0.8))

# Normalised Gene Expression Counts (Already filtered for low expression counts) 
Norm_transcounts <- merge(targeted.class.files[,c("isoform","associated_gene","associated_transcript","structural_category")], files$input_normalized_matrix.tsv, by.x = "isoform", by.y = 0) %>% reshape2::melt() %>% 
  left_join(., phenotype, by = c("variable" = "sample")) %>% 
  mutate(group = factor(group, levels = c("CONTROL", "CASE"),labels = c("CONTROL", "CASE")),
         structural_category=recode(structural_category, 
                                    `full-splice_match`="FSM",`novel_not_in_catalog`="NNC",`incomplete-splice_match`="ISM",
                                    `novel_in_catalog`="NIC"),
         Isoform = paste0(isoform,"_",structural_category))

GeneExp <- Norm_transcounts %>% group_by(associated_gene,variable) %>% summarise(Exp = sum(value)) %>%
  left_join(., phenotype, by = c("variable" = "sample"))

plot_mergedexp <- function(InputGene,IsoformID){
  if (InputGene != "NA"){
    df <- GeneExp %>% filter(associated_gene == InputGene) 
    df$group <- factor(df$group, levels = c("CONTROL","CASE")) %>% 
      recode(., `CONTROL`="Control",`CASE`="AD")
    plot_title <- InputGene
  }else if(IsoformID != "NA"){
    df <- merge(Norm_transcounts %>% filter(isoform == IsoformID), phenotype[,c("sample","col_names")], by.x = "variable", by.y = "sample")
    colnames(df)[6] <- "Exp"
    df$group <- factor(df$group, levels = c("CONTROL","CASE")) %>% 
      recode(., `CONTROL`="Control",`CASE`="AD")
    plot_title <- paste0(df$associated_gene,": ",df$associated_transcript)
  }else{
    print("2 variables required")
  }
  
  p <- ggplot(df, aes(x = group, y = Exp, colour = group)) + geom_boxplot() + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
    labs(y = "Normalised Gene Expression", x = "", title = paste0(plot_title,"\n\n")) + mytheme +
    # colours in the opposite direction of the levels 
    scale_colour_manual(values = c(label_colour("AD"),label_colour("Control"))) + 
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5, size = 16,face = "italic"))
  
  # subtitles
  if(IsoformID != "NA"){
    p <- p + labs(title = plot_title, subtitle = paste0(df$isoform,"\n\n"), y = "Normalised Isoform Expression") + theme(plot.subtitle = element_text(hjust = 0.5, size = 12,face = "italic"))
  }
  
  return(p)
}

plot_transexp <- function(InputGene){
  df <-  Norm_transcounts  %>% filter(associated_gene == InputGene) 
  
  p <- ggplot(df, aes(x = reorder(isoform,-value), y = value, colour = group)) + geom_boxplot() + 
    #stat_summary(data=df, aes(x=group, y=value, group=Isoform), fun ="mean", geom="line", linetype = "dotted") +
    mytheme + labs(x = "", y = "Normalised Isoform Expression",title = paste0(InputGene,"\n\n")) +
    scale_colour_manual(values = c(label_colour("AD"),label_colour("Control"))) +
    theme(strip.background = element_blank(), legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),
          panel.spacing = unit(2, "lines"), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(~structural_category,  scales = "free", space = "free") +
    scale_y_continuous(trans='log10')
  return(p)
}

IFtrans_plot <- function(InputGene){
  dat <- Norm_transcounts  %>% filter(associated_gene == InputGene)
  
  # group by the samples
  dat <- dat %>% group_by(variable) %>% mutate(IF = value/sum(value) * 100)
  # group by isoform, age and group 
  dat2 <- aggregate(.~Isoform+group, dat[,c("Isoform","group","IF")], mean) %>% mutate(IF = tidyr::replace_na(IF, 0))
  p <- ggplot(dat, aes(x = isoform, y = IF, fill = group)) + geom_boxplot() + 
    labs(x = "", y = "Isoform Fraction (%)", title = paste0(InputGene,"\n\n")) + mytheme + 
    guides(fill=guide_legend(ncol=3,bycol=TRUE)) + 
    theme(legend.position="none", strip.background = element_blank(),plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  scale_fill_manual(values = c(label_colour("AD"),label_colour("Control"))) + facet_grid(~structural_category,  scales = "free", space = "free") + ylim(0,100)
  
  p1 <- ggplot(dat2, aes(x = Isoform, y = IF, fill = group)) + geom_bar(stat = "identity", position = position_dodge()) +  
    labs(x = "Age (months)", y = "Isoform Fraction (%)", title = paste0(InputGene,"\n\n")) + mytheme + 
    guides(fill=guide_legend(ncol=3,bycol=TRUE)) + 
    theme(legend.position="bottom", strip.background = element_blank(),plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  scale_fill_manual(values = c(label_colour("AD"),label_colour("Control"))) 
  
  return(p)
}



# Gene Expression Plots 
geneexp_plots <- lapply(unique(GeneExp$associated_gene), function(gene) plot_mergedexp(gene,"NA"))
geneexp_plotsgrobs <- lapply(geneexp_plots , ggplotGrob)
names(geneexp_plotsgrobs) <- unique(GeneExp$associated_gene)

# Target Expression Plots
transexp_plots <- lapply(unique(GeneExp$associated_gene), function(gene) plot_transexp(gene) + theme(legend.position = "none"))
transexp_plotsgrobs <- lapply(transexp_plots, ggplotGrob)
names(transexp_plotsgrobs) <- unique(GeneExp$associated_gene)


pdf (paste0(output_plot_dir,"/DifferentialAnalysis.pdf"), width = 10, height = 15)
# Gene expression Analysis
plot_grid(geneexp_plotsgrobs$APP, geneexp_plotsgrobs$MAPT, geneexp_plotsgrobs$APOE, geneexp_plotsgrobs$BIN1, geneexp_plotsgrobs$SNCA, geneexp_plotsgrobs$TREM2, 
          labels = c("a","b","c","d","e","f"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
plot_grid(geneexp_plotsgrobs$ABCA7, geneexp_plotsgrobs$CD33, geneexp_plotsgrobs$CLU, geneexp_plotsgrobs$FUS, geneexp_plotsgrobs$FYN, geneexp_plotsgrobs$SORL1, labels = c("g","h","i","j","k","l"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
plot_grid(geneexp_plotsgrobs$TARDBP, geneexp_plotsgrobs$PTK2B, geneexp_plotsgrobs$PICALM,
          geneexp_plotsgrobs$RHBDF2,NULL,NULL,
          labels = c("m","n","o","p"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)

# Transcript 
plot_grid(transexp_plotsgrobs$APP, transexp_plotsgrobs$MAPT,ncol = 1, scale = 0.9)
plot_grid(transexp_plotsgrobs$APOE,transexp_plotsgrobs$BIN1,ncol = 1, scale = 0.9)
plot_grid(transexp_plotsgrobs$SNCA, transexp_plotsgrobs$TREM2, ncol = 1, scale = 0.9)
plot_grid(transexp_plotsgrobs$ABCA7, transexp_plotsgrobs$CD33, ncol = 1, scale = 0.9)
plot_grid(transexp_plotsgrobs$CLU,transexp_plotsgrobs$FUS, ncol = 1, scale = 0.9)
plot_grid(transexp_plotsgrobs$FYN, transexp_plotsgrobs$SORL1, ncol = 1, scale = 0.9)
plot_grid(transexp_plotsgrobs$TARDBP, transexp_plotsgrobs$PTK2B, ncol = 1, scale = 0.9)
plot_grid(transexp_plotsgrobs$PICALM,transexp_plotsgrobs$RHBDF2, ncol = 1, scale = 0.9)
dev.off()


plot_mergedexp("NA","PB.9744.148")
plot_mergedexp("NA","PB.5435.41")


# Isoform Fraction Plots
IF_plots <- lapply(unique(GeneExp$associated_gene), function(gene) IFtrans_plot(gene))
IF_plotsgrobs <- lapply(IF_plots, ggplotGrob)
pdf (paste0(output_plot_dir,"/DifferentialAnalysis_IF.pdf"), width = 10, height = 15)
for(i in seq(1,16,2)){print(plot_grid(IF_plotsgrobs[[i]],IF_plotsgrobs[[i+1]], ncol = 1, scale = 0.9))}
dev.off()


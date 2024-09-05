library(grid)  # for the unit() function

# chr9_104782004_T_C
abca1_phased <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/by_loci/PB.10258_size1136/phased.nopartial.cleaned.human_readable.txt", 
                           sep = "\t", as.is = T, header = T)

abca1_phased <- abca1_phased %>% reshape2::melt(id = "haplotype",
                                                variable.name = "isoform",
                                                value.name = "counts") %>% mutate(allele = ifelse(haplotype %in% c("TTGGG", "TTTAA"), "T", "C"))

gtf <- list(
  iso = rtracklayer::import(paste0(dirnames$targ_iso_root, "9_sqanti3/full/AllBDRTargeted.collapsed_classification.filtered_lite.gtf")),
  ont_merged = rtracklayer::import(paste0(dirnames$targ_ont_root, "5_cupcake/7_sqanti3/ontBDR_collapsed.filtered_counts_filtered.gtf")),
  iso_all = rtracklayer::import(paste0(dirnames$targ_iso_root, "9_sqanti3/full/AllBDRTargeted.collapsed_corrected.gtf"))
)

class.file <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/9_sqanti3/full/AllBDRTargeted.collapsed_classification.filtered_lite_classification.txt", header = T)
head(abca1_phased)

abca1_phased_tallied <- abca1_phased %>% group_by(isoform, allele) %>% tally(counts) %>% tidyr::spread(., key = allele, value = n)
plotOutput <- apply(abca1_phased_tallied, 1, function(x) plot_allele_perTranscript(gtf$iso_all, x[["isoform"]], x[["C"]],x[["T"]]))
names(plotOutput) <- abca1_phased_tallied$isoform

refs <- lapply(plotOutput, function(x) x$ref)
combined_refs <- do.call(rbind, refs)
combined_refs <- combined_refs %>% mutate(transcript = factor(transcript, levels = rev(abca1_phased_tallied$isoform)))

ref_lines <- lapply(plotOutput, function(x) x$ref_lines)
combined_refs_lines <- do.call(rbind, ref_lines)

alts <- lapply(plotOutput, function(x) x$alt)
combined_alts <- do.call(rbind, alts)

alt_lines <- lapply(plotOutput, function(x) x$alt_lines)
combined_alt_lines <- do.call(rbind, alt_lines)

p1 <- ggplot() +
  geom_segment(data = combined_refs, aes(x = start, xend = end, y = as.numeric(value), yend = as.numeric(value), colour = "Ref"), linetype = "solid") +
  geom_segment(data = combined_refs_lines, aes(x = position, xend = position, y = 0, yend = as.numeric(value), colour = "Ref"), linetype = "dotted") +
  geom_segment(data = combined_alt_lines, aes(x = position, xend = position, y = 0, yend = as.numeric(value), colour = "Alt"), linetype = "dotted") +
  geom_segment(data = combined_alts, aes(x = start, xend = end, y = as.numeric(value), yend = as.numeric(value), colour = "Alt"), linetype = "solid") +
  #facet_grid(transcript ~ .) +
  facet_grid(transcript ~ ., scale = "free") +
  scale_colour_manual(name = "Allele Type", 
                      values = c("Ref" = "blue", "Alt" = "orange"), 
                      labels = c("Ref" = "Reference", "Alt" = "Alternative")) +
  labs(x = NULL, y = "Number of FL reads") +
  geom_vline(xintercept=104782004, linetype="solid", color = "grey80") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5, linetype = "dotted"),
    legend.position = "bottom"# Add pale grey border
  ) 


output <- plot_allele_perTranscript(gtf$iso_all, "PB.10258.2", 22, 8)
ggplot() +
  geom_segment(data = output$ref, aes(x = start, xend = end, y = value, yend = value), linetype = "solid", colour = "blue") +
  geom_segment(data = output$ref_lines, aes(x = position, xend = position, y = 0, yend = value), linetype = "dotted", colour = "blue") +
  geom_segment(data = output$alt_lines, aes(x = position, xend = position, y = 0, yend = value), linetype = "dotted", colour = "orange") +
  geom_segment(data = output$alt, aes(x = start, xend = end, y = value, yend = value), linetype = "solid", colour = "orange") +
  labs(x = "Position", y = "Value") +
  theme_classic() +
  theme(axis.text.y = element_text(angle=90)) 

datGtf <- gtf$iso_all %>%
  filter(transcript_id %in% abca1_phased_tallied$isoform & type == "exon") %>%
  mutate(transcript_id = factor(transcript_id, levels = rev(abca1_phased_tallied$isoform)))
p2 <- datGtf %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range() + labs(x = NULL, y = "") +
  geom_intron(data = to_intron(datGtf, "transcript_id"),
              aes(strand = strand)) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text on the left
    axis.ticks.y = element_blank(),  # Remove y-axis ticks on the left
    axis.line.y = element_blank(),  # Remove y-axis line on the left
    axis.text.y.left = element_text(color = "black")  # Add y-axis text on the right
  ) + theme(axis.text.y = element_text(angle = 90),
            plot.margin = unit(c(1, 1, 1, 4), "mm"))
plot_grid(p1,p2, nrow=2)


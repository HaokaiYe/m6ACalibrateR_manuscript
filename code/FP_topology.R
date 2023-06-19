

########################################################################################################################
####        Randomly Capturing High-coverage Consensus Sequences to Reconstruct False Positive m6A Landscapes       ####
########################################################################################################################

################################################################################
####                         Quantifying mapped reads                       ####
################################################################################

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(predictiveFeatures)
library(GenomicAlignments)
library(SummarizedExperiment)
library(BiocParallel)

## Extract motif on hg38
motif_hg38 <- sort(predictiveFeatures::sampleSequence(
  "DRACH",
  exons(TxDb.Hsapiens.UCSC.hg38.knownGene),
  BSgenome.Hsapiens.UCSC.hg38
))

## Prepare sample information
bam_path = "~/merip/bam_human_SYSY"
bam_names = c("SRR14765584.bam", "SRR14765585.bam", "SRR14765586.bam", "SRR14765587.bam", 
              "SRR14765588.bam", "SRR14765589.bam", "SRR14765590.bam", "SRR14765591.bam")
IVT_mRNA = c("IVT","IVT","IVT","IVT","mRNA","mRNA","mRNA","mRNA")
IP_input = c("input","IP","input","IP","input","IP","input","IP")

## Read count
featuresCount <- function(features, bam_dirs, parallel = 1, yield_size = 5e+06) {
  register(SerialParam())
  suppressWarnings(register(MulticoreParam(workers = parallel)))
  register(SnowParam(workers = parallel))
  bam_lst = BamFileList(file = bam_dirs, asMates = TRUE)
  yieldSize(bam_lst) = yield_size
  se <- summarizeOverlaps(features = features, 
                          reads = bam_lst, 
                          mode = "Union", 
                          inter.feature = FALSE, 
                          singleEnd = FALSE, 
                          ignore.strand = TRUE, 
                          fragments = TRUE)
  return(se)
}
se <- featuresCount(motif_hg38, file.path(bam_path, bam_names))
colData(se) <- DataFrame(bam_names = bam_names,
                         IVT_mRNA = IVT_mRNA,
                         IP_input = IP_input)
saveRDS(se, "./rds/se_human_SYSY.rds")


################################################################################
####                     Extract all DRACH motif on exons                   ####
################################################################################

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

exbtx <- exonsBy(txdb, by = "tx")
motif_all <- sort(sampleSequence("DRACH", exbtx, bsgenome) - 2)

################################################################################
####            Using IVT and mRNA data to determine two groups:            ####
####                    true positive and false positive                    ####
################################################################################

peaks_mRNA <- readRDS("./rds/peaks_mRNA.rds")
peaks_IVT <- readRDS("./rds/peaks_IVT.rds")

# Modification sites identified exclusively in the mRNA sample were considered true positives,
# while all sites identified in the IVT sample were deemed false positives
TP = peaks_mRNA[-queryHits(findOverlaps(peaks_mRNA, peaks_IVT))]
FP = peaks_IVT

################################################################################
####                   Establishing high-coverage motifs                    ####
################################################################################

se_human_SYSY <- readRDS("./rds/se_human_SYSY.rds")
motif_all$count = rowMeans(assay(se_human_SYSY)[,se_human_SYSY$IP_input == "input"])

# get true m6a input counts
true_m6A = subsetByOverlaps(motif_all, TP)
non_m6A = subsetByOverlaps(motif_all, TP, invert = T)

threshold = mean(true_m6A$count)
high_cov = non_m6A[non_m6A$count >= threshold]


################################################################################
####                                 metagene                               ####
################################################################################

library(ggplot2)
library(patchwork)

relative_pos_on_region <- function(x, region){
  region_map <- mapToTranscripts(x, region)
  region_width <- sum(width(region))[region_map$transcriptsHits]
  start_on_region <- start(region_map)
  return(start_on_region/region_width)
}
metagene <- function(marker, txdb, region_weights, marker_name){
  u5bytx <- fiveUTRsByTranscript(txdb)
  cdsbytx <- cdsBy(txdb, by = "tx")
  u3bytx <- threeUTRsByTranscript(txdb)
  
  utr5_pos <- relative_pos_on_region(marker, u5bytx) * region_weights[1]
  cds_pos <- (relative_pos_on_region(marker, cdsbytx) * region_weights[2]) + region_weights[1]
  utr3_pos <- (relative_pos_on_region(marker, u3bytx) * region_weights[3]) + region_weights[1] + region_weights[2]
  
  
  pldf <- data.frame(relative_pos = c(utr5_pos, cds_pos, utr3_pos),
                     tx_region = rep(c("5'UTR","CDS","3'UTR"),
                                     c(length(utr5_pos),length(cds_pos),length(utr3_pos))))
  
  pldf$tx_region = factor(pldf$tx_region,levels = c("5'UTR","CDS","3'UTR"))
  pldf$class <- rep(marker_name, nrow(pldf))
  return(pldf)
}

# calculate region_weights
u5bytx <- fiveUTRsByTranscript(txdb)
cdsbytx <- cdsBy(txdb, by = "tx")
u3bytx <- threeUTRsByTranscript(txdb)

UTR5_width = mean(sum(width(u5bytx)))
CDS_width = mean(sum(width(cdsbytx)))
UTR3_width = mean(sum(width(u3bytx)))
sum_width = UTR5_width + CDS_width + UTR3_width

region_weights = c(UTR5_width/sum_width, CDS_width/sum_width, UTR3_width/sum_width) # 0.09078195 0.45142217 0.45779587
region_weights = c(0.09078195, 0.45142217, 0.45779587)

# calculate metagene
topology_TP = metagene(TP, txdb, region_weights, "TP")
topology_FP = metagene(FP, txdb, region_weights, "FP")
topology_motif = metagene(motif_all, txdb, region_weights, "All DRACH")
topology_high_cov = metagene(high_cov, txdb, region_weights, "High Cov. DRACH")

metagene_df = rbind(topology_TP, topology_FP, topology_motif, topology_high_cov)
metagene_df$class <- factor(metagene_df$class, levels = c("TP", "FP", "All DRACH", "High Cov. DRACH"))

metagene_plot = ggplot() +
  geom_density(data = metagene_df[which(metagene_df$class == "All DRACH"),], aes(x = relative_pos), fill = "#f2f2f2", color = NA) +
  geom_line(data = metagene_df, aes(x = relative_pos, y = ..density.., color = class), stat = "density") +
  geom_vline(xintercept = c(region_weights[1], region_weights[1] + region_weights[2]), linetype = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  ylab("Density") +
  scale_color_manual(values = c("#B1182D",
                                "#2464AB",
                                "#91C5DA",
                                "#707070"),
                     breaks = c("TP", "FP", "High Cov. DRACH", "All DRACH"),
                     labels = c(expression(paste("High-conf. m"^6, "A")), "False positives", "High-cov. DRACH", "All DRACH")) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour= "black", size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour= "black", size = 12),
        legend.position = c(0.8, 0.8),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(colour = "black", size = 10))


tx_plot = ggplot(data.frame(x = c(0, 1), y = c(0, 1)), aes(x = x, y = y)) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0), linetype = "solid", size = 1.5) +
  geom_rect(aes(xmin = region_weights[1], xmax = region_weights[1] + region_weights[2], ymin = -0.05, ymax = 0.05), fill = "#7e7e7e", color = "black", size = 0.5) +
  annotate("text", x = region_weights[1]/2, y = -0.15, label = "5'UTR", size = 3.6) +
  annotate("text", x = region_weights[2]/2 + region_weights[1], y = -0.15, label = "CDS", size = 3.6) +
  annotate("text", x = region_weights[3]/2 + region_weights[1] + region_weights[2], y = -0.15, label = "3'UTR", size = 3.6) +
  ylim(-0.25, 0.05) +
  theme_void()

metagene_plot/tx_plot + plot_layout(heights = c(10, 1))

ggsave("~/plots/metagene.pdf", width = 5.0, height = 4.0)




########################################################################################################################
####                                               Analyze top 2 features                                           ####
########################################################################################################################

library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene


################################################################################
####                     Load mRNA peaks and IVT peaks                      ####
################################################################################

extractRegionLength <- function(x,
                                region = NULL,
                                ambiguityMethod = c("mean", "sum", "min", "max"),
                                maxgap = -1L,
                                minoverlap = 0L,
                                type = c("any", "start", "end", "within", "equal"),
                                nomapValue = c("NA", "0", "nearest"),
                                ignore.strand = FALSE) {
  stopifnot(is(x, "GRanges"))
  ambiguityMethod <- match.arg(ambiguityMethod)
  nomapValue <- match.arg(nomapValue)
  type <- match.arg(type)
  
  if (is.null(region)) {
    length_property <- width(x)
  } else if (is(region, "GRanges")) {
    length_property <- width(region)
  } else if (is(region, "GRangesList")) {
    region <- region[elementNROWS(region) != 0]
    length_property <- sum(width(region))
  } else{
    stop("`region` should be either `GRanges` or `GRangesList`")
  }
  
  if (is.null(region)) {
    return(length_property)
  } else{
    region_property <- extractRegionProperty(
      x = x,
      property = length_property,
      region = region,
      ambiguityMethod = ambiguityMethod,
      maxgap = maxgap,
      minoverlap = minoverlap,
      type = type,
      nomapValue = nomapValue,
      ignore.strand = ignore.strand
    )
    return(region_property)
  }
}

plot_2D_density_equal_data <- function(peaks_IVT, peaks_mRNA, txdb, antibody) {
  
  TP = peaks_mRNA[-queryHits(findOverlaps(peaks_mRNA, peaks_IVT))]
  FP = peaks_IVT
  
  # equal size
  if(length(TP) > length(FP)) {
    set.seed(110)
    ind_TP = sample(1:length(TP), length(FP))
    TP = TP[ind_TP,]
  } else {
    set.seed(110)
    ind_FP = sample(1:length(FP), length(TP))
    FP = FP[ind_FP,]
  }
  
  AP_grg <- rbind(TP, FP)
  
  exons <- exons(txdb)
  log2_length_exons <- log2(extractRegionLength(AP_grg, exons) + 1)
  
  exonicTranscripts <- exonsBy(txdb, by = "tx")
  log2_length_exonicTranscripts <- log2(extractRegionLength(AP_grg, exonicTranscripts) + 1)
  
  AP_gfeatures <- as.data.frame(cbind(log2_length_exons, log2_length_exonicTranscripts))
  AP_gfeatures$class <- factor(c(rep(1, length(TP)), rep(0, length(FP))), levels = c("1", "0"), labels = c("TP", "FP"))
  AP_gfeatures$antibody <- antibody
  
  return(AP_gfeatures)
}

peaks_mRNA <- readRDS("./rds/peaks_mRNA.rds")
peaks_IVT <- readRDS("./rds/peaks_IVT.rds")

SYSY_data <- plot_2D_density_equal_data(peaks_IVT, peaks_mRNA, txdb_hg38, "SYSY")

################################################################################
####                         combined density plot                          ####
################################################################################

top2features = readRDS("./rds/top2features.rds")

density = ggplot(top2features, aes(x = log2_length_exons, y = log2_length_exonicTranscripts, color = class)) +
  geom_density_2d(h=0.75) +
  ylab(expression(paste("log"[2], "(mRNA length)"))) +
  xlab(expression(paste("log"[2], "(exon length)"))) +
  labs(color = "Class") +
  annotate('text', x = 7, y = 13.15, label = "False\npositives", size = 3.6, color = '#0076ba', hjust = 0.5, lineheight = 0.8) +
  annotate('text', x = 12, y = 10.2, label = expression(paste("High-conf. m"^6, "A")), size = 3.6, color = '#f27200', lineheight = 0.8) +
  scale_color_manual(values = c("#0076ba", "#f27200")) +
  theme_classic() + 
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12),
        legend.position = "none",
        plot.margin = margin(15, 5, 0, 5),
        strip.background = element_rect(fill = "#d9d9d9", size = 0.8),
        strip.text = element_text(colour = "black", size = 10)) +
  facet_wrap(facet = antibody ~ ., scales = "free")


################################################################################
####                                  LR fit                                ####
################################################################################

top2features$prob <- ifelse(top2features$class == "TP", 1, ifelse(top2features$class == "FP", 0, NA))
top2features$antibody <- factor(top2features$antibody, levels = c("Abcam", "NEB", "SYSY", "Mouse"), labels = c("Abcam", "NEB", "SYSY", "Mouse"))

exon = ggplot(top2features, aes(x = log2_length_exons, y = prob, color = antibody)) +
  geom_smooth(formula = y ~ splines::ns(x, 5), 
              method = "glm", method.args = list(family = "binomial")) + 
  geom_hline(yintercept = 0.5, linetype = 3) + 
  ylim(0,1) +
  xlim(6,13) +
  labs(x = expression(paste("log"[2], "(exon length)")), y = expression(paste("Prob of true m"^6, "A"))) + 
  scale_color_manual(values = c("dodgerblue4",
                                "darkorchid3",
                                "#d75427",
                                "darkolivegreen4")) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  theme_classic() + 
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour= "black", size = 12, vjust = 2),
        legend.title = element_blank(),
        legend.text = element_text(colour= "black", size = 10),
        legend.position = c(0.18, 0.7),
        legend.background = element_blank()) +
  annotate("segment", x = 11, xend = 13, y = 0.25, yend = 0.25, colour = "black",
           arrow = arrow(length = unit(0.2,"cm")), linetype = 1, size = 0.6) +
  annotate('text', x = 12, y = 0.13, label = expression(paste("log"[2], "(exon length)")), size = 3.6, color = "black")


mRNA = ggplot(top2features, aes(x = log2_length_exonicTranscripts, y = prob, color = antibody)) +
  geom_smooth(formula = y ~ splines::ns(x, knots = seq(10, 13, by = 0.5)), 
              method = "glm", method.args = list(family = "binomial")) + 
  geom_hline(yintercept = 0.5, linetype = 3) + 
  ylim(0,1) +
  xlim(9,13.6) +
  labs(x = expression(paste("log"[2], "(exon length)")), y = expression(paste("Prob of true m"^6, "A"))) + 
  scale_color_manual(values = c("dodgerblue4",
                                "darkorchid3",
                                "#d75427",
                                "darkolivegreen4")) +
  theme_classic() + 
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none') +
  annotate("segment", x = 12, xend = 13.6, y = 0.25, yend = 0.25, colour = "black",
           arrow = arrow(length = unit(0.2,"cm")), linetype = 1, size = 0.6) +
  annotate('text', x = 12.85, y = 0.13, label = expression(paste("log"[2], "(mRNA length)")), size = 3.6, color = "black")

(exon | mRNA)/density + plot_layout(heights = c(1, 4))

ggsave("./plots/top2features.pdf", width = 7.3, height = 7.6)


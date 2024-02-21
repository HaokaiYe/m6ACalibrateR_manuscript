
########################################################################################################################
####                                      available for peak calling methods                                        ####
########################################################################################################################

library(scales)
library(ggplot2)
library(m6ACalibrateR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

motif_all <- readRDS("~/rds/motif_all.rds") - 2
pred_df <- readRDS("~/rds/pred_df.rds")
TP = motif_all[which(pred_df$mean_p0 <= 0.5),]

##############################################################
####                     ground truth                     ####
##############################################################

single <- readRDS("~/rds/m6A_Atlas2_high_resolution_hg38.rds")
groundTrt = single[single$Technique_Num >= 2]
free = groundTrt

################################################################################
####                                Curve                                   ####
################################################################################
##############################################################
####                    GSE112795;HeLa;EMT                ####
##############################################################

calculate_cumulative_overlap <- function(peaks, free) {
  overlap_counts_all <- countOverlaps(peaks, free)
  
  cumulative_overlap_counts <- cumsum(overlap_counts_all)
  
  cumulative_peak_length <- cumsum(width(peaks))
  
  df <- data.frame(cumulative_peak_length, cumulative_overlap_counts)
  return(df)
}

peaks <- readRDS("./rds/GSE112795_HeLa_EMT_exomePeak2.rds")
peaks_aft = m6ACalibrate(peaks, txdb = txdb, genome = bsgenome, FP_threshold = 0.5)


befo <- calculate_cumulative_overlap(peaks, free)
aft <- calculate_cumulative_overlap(peaks_aft, free)
df_exomePeak2 = rbind(befo, aft)
df_exomePeak2$class = c(rep("Before", nrow(befo)), rep("After", nrow(aft)))


##############################################################
####                GSE112795;HeLa;EMT (MACS2)            ####
##############################################################

peaks <- readRDS("./rds/GSE112795_HeLa_EMT_MACS2.rds")
peaks_aft = m6ACalibrate(peaks, txdb = txdb, genome = bsgenome, FP_threshold = 0.5)


befo <- calculate_cumulative_overlap(peaks, free)
aft <- calculate_cumulative_overlap(peaks_aft, free)
df_MACS2 = rbind(befo, aft)
df_MACS2$class = c(rep("Before", nrow(befo)), rep("After", nrow(aft)))


df_exomePeak2$method = "exomePeak2"
df_MACS2$method = "MACS2"
df_combined = rbind(df_exomePeak2, df_MACS2)

ggplot(df_combined, aes(x=cumulative_peak_length/1000000, y=cumulative_overlap_counts, color = class)) +
  geom_line() +
  xlab("Cumulative width of top peaks (million bases)") +
  ylab("m6A-Atlas2 sites (#)") +
  scale_x_continuous(limits = c(0, 2), labels = comma) +
  scale_y_continuous(limits = c(250, 20000), labels = comma, trans = log10_trans()) +
  scale_color_manual(values = c("#0076ba", "#cb3335"),
                     name = NULL,
                     breaks = c("Before", "After"),
                     labels = c("Pre-Calibration", "Post-Calibration")) +
  theme_classic() + 
  theme(axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 10),
        axis.text.y = element_text(colour = "black", size = 10, angle = 90, hjust = 0.5),
        axis.title = element_text(colour = "black", size = 12),
        axis.title.y = element_text(colour = "black", size = 12, vjust = 2),
        legend.title = element_blank(),
        legend.position = c(0.87, 0.2),
        legend.background = element_blank(),
        legend.text = element_text(colour = "black", size = 10),
        strip.text = element_text(colour = "black", size = 11, hjust = 0),
        strip.background = element_blank(),
        plot.margin = margin(5, 10, 5, 7),
        panel.spacing = unit(0.5, "cm")) +
  facet_wrap(facet = . ~ method, scales = "free")

ggsave("./plots/callmethods.pdf", width = 5.8, height = 2.9)


################################################################################
####                          Area under the curve                          ####
################################################################################

# top width
calculate_cumulative_overlap <- function(peaks, free) {
  # 使用countOverlaps直接计算每个peak的重叠数量
  overlap_counts_all <- countOverlaps(peaks, free)
  
  # 计算累积的重叠数量
  cumulative_overlap_counts <- cumsum(overlap_counts_all)
  
  # 计算累积的peak长度
  cumulative_peak_length <- cumsum(width(peaks))
  
  df <- data.frame(cumulative_peak_length, cumulative_overlap_counts)
  return(df)
}

library(pracma)
getSciPval <- function(p_value) {
  scientific_p <- sprintf("%.1e", p_value) 
  split_p <- strsplit(scientific_p, "e")[[1]]
  p_main <- split_p[1]
  p_exp <- sub("^\\+0*", "", split_p[2])  # Remove leading zeros and plus
  p_exp <- sub("^\\-", "-", p_exp)  # Ensure the minus is preserved
  formatted_p <- paste0(p_main, " %*% 10^", p_exp)
  return(formatted_p)
}
cumWidthAUC = function(df, cumWidth) {
  indx = which(df$cumulative_peak_length <= cumWidth)
  AUC <- trapz(df$cumulative_peak_length[indx], df$cumulative_overlap_counts[indx])
  return(AUC)
}
scale_to_range <- function(df, min_val, max_val) {
  rng <- range(df, na.rm = TRUE)
  scaled_df <- (df - rng[1]) / (rng[2] - rng[1]) * (max_val - min_val) + min_val
  return(scaled_df)
}

human <- readRDS("~/merip/human.rds")
peakData = human[which(human$CallPeak_Methods == "MACS2")]
peakData = human[which(human$CallPeak_Methods == "exomePeak2")]
conditions = unique(peakData$condition)

AUCDf = data.frame()
filterSample = c()
cumWidth = 2000000
for (i in 1:length(conditions)) {
  peaks = peakData[which(peakData$condition == conditions[i])]
  peaks = peaks[order(peaks$pvalue, decreasing = T)]
  peaks_aft = subsetByOverlaps(peaks, TP)
  
  befo <- calculate_cumulative_overlap(peaks, free)
  aft <- calculate_cumulative_overlap(peaks_aft, free)
  
  if(max(aft$cumulative_peak_length) <= cumWidth) {
    filterSample = c(filterSample, i)
    next
  }
  
  tmp = data.frame("Before" = cumWidthAUC(befo, cumWidth),
                   "After" = cumWidthAUC(aft, cumWidth))
  AUCDf = rbind(AUCDf, tmp)
}

AUCDfScaled = scale_to_range(AUCDf, 0.02, 0.98)


########################################
####            dot box plot        ####
########################################

# dot plot
dot_plot = ggplot(data = AUCDfScaled, aes(x = Before, y = After)) +
  geom_point(color = "#b31a2f", size = 2.5, alpha = 0.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#477ca8", size = 0.8) +
  xlab("Pre-Calibration AUC") +
  ylab("Post-Calibration AUC") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
  theme_classic() + 
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12))


# box plot
pDCBoxDf <- reshape2::melt(AUCDfScaled, id.vars = NULL, variable.name = "corType", value.name = "corValue")

wilcoxTest = wilcox.test(AUCDfScaled$After, AUCDfScaled$Before, alternative = "greater")
maxValue = round(boxplot(pDCBoxDf$corValue[which(pDCBoxDf$corType == "After")])$stats[5, ], 2)
maxValue = maxValue + 0.02
shiftVal = 0.02
minValue = 0.01
boxplot_inset = ggplot(pDCBoxDf, aes(x = corType, y = corValue, fill = corType)) + 
  geom_boxplot(outlier.colour = NA, width = 0.55) +
  scale_fill_manual(values = c("#477ca8", "#cb3335")) +
  stat_boxplot(geom = "errorbar", aes_string(ymin = "..ymax.."), width = 0.25, size = .3) +
  stat_boxplot(geom = "errorbar", aes_string(ymax = "..ymin.."), width = 0.25, size = .3) +
  labs(y = "AUC") +
  annotate("segment", x = 1, xend = 1, y = maxValue + shiftVal, yend = maxValue + 2*shiftVal) +
  annotate("segment", x = 2, xend = 2, y = maxValue + shiftVal, yend = maxValue + 2*shiftVal) +
  annotate("segment", x = 1, xend = 2, y = maxValue + 2*shiftVal, yend = maxValue + 2*shiftVal) +
  annotate("text", x = 1.5, y = 1.14, label = getSciPval(wilcoxTest$p.value), parse = TRUE, size = 3.3) +
  coord_cartesian(ylim = c(minValue, maxValue - shiftVal), clip = "off") +
  scale_x_discrete(labels = c("Pre", "Post")) + 
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 9),
        axis.title = element_text(colour = "black", size = 11),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.background = element_blank(),
        plot.margin = margin(10, 5, 5, 5))

dot_plot + 
  inset_element(p = boxplot_inset, left = 0.48, bottom = 0, right = 1.02, top = 0.58)
ggsave("~/plots/AUCDotBox.pdf", width = 3.5, height = 3.3)



################################################################################
####                                Precision                               ####
################################################################################

##############################################################
####                    GSE112795;HeLa;EMT                ####
##############################################################

getPrecision <- function(peaks, cutoffs) {
  precision = round(length(peaks[peaks$class == 1]) / length(peaks), 4)
  precision_df = data.frame(cutoff = 1,
                            precision = precision,
                            peakNum = length(peaks))
  
  for (cutoff in cutoffs) {
    TP = motif_all[which(pred_df$mean_p0 <= cutoff),]
    peaks_aft = subsetByOverlaps(peaks, TP)
    
    precision = round(length(peaks_aft[peaks_aft$class == 1]) / length(peaks_aft), 4)
    
    temp_df = data.frame(cutoff = cutoff,
                         precision = precision,
                         peakNum = length(peaks_aft))
    
    precision_df = rbind(precision_df, temp_df)
  }
  
  return(precision_df)
}
plotPrecision = function(peaks, title) {
  peaks$class = 0
  overlap = findOverlaps(peaks, groundTrt)
  peaks$class[queryHits(overlap)] = 1
  
  precision_df = getPrecision(peaks, cutoffs)
  precision_df$cutoff = factor(precision_df$cutoff, levels = c(1, rev(cutoffs)))
  precision_df$class = title
  return(precision_df)
}

cutoffs <- c(0.1, 0.3, 0.5)
peaks_exo = readRDS("./rds/GSE112795_HeLa_EMT_exomePeak2.rds")
exo_df = plotPrecision(peaks_exo, "exomePeak2")
peaks_macs = readRDS("./rds/GSE112795_HeLa_EMT_MACS2.rds")
macs_df = plotPrecision(peaks_macs, "MACS2")

precision_df = rbind(exo_df, macs_df)

ggplot(precision_df, aes(x = cutoff, y = precision)) +
  geom_bar(stat = "identity", fill = "#f2756d", width = 0.8) +
  geom_bar(data = precision_df[which(precision_df$cutoff == 1),], stat = "identity", fill = "#1ebdc1", width = 0.8) +
  geom_text(aes(label = sprintf("%.2f%%", precision * 100)), vjust = -0.5, size = 3.2) +
  scale_y_continuous(labels = scales::percent_format(suffix = "")) +
  scale_x_discrete(breaks = c(1, rev(cutoffs)),
                   labels = c("Baseline", "Flexible", "Balanced", "Rigorous")) +
  labs(y = "m6A-Atlas2 sites (%)", x = "Calibration Threshold") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.text.x = element_text(angle = 30, vjust = 0.7, hjust = 0.5),
        axis.title = element_text(colour = "black", size = 12),
        plot.title = element_text(colour = "black", size = 12),
        strip.text = element_text(colour = "black", size = 11, hjust = 0),
        strip.background = element_blank()) +
  facet_wrap(facet = . ~ class, scales = "free")

ggsave("~/plots/precisionBar.pdf", width = 5.8, height = 4.1)


##############################################################
####                       All sample                     ####
##############################################################

human <- readRDS("~/merip/human.rds")
peakData = human[which(human$CallPeak_Methods == "MACS2")]
peakData = human[which(human$CallPeak_Methods == "exomePeak2")]
conditions = unique(peakData$condition)

precisionDf = data.frame()
for (i in 1:length(conditions)) {
  peaks = peakData[which(peakData$condition == conditions[i])]
  peaks$class = 0
  overlap = findOverlaps(peaks, groundTrt)
  peaks$class[queryHits(overlap)] = 1
  
  peaks_aft = subsetByOverlaps(peaks, TP)
  precision = round(length(peaks[peaks$class == 1]) / length(peaks), 4)
  precision_aft = round(length(peaks_aft[peaks_aft$class == 1]) / length(peaks_aft), 4)
  
  tmp = data.frame("Before" = precision, "After" = precision_aft, "lengthPre" = length(peaks), "lengthPost" = length(peaks_aft))
  precisionDf = rbind(precisionDf, tmp)
}


library(ggplot2)
library(patchwork)
library(ggpointdensity)
dot = ggplot(data = precisionDf, aes(x = Before, y = After)) + 
  #geom_hex(bins = 50) +
  geom_pointdensity(size = 2.2, alpha = 0.7, show.legend = T) +
  scale_color_gradientn(colours = c("#4a95c5", "#2462ab", "#b31a2f", "#b31a2f")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.8) +
  xlab("Pre-Calibration m6A-Atlas2 Sites (%)") +
  ylab("Post-Calibration m6A-Atlas2 Sites (%)") +
  scale_y_continuous(labels = percent_format(suffix = ""), breaks = seq(0, 1, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(labels = percent_format(suffix = ""), breaks = seq(0, 1, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  annotate('text', x = 0.82, y = 0.59, label = "Sample Pairs", size = 4.0, color = 'black', hjust = 0.5) +
  theme_classic() + 
  theme(axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 11),
        legend.position = c(0.86, 0.3),
        #legend.title = element_text(colour= "black", size = 12),
        legend.title = element_blank(),
        legend.text = element_text(colour= "black", size = 10),
        legend.background = element_blank(),
        plot.margin = margin(5, 12, 5, 5))

# histogram + density

density_befo = ggplot(data = precisionDf, aes(x = Before)) +
  geom_histogram(bins = 40, alpha = 0.7, color = "white", fill = "#96c4db") +
  geom_density(aes(y = ..density.. * 9.2), alpha = 0.1, color = "#2164ab") +
  scale_y_continuous(breaks = c(0, 30), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  theme_classic() + 
  theme(axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour= "black", size = 10),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank())

density_aft = ggplot(data = precisionDf, aes(x = After)) +
  geom_histogram(bins = 40, alpha = 0.7, color = "white", fill = "#d66255") +
  geom_density(aes(y = ..density.. * 9.2), alpha = 0.1, color = "#62061d") +
  scale_y_continuous(breaks = c(0, 30), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  theme_classic() + 
  theme(axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour= "black", size = 10),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank()) +
  coord_flip()

layout <- c(
  area(t = 0, l = 0, b = 15, r = 50),
  area(t = 16, l = 0, b = 100, r = 50),
  area(t = 16, l = 51, b = 100, r = 60)
)
density_befo + dot + density_aft + plot_layout(design = layout)

ggsave("~/plots/precisionDot.pdf", width = 4.25, height = 3.80)

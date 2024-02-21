
################################################################################
####                                Curve                                   ####
################################################################################

library(scales)
library(ggplot2)

motif_all <- readRDS("~/rds/motif_all.rds") - 2
pred_df <- readRDS("~/rds/pred_df.rds")
TP = motif_all[which(pred_df$mean_p0 <= 0.5),]

##############################################################
####                     ground truth                     ####
##############################################################

single <- readRDS("~/rds/m6A_Atlas2_high_resolution_hg38.rds")
groundTrt = single[single$Technique_Num >= 2]
free = groundTrt

##############################################################
####                      plot curve                      ####
##############################################################

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

human <- readRDS("~/merip/human.rds")
peakData = human[which(human$CallPeak_Methods == "MACS2")]
conditions = unique(peakData$condition)
TP = motif_all[which(pred_df$mean_p0 <= 0.5),]


plotTopWidth = function(i) {
  sample_info <- conditions[i]
  info_parts <- strsplit(sample_info, ";")[[1]]
  formatted_sample_info <- paste(info_parts[2], " (", info_parts[3], ")", sep = "")
  
  peaks = peakData[which(peakData$condition == conditions[i])]
  peaks = peaks[order(peaks$pvalue, decreasing = T)]
  peaks_aft = subsetByOverlaps(peaks, TP)
  
  befo <- calculate_cumulative_overlap(peaks, free)
  aft <- calculate_cumulative_overlap(peaks_aft, free)
  df_exomePeak2 = rbind(befo, aft)
  df_exomePeak2$class = c(rep("Before", nrow(befo)), rep("After", nrow(aft)))
  
  df_exomePeak2 = df_exomePeak2[which(df_exomePeak2$cumulative_peak_length < max(aft$cumulative_peak_length)),]
  
  ggplot(df_exomePeak2, aes(x=cumulative_peak_length/1000000, y=cumulative_overlap_counts, color = class)) +
    geom_line() +
    xlab("Cumulative width of top peaks (million bases)") +
    ylab("m6A-Atlas2 sites (#)") +
    ggtitle(formatted_sample_info) +
    scale_color_manual(values = c("#0076ba", "#cb3335"),
                       name = NULL,
                       breaks = c("Before", "After"),
                       labels = c("Pre-Calibration", "Post-Calibration")) +
    theme_classic() + 
    theme(axis.ticks = element_line(color = 'black'),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 11),
          legend.title = element_blank(),
          legend.position = c(0.75, 0.2),
          legend.background = element_blank(),
          legend.text = element_text(colour = "black", size = 10),
          plot.title = element_text(colour = "black", size = 12))
}

plotTopWidth(1) + scale_x_continuous(limits = c(0, 2)) + scale_y_continuous(limits = c(100, 18000), labels = comma, trans = log10_trans())
ggsave("~/plots/AUCCurve1.pdf", width = 3.8, height = 3.3)

plotTopWidth(2) + scale_x_continuous(limits = c(0, 2)) + scale_y_continuous(limits = c(100, 15000), labels = comma, trans = log10_trans())
ggsave("~/plots/AUCCurve2.pdf", width = 3.8, height = 3.3)

plotTopWidth(3) + scale_x_continuous(limits = c(0, 2)) + scale_y_continuous(limits = c(100, 17000), labels = comma, trans = log10_trans())
ggsave("~/plots/AUCCurve3.pdf", width = 3.8, height = 3.3)

plotTopWidth(4) + scale_x_continuous(limits = c(0, 2)) + scale_y_continuous(limits = c(100, 16000), labels = comma, trans = log10_trans())
ggsave("~/plots/AUCCurve4.pdf", width = 3.8, height = 3.3)


################################################################################
####                          Area under the curve                          ####
################################################################################

library(pracma)
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


###################################################
####               plot all samples            ####
###################################################

########################################
####            dot box plot        ####
########################################

# dot plot
dot_plot = ggplot(data = AUCDfScaled, aes(x = Before, y = After)) +
  geom_point(color = "#b31a2f", size = 2.5, alpha = 0.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#477ca8", size = 0.8) +
  xlab("Pre-calibration AUC") +
  ylab("Post-calibration AUC") +
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
ggsave("~/plots/AUCDotBox.pdf", width = 3.3, height = 3.2)



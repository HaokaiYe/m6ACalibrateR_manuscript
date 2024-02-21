
########################################################################################################################
####                    negative correlation between m6A level and gene expression level                            ####
########################################################################################################################

library(ggplot2)
motif_all <- readRDS("~/rds/motif_all.rds") - 2
pred_df <- readRDS("~/rds/pred_df.rds")
TP = motif_all[which(pred_df$mean_p0 <= 0.5),]

################################################################################
####                           neg cor in normPeak                          ####
################################################################################

# atlas ctrl samples (new motif-based)
merip_hg38_peaks = list.files("~/merip/siteRDS", full.names = T)
merip_hg38_peaks = grep("Ctrl", merip_hg38_peaks, value = TRUE)

corDf = data.frame()
for (i in 1:length(merip_hg38_peaks)) {
  peaks <- readRDS(merip_hg38_peaks[i])
  peaksInfo = as.data.frame(mcols(peaks))
  
  preTestPrs = cor.test(peaksInfo$log2FC, peaksInfo$RPM.input)
  preTestSpm = cor.test(peaksInfo$log2FC, peaksInfo$RPM.input, method = "spearman")
  
  preMetExpDf = data.frame(metLevel = peaksInfo$log2FC, expLevel = peaksInfo$RPM.input)
  prelm = lm(expLevel ~ metLevel, data = preMetExpDf)
  preR2 = summary(prelm)$r.squared

  cutoff = 0.1
  TP = motif_all[which(pred_df$mean_p0 <= cutoff),]
  overlapTP = findOverlaps(peaks, TP)
  
  highConfIndx = unique(queryHits(overlapTP))
  if (length(highConfIndx) < 2) {
    next
  }
  postTestPrs = cor.test(peaksInfo$log2FC[highConfIndx], peaksInfo$RPM.input[highConfIndx])
  postTestSpm = cor.test(peaksInfo$log2FC[highConfIndx], peaksInfo$RPM.input[highConfIndx], method = "spearman")
  
  postMetExpDf = data.frame(metLevel = peaksInfo$log2FC[highConfIndx], expLevel = peaksInfo$RPM.input[highConfIndx])
  postlm = lm(expLevel ~ metLevel, data = postMetExpDf)
  postR2 = summary(postlm)$r.squared
  
  
  tmp = data.frame(preCor = preTestPrs$estimate,
                   postCor = postTestPrs$estimate,
                   preR2 = preR2,
                   postR2 = postR2)
  corDf = rbind(corDf, tmp)
  print(i)
}


###################################################
####               plot all samples            ####
###################################################

library(patchwork)
getSciPval <- function(p_value) {
  scientific_p <- sprintf("%.1e", p_value) 
  split_p <- strsplit(scientific_p, "e")[[1]]
  p_main <- split_p[1]
  p_exp <- sub("^\\+0*", "", split_p[2])  # Remove leading zeros and plus
  p_exp <- sub("^\\-", "-", p_exp)  # Ensure the minus is preserved
  formatted_p <- paste0(p_main, " %*% 10^", p_exp)
  return(formatted_p)
}
createEcdfPlot <- function(data, x_var, x_label) {
  ggplot(data, aes_string(x = x_var, color = "corType")) + 
    stat_ecdf(geom = "step", size = 0.8) +
    scale_color_manual(values = c("#477ca8", "#cb3335"), labels = c("Pre-Calibration  ", "Post-Calibration")) +
    xlab(x_label) +
    ylab("Cumulative fraction") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0))) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 12),
          legend.position = "top",
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.box.margin = margin(-8, 0, -8, 0),
          legend.text = element_text(colour= "black", size = 10))
}
createBoxplotInset <- function(data, y_var, y_label, maxValue, shiftVal, minValue = 0) {
  ggplot(data, aes_string(x = "corType", y = y_var, fill = "corType")) + 
    geom_boxplot(outlier.colour = NA, width = 0.55) +
    scale_fill_manual(values = c("#477ca8", "#cb3335")) +
    stat_boxplot(geom = "errorbar", aes_string(ymin = "..ymax.."), width = 0.25, size = .3) +
    stat_boxplot(geom = "errorbar", aes_string(ymax = "..ymin.."), width = 0.25, size = .3) +
    labs(y = y_label) +
    annotate("segment", x = 1, xend = 1, y = maxValue + shiftVal, yend = maxValue + 2*shiftVal) +
    annotate("segment", x = 2, xend = 2, y = maxValue + shiftVal, yend = maxValue + 2*shiftVal) +
    annotate("segment", x = 1, xend = 2, y = maxValue + 2*shiftVal, yend = maxValue + 2*shiftVal) +
    coord_cartesian(ylim = c(minValue, maxValue - shiftVal), clip = "off") +
    theme_classic() +
    theme(axis.text = element_text(colour = "black", size = 9),
          axis.title = element_text(colour = "black", size = 11),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          plot.background = element_blank())
}

# box plot
corBoxDf <- reshape2::melt(corDf, measure.vars = c("preCor", "postCor"), 
                           id.vars = NULL, variable.name = "corType", value.name = "corValue")

corWilcoxTest = wilcox.test(corDf$preCor, corDf$postCor, alternative = "greater")
corMaxValue = round(boxplot(-corBoxDf$corValue[which(corBoxDf$corType == "postCor")])$stats[5, ], 2)

ecdf_plot <- createEcdfPlot(corBoxDf, "-corValue", "Correlation with gene expression (|cor|)")
boxplot_inset <- createBoxplotInset(corBoxDf, "-corValue", "|cor|", corMaxValue, 0.02) +
  annotate("text", x = 1.5, y = 0.74, label = getSciPval(corWilcoxTest$p.value), parse = TRUE, size = 3.3)

ecdf_plot + 
  inset_element(p = boxplot_inset, left = 0.45, bottom = 0, right = 1, top = 0.65)
ggsave("~/plots/negCorEcdfBox.pdf", width = 3.5, height = 3.5)


########################################
####          box plot MACS2        ####
########################################

# atlas MACS2
human <- readRDS("~/merip/human.rds")
peakData = human[which(human$CallPeak_Methods == "MACS2")]
conditions = unique(peakData$condition)
conditions = grep("Ctrl", conditions, value = TRUE)

corDf = data.frame()
for (i in 1:length(conditions)) {
  peaksInfo = peakData[which(peakData$GSE_CellLine_Treatment == conditions[i]),]
  peaks = makeGRangesFromDataFrame(peaksInfo)
  
  preTestPrs = cor.test(peaksInfo$log2FC, peaksInfo$TPM)
  preTestSpm = cor.test(peaksInfo$log2FC, peaksInfo$TPM, method = "spearman")
  
  preMetExpDf = data.frame(metLevel = peaksInfo$log2FC, expLevel = peaksInfo$TPM)
  prelm = lm(expLevel ~ metLevel, data = preMetExpDf)
  preR2 = summary(prelm)$r.squared
  
  cutoff = 0.1
  TP = motif_all[which(pred_df$mean_p0 <= cutoff),]
  overlapTP = findOverlaps(peaks, TP)
  
  highConfIndx = unique(queryHits(overlapTP))
  if (length(highConfIndx) < 2) {
    next
  }
  postTestPrs = cor.test(peaksInfo$log2FC[highConfIndx], peaksInfo$TPM[highConfIndx])
  postTestSpm = cor.test(peaksInfo$log2FC[highConfIndx], peaksInfo$TPM[highConfIndx], method = "spearman")
  
  postMetExpDf = data.frame(metLevel = peaksInfo$log2FC[highConfIndx], expLevel = peaksInfo$TPM[highConfIndx])
  postlm = lm(expLevel ~ metLevel, data = postMetExpDf)
  postR2 = summary(postlm)$r.squared
  
  
  tmp = data.frame(preCor = preTestPrs$estimate,
                   postCor = postTestPrs$estimate,
                   preR2 = preR2,
                   postR2 = postR2)
  corDf = rbind(corDf, tmp)
}

# box plot
corBoxDf <- reshape2::melt(corDf, measure.vars = c("preCor", "postCor"), 
                           id.vars = NULL, variable.name = "corType", value.name = "corValue")

corWilcoxTest = wilcox.test(corDf$preCor, corDf$postCor, alternative = "greater")
corMaxValue = round(boxplot(-corBoxDf$corValue[which(corBoxDf$corType == "postCor")])$stats[5, ], 2)

ecdf_plot <- createEcdfPlot(corBoxDf, "-corValue", "Correlation with gene expression (|cor|)")
boxplot_inset <- createBoxplotInset(corBoxDf, "-corValue", "|cor|", corMaxValue, 0.02, minValue = -0.1) +
  annotate("text", x = 1.5, y = 0.66, label = getSciPval(corWilcoxTest$p.value), parse = TRUE, size = 3.3)

ecdf_plot + 
  inset_element(p = boxplot_inset, left = 0.45, bottom = 0, right = 1, top = 0.65)
ggsave("~/plots/negCorMACS2EcdfBox.pdf", width = 3.5, height = 3.5)


###################################################
####               plot best sample            ####
###################################################

motif_all$p0 = pred_df$mean_p0

i = 44
peaks <- readRDS(merip_hg38_peaks[i])
peaksInfo = as.data.frame(mcols(peaks))

overlap = findOverlaps(peaks, motif_all)
peaksInfo$p0 = NA
peaksInfo$p0[queryHits(overlap)] = motif_all$p0[subjectHits(overlap)]
peaksInfo$p1 = 1 - peaksInfo$p0

preMetExpDf = data.frame(metLevel = peaksInfo$log2FC,
                         expLevel = peaksInfo$RPM.input,
                         p1 = peaksInfo$p1)
preMetExpDf = preMetExpDf[!is.na(preMetExpDf$p1),]
set.seed(123)
preIndx = sample(nrow(preMetExpDf), 1000)
preMetExpDf = preMetExpDf[preIndx,]

indx = which(preMetExpDf$metLevel >= 0 & preMetExpDf$metLevel <= 6 & log2(preMetExpDf$expLevel) >= -3 & log2(preMetExpDf$expLevel) <= 6.3)
ggplot(data = preMetExpDf, aes(x = metLevel, y = log2(expLevel), color = p1)) + 
  geom_point(alpha = 0.6) +
  scale_color_gradientn(colors = c("#2462ab", "#4a95c5", "#ffffff", "#b31a2f", "#62061d"),
                        values = c(0, 0.25, 0.5, 0.75, 1),
                        limits = c(0, 1)) +
  labs(x = expression(paste("log"[2], "(m6A level)")), y = expression(paste("log"[2], "(gene expression level)"))) +
  scale_x_continuous(breaks = c(1, 3, 5)) +
  scale_y_continuous(breaks = c(-2.5, 0, 2.5, 5)) +
  geom_smooth(data = preMetExpDf[indx,], method = "lm", se = FALSE, color = "blue") +
  geom_smooth(data = preMetExpDf[indx,][which(preMetExpDf[indx,]$p1 > 0.9),], method = "lm", se = FALSE, color = "red") +
  annotate("segment", x = c(0.4, 3.3), xend = c(0.9, 3.8), 
           y = 7.95, yend = 7.95, size = 0.8, color = c("blue", "red")) +
  annotate("text", x = c(1.0, 3.9), y = 7.95, 
           label = c("Pre-Calibration", "Post-Calibration"), 
           size = 3.7, hjust = 0, color = "black") +
  annotate("text", x = c(1.25, 4.2), y = 7.30, 
           label = c("cor = -0.17", "cor = -0.46"), 
           size = 3.7, hjust = 0, color = "black") +
  annotate("text", x = 7.55, y = 2.3, label = expression(paste("Predicted high-conf. m"^6, "A probability")), angle = -90, size = 3.7) +
  coord_cartesian(xlim = c(0.5, 5.6), ylim = c(-3.2, 6.3), clip = "off") +
  theme_classic() + 
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(1.10, "cm"),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.margin = margin(35, 25, 5, 5))
ggsave("~/plots/negCorDot.pdf", width = 4.8, height = 3.8)



################################################################################
####                           neg cor in diffPeak                          ####
################################################################################

target = c(1:16)
corDf = data.frame()
for (i in target) {
  peaks <- readRDS(paste0("~/merip/diffPeakAll/diffPeakAll_", i, "/diffPeaks.rds"))
  peaksInfo = as.data.frame(mcols(peaks))
  
  diffFC = peaksInfo$diff.log2FC
  diffExp = (peaksInfo$RPM.input.Treated + 0.001) / (peaksInfo$RPM.input.Control + 0.001)
  preTestPrs = cor.test(diffFC, diffExp)
  preTestSpm = cor.test(diffFC, diffExp, method = "spearman")
  
  preMetExpDf = data.frame(metLevel = diffFC, expLevel = diffExp)
  prelm = lm(expLevel ~ metLevel, data = preMetExpDf)
  preR2 = summary(prelm)$r.squared
  
  cutoff = 0.1
  TP = motif_all[which(pred_df$mean_p0 <= cutoff),]
  overlapTP = findOverlaps(peaks, TP)
  
  highConfIndx = unique(queryHits(overlapTP))
  if (length(highConfIndx) < 2) {
    next
  }
  postTestPrs = cor.test(diffFC[highConfIndx], diffExp[highConfIndx])
  postTestSpm = cor.test(diffFC[highConfIndx], diffExp[highConfIndx], method = "spearman")
  
  postMetExpDf = data.frame(metLevel = diffFC[highConfIndx], expLevel = diffExp[highConfIndx])
  postlm = lm(expLevel ~ metLevel, data = postMetExpDf)
  postR2 = summary(postlm)$r.squared
  
  
  tmp = data.frame(preCor = preTestPrs$estimate,
                   postCor = postTestPrs$estimate,
                   preR2 = preR2,
                   postR2 = postR2)
  corDf = rbind(corDf, tmp)
  print(i)
}

###################################################
####               plot all samples            ####
###################################################
########################################
####          dot box plot cor      ####
########################################

# dot plot
corDf$diff = corDf$postCor - corDf$preCor
dot_plot = ggplot(data = corDf, aes(x = -preCor, y = -postCor, color = -diff, size = -diff)) +
  geom_point(color = "#b31a2f") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#477ca8", size = 0.8) +
  xlab("Pre-calibration correlation |cor|") +
  ylab("Post-calibration correlation |cor|") +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  theme_classic() + 
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12),
        legend.position = "none")


# box plot
corBoxDf <- reshape2::melt(corDf, measure.vars = c("preCor", "postCor"), 
                           id.vars = NULL, variable.name = "corType", value.name = "corValue")

corWilcoxTest = wilcox.test(corDf$preCor, corDf$postCor, alternative = "greater")
corMaxValue = round(boxplot(-corBoxDf$corValue[which(corBoxDf$corType == "postCor")])$stats[5, ], 2)
maxValue = corMaxValue
shiftVal = 0.02
minValue = 0.06
boxplot_inset = ggplot(corBoxDf, aes(x = corType, y = -corValue, fill = corType)) + 
  geom_boxplot(outlier.colour = NA, width = 0.55) +
  scale_fill_manual(values = c("#477ca8", "#cb3335")) +
  stat_boxplot(geom = "errorbar", aes_string(ymin = "..ymax.."), width = 0.25, size = .3) +
  stat_boxplot(geom = "errorbar", aes_string(ymax = "..ymin.."), width = 0.25, size = .3) +
  labs(y = "|cor|") +
  annotate("segment", x = 1, xend = 1, y = maxValue + shiftVal, yend = maxValue + 2*shiftVal) +
  annotate("segment", x = 2, xend = 2, y = maxValue + shiftVal, yend = maxValue + 2*shiftVal) +
  annotate("segment", x = 1, xend = 2, y = maxValue + 2*shiftVal, yend = maxValue + 2*shiftVal) +
  annotate("text", x = 1.5, y = 0.85, label = getSciPval(corWilcoxTest$p.value), parse = TRUE, size = 3.3) +
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
ggsave("~/plots/negCorDiffDotBox.pdf", width = 3.5, height = 3.5)


########################################
####            dumbbell cor        ####
########################################

data_dir <- "~/merip/Data/"
writer_list <- c("KIAA1429", "METTL3", "METTL14", "WTAP", "METTL16", "M3M4")
peak_dir_all <- grep("PEAK", list.files(data_dir), value = TRUE)
peak_dir <- unlist(lapply(writer_list, function(writer) {
  grep(writer, peak_dir_all, value = TRUE)
}))
peak_dir = peak_dir[target]

peaksWriter = sapply(strsplit(peak_dir, "_"), function(peak) strsplit(peak, "_")[[3]][1])
peaksTissue = sapply(strsplit(peak_dir, "_"), function(peak) strsplit(peak, "_")[[2]][1])

corDf$writer = peaksWriter
corDf$tissue = peaksTissue
corDf$info = paste0(peaksWriter, " KD - ", peaksTissue)
corDf$writer <- ave(corDf$writer, corDf$writer, FUN = function(x) {
  if(length(x) > 1) paste0(x, "-", seq_along(x)) else x
})
corDf = corDf[order(corDf$writer),]



corDf2 <- reshape2::melt(corDf, measure.vars = c("preCor", "postCor"), 
                           id.vars = "writer", variable.name = "corType", value.name = "corValue")
ggplot(corDf2, aes(x = -corValue, y = writer)) +
  geom_line(aes(group = writer), color = "#fbdec8", size = 0.8, alpha = 0.8) +
  geom_point(aes(color = corType), size = 2) +
  xlab("Correlation with gene expression (|cor|)") +
  scale_y_discrete(breaks = corDf$writer, labels = corDf$info) + 
  scale_color_manual(breaks = c("preCor", "postCor"),
                     values = c("#477ca8", "#cb3335"), labels = c("Pre-Calibration   ", "Post-Calibration")) +
  theme_classic() +
  theme(axis.text = element_text(colour= "black", size = 10),
        axis.title = element_text(colour= "black", size = 11),
        axis.title.y = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.05, "cm"),
        legend.background = element_blank(),
        legend.margin = margin(5, 5, 2, -6),
        legend.box.margin = margin(-8, 0, -8, 0),
        legend.text = element_text(colour= "black", size = 10))

ggsave("~/plots/negCorDiffDumbbell.pdf", width = 4.5, height = 3.7)





###################################################
####               plot best sample            ####
###################################################

i = 14
peaks <- readRDS(paste0("~/merip/diffPeakAll/diffPeakAll_", i, "/diffPeaks.rds"))
peaksInfo = as.data.frame(mcols(peaks))

diffFC = peaksInfo$diff.log2FC
diffExp = (peaksInfo$RPM.input.Treated + 0.001) / (peaksInfo$RPM.input.Control + 0.001)

motif_all$p0 = pred_df$mean_p0
overlap = findOverlaps(peaks, motif_all)
peaksInfo$p0 = NA
peaksInfo$p0[queryHits(overlap)] = motif_all$p0[subjectHits(overlap)]
peaksInfo$p1 = 1 - peaksInfo$p0

preMetExpDf = data.frame(metLevel = diffFC,
                         expLevel = diffExp,
                         p1 = peaksInfo$p1)
preMetExpDf = preMetExpDf[!is.na(preMetExpDf$p1),]
set.seed(1234)
preIndx = sample(nrow(preMetExpDf), 1000)
preMetExpDf = preMetExpDf[preIndx,]

indx = which(preMetExpDf$metLevel >= -5.35 & preMetExpDf$metLevel <= 2.5 & log2(preMetExpDf$expLevel) >= -1.65 & log2(preMetExpDf$expLevel) <= 3.5)
ggplot(data = preMetExpDf, aes(x = metLevel, y = log2(expLevel), color = p1)) + 
  geom_point(alpha = 0.6) +
  scale_color_gradientn(colors = c("#2462ab", "#4a95c5", "#ffffff", "#b31a2f", "#62061d"),
                        values = c(0, 0.25, 0.5, 0.75, 1),
                        name = "",
                        limits = c(0, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#707070") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#707070") +
  labs(x = expression(paste("log"[2], "FC of m6A level after KIAA1429 KD")), y = expression(paste("log"[2], "FC of gene expression level (FPKM)"))) +
  geom_smooth(data = preMetExpDf[indx,], method = "lm", se = FALSE, color = "blue") +
  geom_smooth(data = preMetExpDf[indx,][which(preMetExpDf[indx,]$p1 > 0.9),], method = "lm", se = FALSE, color = "red") +
  annotate("segment", x = c(-5.4, -1.25), xend = c(-4.7, -0.55), 
           y = 4.4, yend = 4.4, size = 0.8, color = c("blue", "red")) +
  annotate("text", x = c(-4.55, -0.4), y = 4.4, 
           label = c("Pre-Calibration", "Post-Calibration"), 
           size = 3.7, hjust = 0, color = "black") +
  annotate("text", x = c(-4.15, 0.10), y = 4.02, 
           label = c("cor = -0.19", "cor = -0.52"), 
           size = 3.7, hjust = 0, color = "black") +
  annotate("text", x = 5.4, y = 1.33, label = expression(paste("Predicted high-conf. m"^6, "A probability")), angle = -90, size = 3.7) +
  coord_cartesian(xlim = c(-5.35, 2.5), ylim = c(-1.65, 3.5), clip = "off") +
  theme_classic() + 
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(1.10, "cm"),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.margin = margin(35, 25, 5, 5))
ggsave("~/plots/corDiffDot.pdf", width = 4.8, height = 3.8)



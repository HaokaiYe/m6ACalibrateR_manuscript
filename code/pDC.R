
################################################################################
####                 proportion of directional consistency                  ####
################################################################################

library(scales)
library(ggplot2)
library(patchwork)
library(SummarizedExperiment)
source("~/reviewerComments/surf/modules_v0.3.R")

motif_all <- readRDS("~/rds/motif_all.rds") - 2
pred_df <- readRDS("~/rds/pred_df.rds")
TP = motif_all[which(pred_df$mean_p0 <= 0.5),]

data_dir <- "~/merip/Data/"
writer_list <- c("KIAA1429", "METTL3", "METTL14", "WTAP", "METTL16", "M3M4")
peak_dir_all <- grep("PEAK", list.files(data_dir), value = TRUE)
peak_dir <- unlist(lapply(writer_list, function(writer) {
  grep(writer, peak_dir_all, value = TRUE)
}))

##############################################################
####                    exomePeak2 data                   ####
##############################################################

pDCDf = data.frame()
target = c(1:16)
for (i in target) {
  diffPeakData <- read.csv(paste0("~/merip/diffPeakAll/diffPeakAll_", i, "/diffPeaks.csv"))
  diffPeakGrgList <- readRDS(paste0("~/merip/diffPeakAll/diffPeakAll_", i, "/diffPeaks.rds"))
  
  lowCovIndx = which(diffPeakData$RPM.IP.Control < 1 | diffPeakData$RPM.input.Control < 1
                     | diffPeakData$RPM.IP.Treated < 1 | diffPeakData$RPM.input.Treated < 1)
  diffPeakData = diffPeakData[-lowCovIndx,]
  diffPeakGrgList = diffPeakGrgList[-lowCovIndx]
  
  overlapTP = findOverlaps(diffPeakGrgList, TP)
  highConfIndx = unique(queryHits(overlapTP))
  
  lfc = diffPeakData$diff.log2FC
  pDC = sum(lfc < 0) / length(diffPeakGrgList)
  pDC_Cali = sum(lfc[highConfIndx] < 0) / length(highConfIndx)
  
  tmp = data.frame(pDC = pDC,
                   pDC_Cali= pDC_Cali)
  
  pDCDf = rbind(pDCDf, tmp)
  print(i)
}


###################################################
####               plot all samples            ####
###################################################
########################################
####          dot box plot pDC      ####
########################################

# dot plot
pDCDf$diff = pDCDf$pDC_Cali - pDCDf$pDC
dot_plot = ggplot(data = pDCDf, aes(x = pDC, y = pDC_Cali, size = diff)) +
  geom_point(color = "#b31a2f") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#477ca8", size = 0.8) +
  xlab("Pre-calibration consistency (%)") +
  ylab("Post-calibration consistency (%)") +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1), expand = c(0, 0), labels = percent_format(suffix = "")) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1), expand = c(0, 0), labels = percent_format(suffix = "")) +
  theme_classic() + 
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12),
        legend.position = "none")


# box plot
pDCBoxDf <- reshape2::melt(pDCDf, measure.vars = c("pDC", "pDC_Cali"),
                           id.vars = NULL, variable.name = "corType", value.name = "corValue")

wilcoxTest = wilcox.test(pDCDf$pDC_Cali, pDCDf$pDC, alternative = "greater")
maxValue = round(boxplot(pDCBoxDf$corValue[which(pDCBoxDf$corType == "pDC_Cali")])$stats[5, ], 2)
maxValue = maxValue + 0.02
shiftVal = 0.02
minValue = 0.01
boxplot_inset = ggplot(pDCBoxDf, aes(x = corType, y = corValue, fill = corType)) + 
  geom_boxplot(outlier.colour = NA, width = 0.55) +
  scale_fill_manual(values = c("#477ca8", "#cb3335")) +
  stat_boxplot(geom = "errorbar", aes_string(ymin = "..ymax.."), width = 0.25, size = .3) +
  stat_boxplot(geom = "errorbar", aes_string(ymax = "..ymin.."), width = 0.25, size = .3) +
  labs(y = "Consistency (%)") +
  annotate("segment", x = 1, xend = 1, y = maxValue + shiftVal, yend = maxValue + 2*shiftVal) +
  annotate("segment", x = 2, xend = 2, y = maxValue + shiftVal, yend = maxValue + 2*shiftVal) +
  annotate("segment", x = 1, xend = 2, y = maxValue + 2*shiftVal, yend = maxValue + 2*shiftVal) +
  annotate("text", x = 1.5, y = 1.12, label = getSciPval(wilcoxTest$p.value), parse = TRUE, size = 3.3) +
  coord_cartesian(ylim = c(minValue, maxValue - shiftVal), clip = "off") +
  scale_x_discrete(labels = c("Pre", "Post")) + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = percent_format(suffix = "")) + 
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 9),
        axis.title = element_text(colour = "black", size = 11),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.background = element_blank(),
        plot.margin = margin(10, 5, 5, 5))

dot_plot + 
  inset_element(p = boxplot_inset, left = 0.48, bottom = 0, right = 1.02, top = 0.58)
ggsave("~/plots/pDCDotBox.pdf", width = 3.5, height = 3.5)



########################################
####           dumbbell pDC         ####
########################################

peaksWriter = sapply(strsplit(peak_dir, "_"), function(peak) strsplit(peak, "_")[[3]][1])
peaksTissue = sapply(strsplit(peak_dir, "_"), function(peak) strsplit(peak, "_")[[2]][1])

pDCDf$writer = peaksWriter
pDCDf$tissue = peaksTissue
pDCDf$info = paste0(peaksWriter, " KD - ", peaksTissue)
pDCDf$writer <- ave(pDCDf$writer, pDCDf$writer, FUN = function(x) {
  if(length(x) > 1) paste0(x, "-", seq_along(x)) else x
})
pDCDf = pDCDf[order(pDCDf$writer),]

pDCDf2 <- reshape2::melt(pDCDf, measure.vars = c("pDC", "pDC_Cali"), 
                         id.vars = "writer", variable.name = "corType", value.name = "corValue")
ggplot(pDCDf2, aes(x = corValue, y = writer)) +
  geom_line(aes(group = writer), color = "#fbdec8", size = 0.8, alpha = 0.8) +
  geom_point(aes(color = corType), size = 2) +
  xlab("Directional consistency (%)") +
  scale_y_discrete(breaks = pDCDf$writer, labels = pDCDf$info) + 
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = percent_format(suffix = "")) + 
  scale_color_manual(breaks = c("pDC", "pDC_Cali"),
                     values = c("#477ca8", "#cb3335"), labels = c("Pre-Calibration   ", "Post-Calibration")) +
  theme_classic() +
  theme(axis.text = element_text(colour= "black", size = 10),
        axis.title = element_text(colour= "black", size = 12),
        axis.title.y = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.05, "cm"),
        legend.background = element_blank(),
        legend.margin = margin(5, 5, 2, -6),
        legend.box.margin = margin(-8, 0, -8, 0),
        legend.text = element_text(colour= "black", size = 10))

ggsave("~/plots/pDCDumbbell.pdf", width = 4.5, height = 3.7)



###################################################
####               plot best sample            ####
###################################################

# pDC curve
library(scales)
TP = motif_all[which(pred_df$mean_p0 <= 0.5),]
getDf_pDC = function(i, covCut = 1) {
  diffPeakData <- read.csv(paste0("~/diffPeakAll/diffPeakAll_", i, "/diffPeaks.csv"))
  diffPeakGrgList <- readRDS(paste0("~/diffPeakAll/diffPeakAll_", i, "/diffPeaks.rds"))
  
  lowCovIndx = which(diffPeakData$RPM.IP.Control < covCut | diffPeakData$RPM.input.Control < covCut
                     | diffPeakData$RPM.IP.Treated < covCut | diffPeakData$RPM.input.Treated < covCut)
  diffPeakData = diffPeakData[-lowCovIndx,]
  diffPeakGrgList = diffPeakGrgList[-lowCovIndx]
  
  # TP
  overlapTP = findOverlaps(diffPeakGrgList, TP)
  highConfIndx = unique(queryHits(overlapTP))
  
  top_sites_used = length(highConfIndx)
  # before
  pval_pre <- abs(diffPeakData$diff.log2FC)
  lfc_pre <- diffPeakData$diff.log2FC
  DC <- cumsum((lfc_pre < 0)[order(pval_pre, decreasing = T)][seq_len(top_sites_used)])
  prop_DC <- DC/seq_along(DC)
  
  # after
  pval_Cali <- pval_pre[highConfIndx]
  lfc_Cali <- lfc_pre[highConfIndx]
  DC_Cali <- cumsum((lfc_Cali < 0)[order(pval_Cali, decreasing = T)][seq_len(top_sites_used)])
  prop_DC_Cali <- DC_Cali/seq_along(DC_Cali)
  
  # Create dataframe for the plot
  plot_df <- data.frame(pDC = c(prop_DC, prop_DC_Cali),
                        topSiteNum = c(seq_along(DC), seq_along(DC_Cali)),
                        class = rep(c("Before", "After"), c(length(DC), length(DC_Cali))))
  return(plot_df)
}
pDC_curveDf = getDf_pDC(1, covCut = 1)

i = 1
split_title <- strsplit(peak_dir[i], "_")[[1]]
sampleTitle <- paste(split_title[2], "Cells with", split_title[3], "Knockdown")

ggplot(data = pDC_curveDf) + 
  geom_path(aes(y = pDC, x = topSiteNum, colour = class)) +
  scale_color_manual(values = c("#0076ba", "#cb3335"),
                     name = NULL,
                     breaks = c("Before", "After"),
                     labels = c("Pre-Calibration", "Post-Calibration")) +
  scale_y_continuous(breaks = seq(0.92, 1, by = 0.02),
                     labels = scales::percent_format(suffix = "")) +
  scale_x_continuous(labels = comma) +
  theme_classic() +  
  labs(x = "Top differentially methylated peaks", y = "Directional consistency (%)", title = sampleTitle) +
  theme(axis.text = element_text(colour= "black", size = 10),
        axis.title = element_text(colour= "black", size = 12),
        legend.position = c(0.3, 0.16),
        legend.title = element_blank(),
        legend.background = element_blank(), 
        legend.text = element_text(colour = "black", size = 10),
        plot.title = element_text(colour = "black", size = 11, hjust = 0.5))

ggsave("~/plots/pDCCurve.pdf", width = 3.3, height = 3.2)





################################################################################
####                compare m6ACali with other prediction tools             ####
################################################################################

# whistle
library(rtracklayer)
hg19ToHg38Chain = import.chain("~/hg19ToHg38.over.chain")
whistle_csv = read.csv("~/PredictedByMatureModel.csv")
whistle_hg19 = makeGRangesFromDataFrame(whistle_csv)
whistle_hg19$whistleProb = whistle_csv$probability
whistle_hg38 <- unlist(liftOver(whistle_hg19, hg19ToHg38Chain))

# iM6A
iM6A <- readRDS("~/iM6A.rds")

# TP from other predictor
motif_all$m6ACali = 1 - pred_df$mean_p0

m6ACaliTP = motif_all[motif_all$m6ACali >= median(motif_all$m6ACali)]
whistleTP = whistle_hg38[whistle_hg38$whistleProb >= median(whistle_hg38$whistleProb)]
iM6ATP = iM6A[iM6A$iM6AProb >= median(iM6A$iM6AProb)]

TPList <- list(m6ACaliTP = m6ACaliTP, whistleTP = whistleTP, iM6ATP = iM6ATP)

##############################################################
####                  different threshold                 ####
##############################################################

getTPSubset <- function(gr, col, quantilePercentage) {
  threshold <- quantile(as.numeric(mcols(gr)[[col]]), quantilePercentage)
  return(gr[mcols(gr)[[col]] >= threshold])
}

thresholds <- seq(0.2, 0.8, 0.2)

getAllpDC = function(percentage, covCut) {
  pDCDf <- data.frame()
  target = setdiff(c(1:19), c(3, 7, 12))
  
  # new TPList
  m6ACaliTP <- getTPSubset(motif_all, "m6ACali", percentage)
  whistleTP <-getTPSubset(whistle_hg38, "whistleProb", percentage)
  iM6ATP <- getTPSubset(iM6A, "iM6AProb", percentage)
  TPList <- list(m6ACaliTP = m6ACaliTP, whistleTP = whistleTP, iM6ATP = iM6ATP)
  
  for (i in target) {
    diffPeakData <- read.csv(paste0("~/diffPeakAll/diffPeakAll_", i, "/diffPeaks.csv"))
    diffPeakGrgList <- readRDS(paste0("~/diffPeakAll/diffPeakAll_", i, "/diffPeaks.rds"))
    
    lowCovIndx = which(diffPeakData$RPM.IP.Control < covCut | diffPeakData$RPM.input.Control < covCut
                       | diffPeakData$RPM.IP.Treated < covCut | diffPeakData$RPM.input.Treated < covCut)
    diffPeakData = diffPeakData[-lowCovIndx,]
    diffPeakGrgList = diffPeakGrgList[-lowCovIndx]
    
    lfc = diffPeakData$diff.log2FC
    pDC = sum(lfc < 0) / length(diffPeakGrgList)
    
    for (TP_name in names(TPList)) {
      TP <- TPList[[TP_name]]
      
      overlapTP = findOverlaps(diffPeakGrgList, TP)
      highConfIndx = unique(queryHits(overlapTP))
      pDC_Cali = sum(lfc[highConfIndx] < 0) / length(highConfIndx)
      
      tmp = data.frame(Sample = i, Threshold = percentage, Type = TP_name, Value = pDC_Cali)
      pDCDf = rbind(pDCDf, tmp)
    }
    
    pDCDf = rbind(pDCDf, data.frame(Sample = i, Threshold = percentage, Type = "pDC", Value = pDC))
    print(i)
  }
  pDCDf$Type = factor(pDCDf$Type, levels = c("pDC", "iM6ATP", "whistleTP", "m6ACaliTP"))
  return(pDCDf)
}
pDCList <- lapply(thresholds, function(x) getAllpDC(x, covCut = 1))
pDCDf <- do.call(rbind, pDCList)

pDCDf$Type = factor(pDCDf$Type, levels = c("pDC", "iM6ATP", "whistleTP", "m6ACaliTP"))
pDCDf$Threshold = factor(pDCDf$Threshold, levels = rev(thresholds))
# 绘制平均值线图
ggplot(pDCDf[which(pDCDf$Type != "pDC"),], aes(x = Threshold, y = Value, color = Type, group = Type)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point") +
  geom_hline(yintercept = 0.46, linetype = "dashed", color = "#0076ba", size = 0.8) +
  scale_color_manual(values = rev(c("#7fc3ab", "#ee804f", "#8b1643")),
                     name = "Calibration tools",
                     breaks = rev(c("iM6ATP", "whistleTP", "m6ACaliTP")),
                     labels = rev(c("iM6A", "WHISTLE", "m6ACali"))) +
  scale_x_discrete(breaks = seq(0.2, 0.8, 0.2), labels = rev(seq(20, 80, 20))) +
  labs(x = "Calibration threshold (%)", y = "Average pDC") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12),
        legend.title = element_text(colour= "black", size = 11, hjust = 0.5),
        legend.background = element_blank(),
        legend.box.margin = margin(0, -5, 0, -15),
        legend.text = element_text(colour = "black", size = 10))
ggsave("~/plots/pDCCompareMean.pdf", width = 4.6, height = 2.65)




##############################################################
####                      all samples                     ####
##############################################################

library(scales)
getAllpDC = function(percentage, covCut) {
  pDCDf <- data.frame()
  target = c(1:16)
  
  # new TPList
  m6ACaliTP <- getTPSubset(motif_all, "m6ACali", percentage)
  whistleTP <-getTPSubset(whistle_hg38, "whistleProb", percentage)
  iM6ATP <- getTPSubset(iM6A, "iM6AProb", percentage)
  TPList <- list(m6ACaliTP = m6ACaliTP, whistleTP = whistleTP, iM6ATP = iM6ATP)
  
  for (i in target) {
    diffPeakData <- read.csv(paste0("~/diffPeakAll/diffPeakAll_", i, "/diffPeaks.csv"))
    diffPeakGrgList <- readRDS(paste0("~/diffPeakAll/diffPeakAll_", i, "/diffPeaks.rds"))
    
    lowCovIndx = which(diffPeakData$RPM.IP.Control < covCut | diffPeakData$RPM.input.Control < covCut
                       | diffPeakData$RPM.IP.Treated < covCut | diffPeakData$RPM.input.Treated < covCut)
    diffPeakData = diffPeakData[-lowCovIndx,]
    diffPeakGrgList = diffPeakGrgList[-lowCovIndx]
    
    lfc = diffPeakData$diff.log2FC
    pDC = sum(lfc < 0) / length(diffPeakGrgList)
    
    for (TP_name in names(TPList)) {
      TP <- TPList[[TP_name]]
      
      overlapTP = findOverlaps(diffPeakGrgList, TP)
      highConfIndx = unique(queryHits(overlapTP))
      pDC_Cali = sum(lfc[highConfIndx] < 0) / length(highConfIndx)
      
      tmp = data.frame(Sample = i, Threshold = percentage, Type = TP_name, Value = pDC_Cali)
      pDCDf = rbind(pDCDf, tmp)
    }
    
    pDCDf = rbind(pDCDf, data.frame(Sample = i, Threshold = percentage, Type = "pDC", Value = pDC))
    print(i)
  }
  pDCDf$Type = factor(pDCDf$Type, levels = c("pDC", "iM6ATP", "whistleTP", "m6ACaliTP"))
  return(pDCDf)
}
pDCDf = getAllpDC(0.5, covCut = 1)

# box plot
ggplot(pDCDf, aes(x = Type, y = Value, fill = Type)) + 
  geom_boxplot(outlier.colour = NA, width = 0.75) +
  scale_fill_manual(values = rev(c("#0076ba", "#7fc3ab", "#ee804f", "#8b1643")),
                     name = "Calibration tools",
                     breaks = rev(c("pDC", "iM6ATP", "whistleTP", "m6ACaliTP")),
                     labels = rev(c("Baseline", "iM6A", "WHISTLE", "m6ACali"))) +
  stat_boxplot(geom = "errorbar",
               aes(ymin = ..ymax..),
               width = 0.25, size = .3) +
  stat_boxplot(geom = "errorbar",
               aes(ymax = ..ymin..),
               width = 0.25, size = .3) +
  scale_y_continuous(labels = scales::percent_format(suffix = "")) +
  scale_x_discrete(breaks = c("pDC", "iM6ATP", "whistleTP", "m6ACaliTP"), labels = c("Baseline", "iM6A", "WHISTLE", "m6ACali")) +
  labs(x = "Calibration tools", y = "Directional consistency (%)") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12),
        legend.position = "none",
        plot.background = element_blank())
ggsave("~/plots/pDCCompareBox.pdf", width = 3.3, height = 3.3)



##############################################################
####                      best sample                     ####
##############################################################

compareDf_pDC = function(i, covCut = 1, TPList = TPList) {
  diffPeakData <- read.csv(paste0("~/diffPeakAll/diffPeakAll_", i, "/diffPeaks.csv"))
  diffPeakGrgList <- readRDS(paste0("~/diffPeakAll/diffPeakAll_", i, "/diffPeaks.rds"))
  
  lowCovIndx = which(diffPeakData$RPM.IP.Control < covCut | diffPeakData$RPM.input.Control < covCut
                     | diffPeakData$RPM.IP.Treated < covCut | diffPeakData$RPM.input.Treated < covCut)
  diffPeakData <- diffPeakData[-lowCovIndx,]
  diffPeakGrgList <- diffPeakGrgList[-lowCovIndx]
  
  # Calculate the minimum length of highConfIndx across all predictors
  highConfIndxLengths <- sapply(TPList, function(TP) {
    overlapTP <- findOverlaps(diffPeakGrgList, TP)
    length(unique(queryHits(overlapTP)))
  })
  minHighConfIndxLength <- min(highConfIndxLengths)
  
  # Calculate 'Before' once for all predictors
  pval_pre <- abs(diffPeakData$diff.log2FC)
  lfc_pre <- diffPeakData$diff.log2FC
  
  DC <- cumsum((lfc_pre < 0)[order(pval_pre, decreasing = TRUE)])
  prop_DC <- DC / seq_along(DC)
  
  before_df <- data.frame(
    pDC = prop_DC[seq_len(minHighConfIndxLength)],
    topSiteNum = seq_len(minHighConfIndxLength),
    class = rep("Before", minHighConfIndxLength)
  )
  
  after_dfs <- lapply(names(TPList), function(predictorName) {
    TP <- TPList[[predictorName]]
    
    # TP
    overlapTP <- findOverlaps(diffPeakGrgList, TP)
    highConfIndx <- unique(queryHits(overlapTP))
    
    pval_Cali <- pval_pre[highConfIndx]
    lfc_Cali <- lfc_pre[highConfIndx]
    DC_Cali <- cumsum((lfc_Cali < 0)[order(pval_Cali, decreasing = TRUE)])
    prop_DC_Cali <- DC_Cali / seq_along(DC_Cali)
    
    after_df <- data.frame(
      pDC = prop_DC_Cali[seq_len(minHighConfIndxLength)],
      topSiteNum = seq_len(minHighConfIndxLength),
      class = rep(paste0("After_", predictorName), minHighConfIndxLength)
    )
    
    return(after_df)
  })
  
  # Combine 'Before' dataframe with all 'After' dataframes
  final_plot_df <- rbind(before_df, do.call(rbind, after_dfs))
  return(final_plot_df)
}

# Example usage
pDC_compareDf <- compareDf_pDC(i = 1, covCut = 1, TPList = TPList)

i = 1
split_title <- strsplit(peak_dir[i], "_")[[1]]
sampleTitle <- paste(split_title[2], "Cells with", split_title[3], "Knockdown")

ggplot(data = pDC_compareDf) + 
  geom_path(aes(y = pDC, x = topSiteNum, colour = class)) +
  scale_color_manual(values = rev(c("#0076ba", "#7fc3ab", "#ee804f", "#8b1643")),
                     name = "Calibration tools",
                     breaks = rev(c("Before", "After_iM6ATP", "After_whistleTP", "After_m6ACaliTP")),
                     labels = rev(c("Baseline", "iM6A", "WHISTLE", "m6ACali"))) +
  scale_y_continuous(labels = scales::percent_format(suffix = "")) +
  scale_x_continuous(labels = comma) +
  theme_classic() +  
  labs(x = "Top differentially methylated peaks", y = "Directional consistency (%)", title = sampleTitle) +
  theme(#panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.text = element_text(colour= "black", size = 10),
    axis.title = element_text(colour= "black", size = 12),
    legend.position = c(0.3, 0.28),
    legend.title = element_text(colour= "black", size = 11, hjust = 0.5),
    legend.background = element_blank(), 
    legend.text = element_text(colour = "black", size = 10),
    plot.title = element_text(colour = "black", size = 11, hjust = 0.5))

ggsave("~/plots/pDCCompare.pdf", width = 3.3, height = 3.2)






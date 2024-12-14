
library(ggplot2)
library(patchwork)
library(SummarizedExperiment)

#homoSE <- readRDS("~/rds/homoSE.rds")
homoSE_IP <- readRDS("~/rds/homoSE_IP.rds")
homoSE_input <- readRDS("~/rds/homoSE_input.rds")

################################################################################
####                          variance of m6A level                         ####
################################################################################

IPCount <- assays(homoSE_IP)[["IP"]]
inputCount <- assays(homoSE_input)[["input"]]

# count + 1
IPCountNorm = t(t(IPCount + 1)/colMeans(IPCount))
inputCountNorm = t(t(inputCount + 1)/colMeans(inputCount))
m6ALevel = log(IPCountNorm/inputCountNorm)

homeGrg = rowRanges(homoSE)
homeGrg$variance = rowVars(m6ALevel)

varianceDf = data.frame(variance = c(homeGrg$variance[homeGrg$p0 <= 0.5],
                                     homeGrg$variance[homeGrg$p0 > 0.5],
                                     homeGrg$variance[homeGrg$p0 > 0.9]),
                        class = c(rep("TP", sum(homeGrg$p0 <= 0.5)),
                                  rep("FP", sum(homeGrg$p0 > 0.5)),
                                  rep("Top FP", sum(homeGrg$p0 > 0.9))))
varianceDf$class = factor(varianceDf$class, levels = c("TP", "FP", "Top FP"))

##wilcoxTest1 = wilcox.test(homeGrg$variance[homeGrg$p0 <= 0.5], homeGrg$variance[homeGrg$p0 > 0.5], alternative = "less")
##wilcoxTest2 = wilcox.test(homeGrg$variance[homeGrg$p0 <= 0.5], homeGrg$variance[homeGrg$p0 > 0.9], alternative = "less")

densityPlot = ggplot(varianceDf) + 
  geom_density(aes(x = variance, fill = class), alpha = 0.7, color = NA) +
  scale_fill_manual(values = c("#f7a698", "#69a4d2", "#002D52"),
                    breaks = c("TP", "FP", "Top FP"), 
                    labels = c(expression(paste("High-conf. m"^6, "A   ")), "Off-targets   ", "Strong off-targets")) +
  labs(x = "Variance of methylation level", y = "Density") +
  xlim(0, 5) +
  theme_classic() +
  theme(axis.text = element_text(colour= "black", size = 10),
        axis.title = element_text(colour= "black", size = 12),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.spacing.x = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_blank(), 
        legend.text = element_text(colour = "black", size = 10))

boxPlot = ggplot(varianceDf, aes(x = class, y = variance, fill = class)) + 
  geom_boxplot(outlier.colour = NA, width = 0.55, alpha = 0.7) +
  scale_fill_manual(values = c("#f7a698", "#69a4d2", "#002D52"),
                     breaks = c("TP", "FP", "Top FP"), 
                     labels = c(expression(paste("High-conf. m"^6, "A   ")), "Off-targets   ", "Strong off-targets")) +
  stat_boxplot(geom = "errorbar",
               aes(ymin = ..ymax..),
               width = 0.25, size = .3) +
  stat_boxplot(geom = "errorbar",
               aes(ymax = ..ymin..),
               width = 0.25, size = .3) +
  labs(y = "Variance of methylation level") +
  ylim(0, 6) +
  annotate("segment", x = 1, xend = 1, y = 5.0, yend = 5.1) +
  annotate("segment", x = 2, xend = 2, y = 5.0, yend = 5.1) +
  annotate("segment", x = 1, xend = 2, y = 5.1, yend = 5.1) +
  annotate("text", x = 1.5, y = 5.2, label = "****", size = 5.2) +
  annotate("segment", x = 1, xend = 1, y = 5.7, yend = 5.8) +
  annotate("segment", x = 3, xend = 3, y = 5.7, yend = 5.8) +
  annotate("segment", x = 1, xend = 3, y = 5.8, yend = 5.8) +
  annotate("text", x = 2, y = 5.9, label = "****", size = 5.2) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.background = element_blank())

densityBoxPlot = densityPlot + boxPlot + plot_layout(widths = c(2.0, 1))
guide_area()/densityBoxPlot + plot_layout(heights = c(1, 10), guides = 'collect')

ggsave("~/plots/metLevelVar.pdf", width = 5.8, height = 3.3)



ggplot(varianceDf, aes(x = class, y = variance, fill = class)) + 
  geom_boxplot(outlier.colour = NA, width = 0.7, alpha = 0.7) +
  scale_fill_manual(values = c("#f7a698", "#69a4d2", "#002D52"),
                    breaks = c("TP", "FP", "Top FP")) +
  stat_boxplot(geom = "errorbar",
               aes(ymin = ..ymax..),
               width = 0.25, size = .3) +
  stat_boxplot(geom = "errorbar",
               aes(ymax = ..ymin..),
               width = 0.25, size = .3) +
  labs(y = "Variance of methylation level") +
  scale_x_discrete(breaks = c("TP", "FP", "Top FP"),
                   labels = c(expression(paste("High-conf. m"^6, "A")), "Off-target", "Strong off-target")) +
  scale_y_continuous(breaks = c(0, 2, 4), limits = c(0, 4.7)) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12),
        axis.text.x = element_text(vjust = 0.65, angle = 20),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.background = element_blank())

ggsave("~/plots/metLevelBox.pdf", width = 3.3, height = 3.3)


################################################################################
####                  selected top variance m6A sites                       ####
################################################################################

filterSE = function (RM_se, missing_threshold_IP = 10, missing_threshold_input = 10, 
          row_missing_p = 0.5, col_missing_p = 0.95) {
  rownames(RM_se) <- NULL
  indx_missing <- !(assays(RM_se)$IP >= missing_threshold_IP & 
                      assays(RM_se)$input >= missing_threshold_input)
  indx_col_sub <- colMeans(indx_missing) <= col_missing_p
  if (sum(indx_col_sub) < 10) {
    warning("The sample filtering is not applied because too little samples (< 10) are left afterward.")
    indx_col_sub = TRUE
  }
  indx_missing <- indx_missing[, indx_col_sub]
  indx_row_sub <- rowMeans(indx_missing) <= row_missing_p
  RM_se <- RM_se[indx_row_sub, indx_col_sub]
  indx_missing <- indx_missing[indx_row_sub, ]
  assays(RM_se)$NAindx <- indx_missing
  return(RM_se)
}

reduceSE = function (RM_se, gap_bp = 100) {
  site_gr <- rowRanges(RM_se)
  motif_reduced <- reduce(site_gr + gap_bp)
  fol <- findOverlaps(motif_reduced, site_gr)
  non_missing_num <- rowSums(!assays(RM_se)$NAindx)
  subj_nmn <- non_missing_num[subjectHits(fol)]
  names(subj_nmn) <- subjectHits(fol)
  hit_indx <- tapply(subj_nmn, queryHits(fol), function(x) names(which.max(x)))
  return(RM_se[as.numeric(hit_indx), ])
}

homoSE = filterSE(homoSE, missing_threshold_IP = 40, missing_threshold_input = 40, row_missing_p = 0.35)
homoSE = reduceSE(homoSE, gap_bp = 100)


IPCount <- assays(homoSE)[["IP"]]
inputCount <- assays(homoSE)[["input"]]

# count + 1
IPCountNorm = t(t(IPCount + 1)/colMeans(IPCount))
inputCountNorm = t(t(inputCount + 1)/colMeans(inputCount))
m6ALevel = log(IPCountNorm/inputCountNorm)

homeGrg = rowRanges(homoSE)
homeGrg$variance = rowVars(m6ALevel)


homeGrgData = as.data.frame(mcols(homeGrg))
homeGrgData = homeGrgData[order(homeGrgData$variance, decreasing = T),]
homeGrgData$m6ACali = ifelse(homeGrgData$p0 <= 0.5, 1, 0)

homeGrgData$cumFP = cumsum(homeGrgData$m6ACali == 0) / 1:nrow(homeGrgData)
cumFPDf = data.frame(cumFP = homeGrgData$cumFP,
                     class = ifelse(homeGrgData$m6ACali == 0, "FP", "TP"),
                     topSiteNum = 1:nrow(homeGrgData))

library(zoo)
library(scales)
bin = 50
bin_class <- rollapply(cumFPDf$class, width = bin, FUN = function(x) mean(x == "FP") > 27/bin, by = bin, align = "left")
tileDf <- data.frame(class = ifelse(bin_class, "FP", "TP"),
                     topSiteNum = seq(bin, nrow(homeGrgData), by = bin))

ggplot(data = cumFPDf, aes(x = topSiteNum, y = cumFP)) + 
  geom_path(size = 1, color = "steelblue") +
  geom_tile(data = tileDf, aes(x = topSiteNum, y = 0.42, fill = factor(class)), width = bin, color = NA, height = 0.02) +
  scale_fill_manual(values = c("FP" = "steelblue", "TP" = "#ffffff"), guide = FALSE) +
  theme_classic() +
  scale_y_continuous(limits = c(0.4, 0.7), 
                     breaks = seq(0.4, 0.7, by = 0.1),
                     labels = percent_format(accuracy = 1)) +
  scale_x_continuous(limits = c(380, 11000),
                     breaks = seq(0, 10000, by = 2000),
                     labels = comma, expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Number of top variance sites", 
       y = "Predicted off-targets by m6ACali") +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12))

ggsave("~/plots/metLevelVarTop.pdf", width = 3.6, height = 2.6)



################################################################################
####                           top variance GO                              ####
################################################################################

IPCount <- assays(homoSE)[["IP"]]
inputCount <- assays(homoSE)[["input"]]

# count + 1
IPCountNorm = t(t(IPCount + 1)/colMeans(IPCount))
inputCountNorm = t(t(inputCount + 1)/colMeans(inputCount))
m6ALevel = log(IPCountNorm/inputCountNorm)

homeGrg = rowRanges(homoSE)
homeGrg$variance = rowVars(m6ALevel)

library(clusterProfiler)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

exbg = exonsBy(txdb, by = "gene")
overlapGene = findOverlaps(homeGrg, exbg)
homeGrg$geneID[queryHits(overlapGene)] = names(exbg)[subjectHits(overlapGene)]

highConfIndx = which(homeGrg$p0 <= 0.5)
homeGrgPost = homeGrg[highConfIndx]

# top 2000 variance
preTopVarIndex = order(homeGrg$variance, decreasing = T)[1:2000]
postTopVarIndex = order(homeGrgPost$variance, decreasing = T)[1:2000]

ego <- enrichGO(gene = homeGrg$geneID[preTopVarIndex],
                OrgDb = org.Hs.eg.db, 
                keyType = 'ENTREZID', 
                ont = "ALL",
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                readable = TRUE)
ego2 <- enrichGO(gene = homeGrgPost$geneID[postTopVarIndex], 
                 OrgDb = org.Hs.eg.db, 
                 keyType = 'ENTREZID', 
                 ont = "ALL",
                 pvalueCutoff = 1,
                 qvalueCutoff = 1,
                 readable = TRUE)

ego_result <- as.data.frame(ego)
ego_result2 <- as.data.frame(ego2)

ego_post = ego_result2[order(ego_result2$p.adjust)[1:10],]
ego_pre = ego_result[ego_result$ID %in% ego_post$ID,]

ego_data = rbind(ego_pre, ego_post)

egoDf = data.frame(ID = ego_data$ID,
                   GeneRatio = apply(ego_data, 1, function(x) {
                     eval(parse(text = x["GeneRatio"]))
                   }),
                   padj = ego_data$p.adjust,
                   class = rep(c("pre", "post"), each = 10))
egoDf$class = factor(egoDf$class, levels = c("pre", "post"))

ggplot(data = egoDf, aes(x = class, y = ID, color = -log10(padj), size = GeneRatio)) + 
  geom_point() +
  scale_size(breaks = c(0.02, 0.04, 0.06)) +
  scale_color_gradientn(colours = c("#000ef9", "#d41f6b", "#f7251c"), name = expression(paste("-log"[10], "(p.adjust)"))) +
  scale_x_discrete(labels = c("Pre-Calibration", "Post-Calibration")) + 
  theme_bw() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_blank())

ggsave("~/plots/GO.pdf", width = 4.8, height = 3.7)


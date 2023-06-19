

########################################################################################################################
####                              Comparison between before and after calibration                                   ####
########################################################################################################################

library(caret)
library(predictiveFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

exbtx <- exonsBy(txdb, by = "tx")
motif_all <- sort(sampleSequence("DRACH", exbtx, bsgenome) - 2)

################################################################################
####                           Topology gradient                            ####
################################################################################

##############################################################
####                        before data                   ####
##############################################################

pred_df <- readRDS("~/rds/pred_df.rds")
mod_on_all_motif = readRDS("~/rds/mod_on_all_motif.rds")

merip = mod_on_all_motif[, 1:24]
single = mod_on_all_motif[, 25:49]
anti = mod_on_all_motif[, 1:49]

create_union_label <- function(union_df, one_cutoff) {
  num_ones <- rowSums(union_df == 1, na.rm = TRUE)
  num_zeros <- rowSums(union_df == 0, na.rm = TRUE)
  label <- as.factor(ifelse(num_ones > one_cutoff, 1, ifelse(num_zeros > 0, 0, NA)))
  return(label)
}

merip_label = create_union_label(merip, 10)
single_label = create_union_label(single, 0)
anti_label = create_union_label(data.frame(merip_label, single_label), 0)

merip_befo = motif_all[which(merip_label == 1)]
single_befo = motif_all[which(single_label == 1)]
anti_befo = motif_all[which(anti_label == 1)]

##############################################################
####                        after data                    ####
##############################################################

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

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

#region_weights = c(0.09078195, 0.45142217, 0.45779587)
region_weights = c(1/3, 1/3, 1/3)

cutoffs <- c(seq(0.2, 0.8, 0.2))

calculate_topology <- function(target, cutoff) {
  metagene_df = metagene(get(paste0(target, "_befo")), txdb, region_weights, paste0(target, "_befo"))
  
  for (cutoff in cutoffs) {
    motif_aft = motif_all[which(get(paste0(target, "_label")) == 1 & pred_df$mean_p0 < cutoff)]
    print(paste0(cutoff, ": ", length(motif_aft)))
    
    topology_aft = metagene(motif_aft, txdb, region_weights, paste0(target, "_", cutoff))
    
    metagene_df = rbind(metagene_df, topology_aft)
  }
  
  return(metagene_df)
}

metagene_anti = calculate_topology("anti", cutoffs)

##############################################################
####                       plot topology                  ####
##############################################################

library(ggplot2)
library(patchwork)

metagene_plot = ggplot() +
  geom_density(data = metagene_anti, aes(x = relative_pos, color = class)) +
  geom_vline(xintercept = c(region_weights[1], c(region_weights[1] + region_weights[2])), linetype = 5, color = "#707070") +
  scale_y_continuous(breaks = seq(0, 4, 0.5)) +
  ylab("Density") +
  scale_color_manual(values = c(colorRampPalette(c("#62061d", "#f0a787"))(4), "#2462ab"),
                     name = "Calibration Threshold",
                     labels = c("Highly Rigorous (n=78614)", "Rigorous (n=163025)", "Flexible (n=196986)",
                                "Highly Flexible (n=220391)", "Baseline (n=254945)")) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour= "black", size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour= "black", size = 12),
        legend.position = c(0.3, 0.75),
        legend.background = element_blank(),
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

ggsave("~/plots/metagene_gradient.pdf", width = 4.8, height = 4.1)


################################################################################
####                            Sensitivity gradient                        ####
################################################################################

##############################################################
####                     union label bar                  ####
##############################################################

library(caret)
pred_df <- readRDS("~/rds/pred_df.rds")
mod_on_all_motif = readRDS("~/rds/mod_on_all_motif.rds")

merip = mod_on_all_motif[, 1:24]
single = mod_on_all_motif[, 25:49]
anti = mod_on_all_motif[, 1:49]
free = mod_on_all_motif[, 50:64]

create_union_label_number <- function(union_df, one_cutoff) {
  num_ones <- rowSums(union_df == 1, na.rm = TRUE)
  num_zeros <- rowSums(union_df == 0, na.rm = TRUE)
  label <- ifelse(num_ones > one_cutoff, 1, ifelse(num_zeros > 0, 0, NA))
  return(label)
}

merip_label = create_union_label_number(merip, 10)
single_label = create_union_label_number(single, 0)
anti_label = create_union_label_number(data.frame(merip_label, single_label), 0)
free_label = create_union_label_number(free, 0)

new_mod = data.frame(merip_label = merip_label,
                     single_label = single_label,
                     anti_label = anti_label,
                     free_label = free_label)
new_mod[new_mod == 0] <- NA

calculate_percentage <- function(prediction, reference, cutoffs) {
  
  perc_df = data.frame()
  
  for (cutoff in cutoffs) {
    mod_aft = new_mod[which(pred_df$mean_p0 <= cutoff),]
    
    ref <- mod_aft[[paste0(reference, "_label")]]
    pred <- mod_aft[[paste0(prediction, "_label")]]
    
    ref <- ifelse(is.na(ref) & !is.na(pred), 0, ref)
    pred <- ifelse(is.na(pred) & !is.na(ref), 0, pred)
    
    cm <- confusionMatrix(reference = factor(ref, levels = c(1, 0)), data = factor(pred, levels = c(1, 0)))
    A = cm$table[1]
    C = cm$table[2]
    ref_1 = A+C
    A_perc <- A/(A+C)
    C_perc <- C/(A+C)
    
    print(cutoff)
    
    temp_df = data.frame(technique = prediction,
                         cutoff = cutoff,
                         perc = c(A_perc, C_perc),
                         class = c("A", "C"),
                         ref_1 = c(NA, paste0("n=", ref_1)))
    
    perc_df = rbind(perc_df, temp_df)
  }
  
  return(perc_df)
}

cutoffs <- c(0.1, 0.3, 0.5, 1)
perf_anti = calculate_percentage("anti", "free", cutoffs)

library(ggplot2)
perf_anti$class = factor(perf_anti$class, levels = c("C", "A"))

ggplot(perf_anti, aes(x = as.factor(cutoff), y = perc, fill = class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  geom_text(aes(label = ref_1, y = 1), vjust = -0.5, size = 3.2) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_x_discrete(breaks = cutoffs,
                   labels = c("Rigorous", "Balanced", "Flexible", "Baseline")) +
  scale_fill_manual(values = c("#1ebdc1", "#f2756d"), 
                    name = NULL,
                    breaks = c("C", "A"),
                    labels = c("No   ", "Yes")) +
  labs(x = "Calibration Threshold", y = expression(paste("Proportion of m"^6, "A sites")), title = "Overlap with Ground Truth") +
  theme_classic() +
  theme(axis.text = element_text(colour= "black", size = 10),
        axis.text.x = element_text(angle = 35, vjust = 0.9, hjust = 0.8),
        axis.title = element_text(colour= "black", size = 12),
        #plot.margin = margin(10, 16, 2, 0),
        legend.position = "top",
        legend.justification = "left",
        legend.box.just = "left",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(colour= "black", size = 12),
        #legend.margin = margin(0, 0, 0, -6),
        legend.box.margin = margin(-8, 0, -8, 0),
        plot.title = element_text(hjust = 0, vjust = 1.5, size = 12))

ggsave("~/plots/recallBar.pdf", width = 3.3, height = 4.3)



################################################################################
####                   FP read count & FP prob (gradient)                   ####
################################################################################

##############################################################
####                       prepare data                   ####
##############################################################

library(SummarizedExperiment)
library(ggplot2)
library(patchwork)

pred_df <- readRDS("~/rds/pred_df.rds")
se_human_SYSY <- readRDS("~/rds/se_human_SYSY.rds")

motif_all = rowRanges(se_human_SYSY) - 2
motif_all$count = rowMeans(assay(se_human_SYSY)[,se_human_SYSY$IP_input == "input"])
motif_all$p0 = pred_df$SYSY_p0


motif_FP = readRDS("./rds/peaks_IVT.rds")
overlaps <- findOverlaps(motif_FP, motif_all)

motif_FP$count[queryHits(overlaps)] <- motif_all$count[subjectHits(overlaps)]
motif_FP$p0[queryHits(overlaps)] <- motif_all$p0[subjectHits(overlaps)]


FP_info = as.data.frame(mcols(motif_FP))

# Using the cut function for grouping
breaks <- quantile(FP_info$count, probs = seq(0, 1, length.out = 11))
FP_info$count_group <- cut(FP_info$count, breaks = breaks, include.lowest = TRUE, labels = FALSE)

p_values <- numeric(length(unique(FP_info$count_group)) - 1)
min_values <- numeric(length(unique(FP_info$count_group)) - 1)

# Conduct wilcox tests for each pair of adjacent groups
for (i in 1:(length(unique(FP_info$count_group)) - 1)) {
  group1 <- FP_info[FP_info$count_group == i, "p0"]
  group2 <- FP_info[FP_info$count_group == (i + 1), "p0"]
  test <- wilcox.test(group1, group2, alternative = "less")
  p_values[i] <- test$p.value
  
  boxplot_values <- boxplot(group1, plot = FALSE)
  min_values[i] <- boxplot_values$stats[1]
}


##############################################################
####                          plot                        ####
##############################################################

boxes = ggplot(FP_info, aes(x = as.factor(count_group), y = p0, fill = as.factor(count_group))) +
  geom_boxplot(outlier.colour = NA, width = 0.55) +
  ylab("Predicted false positives probability") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0.05))) +
  stat_boxplot(geom = "errorbar",
               aes(ymin = ..ymax..),
               width = 0.2, size = .3) +
  stat_boxplot(geom = "errorbar",
               aes(ymax = ..ymin..),
               width = 0.2, size = .3) +
  scale_fill_manual(values = colorRampPalette(c("#d2e2ef", "#265998"))(10)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour= "black", size = 10),
        axis.title = element_text(colour= "black", size = 12),
        legend.position = "none",
        plot.margin = margin(0, 0, 1, 0))


add_p_value <- function(plot_object, x_start, x_end, y_position, p_value) {
  scientific_p <- sprintf("%.1e", p_value) 
  split_p <- strsplit(scientific_p, "e")[[1]]
  p_main <- split_p[1]
  p_exp <- sub("^-0", "-", split_p[2])
  exp <- bquote(.(p_main) ~ " x 10"^.(p_exp))
  
  plot_object <- plot_object + 
    annotate("segment", x = x_start, xend = x_start, y = y_position, yend = y_position + 0.01) +
    annotate("segment", x = x_end, xend = x_end, y = y_position, yend = y_position + 0.01) +
    annotate("segment", x = x_start, xend = x_end, y = y_position, yend = y_position) +
    annotate("text", x = (x_start + 0.95 + x_end) / 2, y = y_position - 0.085, 
             label = exp, size = 3.5, parse = F, angle = -40)
  return(plot_object)
}


for(i in 1:length(p_values)) {
  boxes = add_p_value(boxes, i, i + 1, min_values[i] - 0.02, p_values[i])
}

gradient = ggplot(data.frame(x = 0:100, y = 0.5), aes(x = x, y = y, fill = factor(x))) +
  geom_tile(color = NA, show.legend = F) +
  scale_fill_manual(values = colorRampPalette(c("#d2e2ef", "#265998"))(101)) +
  geom_polygon(data = data.frame(x = c(-0.5, 100.5, 100.5),
                                 y = c(0, 0, 1)),
               aes(x = x, y = y), colour = "black", fill = NA, size = 0.8) +
  xlab("Read coverage") +
  scale_x_continuous(breaks = c(-0.4,100.4), 
                     labels = c("1", "35475"), 
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_classic() +
  theme(axis.ticks.length = unit(0.4,"lines"),
        axis.text.x = element_text(colour = "black", size = 10),
        axis.ticks.x = element_line(colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(colour= "black", size = 12, vjust = 3),
        plot.margin = margin(-30, 0, -30, 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

layout <- c(
  area(t = 1, l = 0, b = 36, r = 50),
  area(t = 37, l = 3, b = 38, r = 48)
)
boxes + gradient + 
  plot_layout(design = layout)
ggsave("~/plots/coverage.pdf", width = 5.2, height = 3.9)



################################################################################
####                  Validation of predicted false positives               ####
################################################################################

##############################################################
####                        truth FP                      ####
##############################################################

motif_all <- sort(sampleSequence("DRACH", exbtx, bsgenome) - 2)
pred_df <- readRDS("~/rds/pred_df.rds")

motif_all$mean_p0 = pred_df$mean_p0
motif_all$mean_predict = as.numeric(as.character(pred_df$mean_predict))
motif_all_meta <- mcols(motif_all)
motif_all_meta$Abcam <- NA
motif_all_meta$NEB <- NA
motif_all_meta$SYSY <- NA

antibodies <- c("Abcam", "NEB", "SYSY")
for (ab in antibodies) {
  # load data
  peaks_IVT <- readRDS(paste0("~/peaks/", ab, "_IVT.rds"))
  peaks_mRNA <- readRDS(paste0("~/peaks/", ab, "_mRNA.rds"))
  
  TP = subsetByOverlaps(peaks_mRNA, peaks_IVT, invert = T)
  FP = peaks_IVT
  
  TP$class = 1  
  FP$class = 0
  
  AP = c(TP, FP)
  
  overlaps <- findOverlaps(motif_all, AP)
  
  motif_all_meta[queryHits(overlaps), ab] <- AP$class[subjectHits(overlaps)]
}
mcols(motif_all) <- motif_all_meta

meta_df = as.data.frame(motif_all_meta)

##############################################################
####                        indep FP                      ####
##############################################################

mod_on_all_motif = readRDS("~/rds/mod_on_all_motif.rds")

MeRIP = mod_on_all_motif[, 1:24]
miCLIP = mod_on_all_motif[, 25:39]
m6A_PA_CLIP = mod_on_all_motif[, 40:43]
m6A_CLIP = mod_on_all_motif[, 44:49]
single = mod_on_all_motif[, 25:49]
anti = mod_on_all_motif[, 1:49]

create_union_label <- function(union_df, one_cutoff) {
  num_ones <- rowSums(union_df == 1, na.rm = TRUE)
  num_zeros <- rowSums(union_df == 0, na.rm = TRUE)
  label <- as.factor(ifelse(num_ones > one_cutoff, 1, ifelse(num_zeros > 0, 0, NA)))
  return(label)
}

apply_label <- function(df, label_name, union_df, one_cutoff, p_cutoff) {
  label = create_union_label(union_df, one_cutoff)
  df[[label_name]] <- ifelse(label == 1 & pred_df$mean_p0 >= p_cutoff, 0, 
                             ifelse(label == 1 & pred_df$mean_p0 < p_cutoff, 1, NA))
  return(df)
}

meta_df <- apply_label(meta_df, "MeRIP", MeRIP, 1, 0.5)
meta_df <- apply_label(meta_df, "miCLIP", miCLIP, 1, 0.5)
meta_df <- apply_label(meta_df, "m6A_PA_CLIP", m6A_PA_CLIP, 1, 0.5)
meta_df <- apply_label(meta_df, "m6A_CLIP", m6A_CLIP, 1, 0.5)
meta_df <- apply_label(meta_df, "single", single, 0, 0.5)
meta_df <- apply_label(meta_df, "anti", anti, 0, 0.5)

##############################################################
####                          plot                        ####
##############################################################

support_percentage_df <- data.frame()
predict_colnames <- c("MeRIP", "miCLIP", "m6A_PA_CLIP", "m6A_CLIP", "single", "anti")

for (predict_colname in predict_colnames) {
  
  pred_rm_na = meta_df[which(meta_df[[predict_colname]] == 0), ]
  support_count <- rowSums(pred_rm_na[, c("Abcam", "NEB", "SYSY")] == pred_rm_na[[predict_colname]], na.rm = T)
  support_percentage <- table(support_count) / nrow(pred_rm_na)
  
  temp_df <- data.frame(Data = predict_colname,
                        SupportNumber = names(support_percentage),
                        Percentage = as.numeric(support_percentage))
  
  support_percentage_df <- rbind(support_percentage_df, temp_df)
  
  print(paste0(predict_colname, ": ", nrow(pred_rm_na)))
}

library(ggplot2)
support_percentage_df$SupportNumber = factor(support_percentage_df$SupportNumber, levels = c("3", "2", "1", "0"))
support_percentage_df$Data = factor(support_percentage_df$Data, levels = rev(c("MeRIP", "miCLIP", "m6A_PA_CLIP", "m6A_CLIP", "single", "anti")))

ggplot(support_percentage_df, aes(x = Data, y = Percentage, fill = SupportNumber)) +
  geom_bar(stat = "identity", position = "stack", color = "#445b83", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
  scale_x_discrete(breaks = c("MeRIP", "miCLIP", "m6A_PA_CLIP", "m6A_CLIP", "single", "anti"),
                   labels = c("MeRIP", "miCLIP", "m6A-PA-CLIP", "m6A-CLIP", "Single-based", "Anti-based")) +
  scale_fill_manual(values = c("#d2e2ef", "#90bfd5", "#488db8", "#265998"), 
                    name = NULL,
                    breaks = c("0", "1", "2", "3"),
                    labels = c("0 dataset        ", "1 dataset", "2 datasets         ", "3 datasets")) +
  labs(x = "", y = "Percentage of predicted false positives", title = "Validated by IVT") +
  coord_flip() +
  theme_classic() +
  guides(fill = guide_legend(nrow = 2, ncol = 2, byrow = TRUE)) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(colour= "black", size = 12),
        axis.text.y = element_text(margin = margin(l = -10)),
        axis.title = element_text(colour= "black", size = 12),
        axis.title.x = element_text(colour= "black", size = 12, vjust = 0),
        plot.margin = margin(10, 16, 5, 0),
        legend.position = "top",
        legend.justification = "left",
        legend.box.just = "left",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(colour= "black", size = 12),
        legend.margin = margin(0, 0, 0, -6),
        legend.box.margin = margin(0, 0, -10, 0),
        plot.title = element_text(hjust = 0, vjust = 1.5, size = 12))

ggsave("~/plots/overlapFP.pdf", width = 4.5, height = 3.4)



##############################################################
####                    plot merip cutoff                 ####
##############################################################

meta_df <- apply_label(meta_df, "MeRIP1", MeRIP, 1, 0.5)
meta_df <- apply_label(meta_df, "MeRIP5", MeRIP, 5, 0.5)
meta_df <- apply_label(meta_df, "MeRIP10", MeRIP, 10, 0.5)
meta_df <- apply_label(meta_df, "MeRIP15", MeRIP, 15, 0.5)
meta_df <- apply_label(meta_df, "MeRIP20", MeRIP, 20, 0.5)

support_percentage_merip_df <- data.frame()
predict_colnames <- c("MeRIP1", "MeRIP5", "MeRIP10", "MeRIP15", "MeRIP20")

for (predict_colname in predict_colnames) {
  
  pred_rm_na = meta_df[which(meta_df[[predict_colname]] == 0), ]
  support_count <- rowSums(pred_rm_na[, c("Abcam", "NEB", "SYSY")] == pred_rm_na[[predict_colname]], na.rm = T)
  support_percentage <- table(support_count) / nrow(pred_rm_na)
  
  temp_df <- data.frame(Data = predict_colname,
                        SupportNumber = names(support_percentage),
                        Percentage = as.numeric(support_percentage),
                        num = c(NA, NA, NA, paste0(nrow(pred_rm_na))))
  
  support_percentage_merip_df <- rbind(support_percentage_merip_df, temp_df)
  
  print(paste0(predict_colname, ": ", nrow(pred_rm_na)))
}

library(ggplot2)
support_percentage_merip_df$SupportNumber = factor(support_percentage_merip_df$SupportNumber, levels = c("0", "1", "2", "3"))
support_percentage_merip_df$Data = factor(support_percentage_merip_df$Data, levels = c("MeRIP1", "MeRIP5", "MeRIP10", "MeRIP15", "MeRIP20"))

ggplot(support_percentage_merip_df, aes(x = Data, y = Percentage, fill = SupportNumber)) +
  geom_bar(stat = "identity", position = "stack", color = "#445b83", width = 0.6) +
  geom_text(aes(label = num, y = 1), vjust = -0.5, size = 3.2) +
  scale_y_continuous(labels = scales::percent_format(), breaks = seq(0, 1, 0.2), limits = c(0, 1), 
                     expand = expansion(mult = c(0, 0.07)), position = "right") +
  scale_x_discrete(breaks = c("MeRIP1", "MeRIP5", "MeRIP10", "MeRIP15", "MeRIP20"),
                   labels = c(">0", ">5", ">10", ">15", ">20")) +
  scale_fill_manual(values = c("#d2e2ef", "#90bfd5", "#488db8", "#265998"), 
                    name = NULL,
                    breaks = c("0", "1", "2", "3"),
                    labels = c("0 dataset        ", "1 dataset", "2 datasets         ", "3 datasets")) +
  labs(x = expression(paste("Number of samples for m"^6, "A sites")), y = "Percentage of predicted false positives", title = "MeRIP (24 samples)") + #Sample Consistency Threshold
  theme_classic() +
  guides(fill = guide_legend(nrow = 2, ncol = 2, byrow = TRUE)) +
  theme(axis.text = element_text(colour= "black", size = 12),
        axis.title = element_text(colour= "black", size = 12),
        plot.margin = margin(10, 5, 5, 5),
        legend.position = "none",
        plot.title = element_text(hjust = 0, vjust = 1.5, size = 12))

ggsave("~/plots/overlapFPMeRIP.pdf", width = 3.2, height = 3.48)


################################################################################
####                                 Dot plot                               ####
################################################################################

library(caret)

pred_df <- readRDS("~/rds/pred_df.rds")
mod_on_all_motif <- readRDS("~/rds/mod_on_all_motif.rds")

mod_only_pos = mod_on_all_motif
mod_only_pos[mod_only_pos == 0] <- NA

cutoffs <- seq(0.1, 1, 0.1)

for (cutoff in cutoffs) {
  mod_aft = mod_only_pos[which(pred_df$mean_p0 <= cutoff),]
  
  F1_scores_aft <- matrix(0, nrow = ncol(mod_aft), ncol = ncol(mod_aft))
  precision_aft <- matrix(0, nrow = ncol(mod_aft), ncol = ncol(mod_aft))
  recall_aft <- matrix(0, nrow = ncol(mod_aft), ncol = ncol(mod_aft))
  accuracy_aft <- matrix(0, nrow = ncol(mod_aft), ncol = ncol(mod_aft))
  
  for(i in 1:ncol(mod_aft)) {
    for(j in 1:ncol(mod_aft)) {
      ref <- mod_aft[,i]
      pred <- mod_aft[,j]
      ref <- ifelse(is.na(ref) & !is.na(pred), 0, ref)
      pred <- ifelse(is.na(pred) & !is.na(ref), 0, pred)
      
      cm <- confusionMatrix(reference = factor(ref, levels = c(1, 0)), data = factor(pred, levels = c(1, 0)))
      F1_scores_aft[i,j] <- cm$byClass['F1']
      precision_aft[i,j] <- cm$byClass['Precision']
      recall_aft[i,j] <- cm$byClass['Recall']
      accuracy_aft[i,j] <- cm$overall['Accuracy']
    }
    print(i)
  }
  
  results = list(F1_scores_aft = F1_scores_aft, precision_aft = precision_aft, recall_aft = recall_aft, accuracy_aft = accuracy_aft)
  saveRDS(results, paste0("~/rds/scores_", cutoff, ".rds"))
}

##############################################################
####                  precision noly pos                  ####
##############################################################

library(ggplot2)
library(patchwork)
library(ggpointdensity)
scores_1 <- readRDS("~/rds/scores_1.rds")
scores_0.1 <- readRDS("~/rds/scores_0.1.rds")

precision_df = data.frame(Before = c(as.numeric(scores_1$recall_aft[c(50:64), c(1:24)])),
                          After = c(as.numeric(scores_0.1$recall_aft[c(50:64), c(1:24)])))

dot = ggplot(data = precision_df, aes(x = After, y = Before)) + 
  geom_pointdensity(size = 2.2, alpha = 0.7, show.legend = T) +
  scale_color_gradientn(colours = c("#4a95c5", "#2462ab", "#b31a2f", "#b31a2f"),
                        values = scales::rescale(c(0, 100, 200, 300))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
  xlab("Post-Calibration Sensitivity") +
  ylab("Pre-Calibration Sensitivity") +
  labs(color = "Sample Pairs Density") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1),expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
  theme_classic() + 
  theme(axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12),
        legend.position = c(0.32, 0.75),
        legend.title = element_text(colour= "black", size = 12),
        legend.text = element_text(colour= "black", size = 12),
        legend.background = element_blank(),
        plot.margin = margin(5, 10, 5, 5))

# histogram + density
precision_befo = data.frame(precision = as.numeric(scores_1$recall_aft[c(50:64), c(1:24)]))
precision_aft = data.frame(precision = as.numeric(scores_0.1$recall_aft[c(50:64), c(1:24)]))

density_befo = ggplot(data = precision_befo, aes(x = precision)) +
  geom_histogram(bins = 40, alpha = 0.7, color = "white", fill = "#96c4db") +
  geom_density(aes(y = ..density.. * 9.2), alpha = 0.1, color = "#2164ab") +
  scale_y_continuous(breaks = c(0, 40), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  theme_classic() + 
  theme(axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour= "black", size = 12),
        axis.text.y = element_blank(),
        axis.title = element_blank()) +
  coord_flip()

density_aft = ggplot(data = precision_aft, aes(x = precision)) +
  geom_histogram(bins = 40, alpha = 0.7, color = "white", fill = "#d66255") +
  geom_density(aes(y = ..density.. * 9.2), alpha = 0.1, color = "#62061d") +
  scale_y_continuous(breaks = c(0, 25), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  theme_classic() + 
  theme(axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour= "black", size = 12),
        axis.text.x = element_blank(),
        axis.title = element_blank())

layout <- c(
  area(t = 0, l = 0, b = 15, r = 50),
  area(t = 16, l = 0, b = 100, r = 50),
  area(t = 16, l = 51, b = 100, r = 60)
)
density_aft + dot + density_befo + plot_layout(design = layout)

ggsave("~/plots/recallDot.pdf", width = 4.8, height = 4.3)



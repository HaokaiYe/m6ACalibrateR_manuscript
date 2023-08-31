
########################################################################################################################
####                                   Thorough analysis across DRACH motifs                                        ####
########################################################################################################################

################################################################################
####                             color bar code                             ####
################################################################################

##############################################################
####                       prepare data                   ####
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

create_union_label_number <- function(union_df) {
  num_ones <- rowSums(union_df == 1, na.rm = TRUE)
  num_zeros <- rowSums(union_df == 0, na.rm = TRUE)
  label <- ifelse(num_ones > 0, 1, ifelse(num_zeros > 0, 0, NA))
  return(label)
}

meta_df$union = create_union_label_number(meta_df[, 3:5])


##############################################################
####                          plot                        ####
##############################################################

library(ggplot2)
library(dplyr)
meta_df_sorted <- meta_df[order(meta_df$mean_p0), ]


meta_df_sorted$union_TP = ifelse(meta_df_sorted$mean_predict == 1 & meta_df_sorted$union == 1, TRUE, FALSE)
meta_df_sorted$union_TP[is.na(meta_df_sorted$union_TP)] <- FALSE
meta_df_sorted$union_FP = ifelse(meta_df_sorted$mean_predict == 0 & meta_df_sorted$union == 0, TRUE, FALSE)
meta_df_sorted$union_FP[is.na(meta_df_sorted$union_FP)] <- FALSE

library(zoo)
bin = 10
bin_TP <- rollapply(meta_df_sorted$union_TP, width = bin, FUN = function(x) mean(x) > 1/bin, by = bin, align = "left")
bin_FP <- rollapply(meta_df_sorted$union_FP, width = bin, FUN = function(x) mean(x) > 1/bin, by = bin, align = "left")


bin_df = data.frame(TP = bin_TP,
                    FP = bin_FP,
                    color = ifelse(bin_FP == TRUE, "FP", ifelse(bin_TP == TRUE, "TP", "None")),
                    index = 1:(nrow(meta_df_sorted)/bin))


ggplot(bin_df, aes(x = index, y = 1, fill = color)) +
  geom_tile(color = NA) +
  scale_fill_manual(values = c("FP" = "#477ca8", "TP" = "#cb3335", "None" = "white"), guide = FALSE) +
  geom_polygon(data = data.frame(x = c(min(bin_df$index), max(bin_df$index), max(bin_df$index)),
                                 y = c(0, 0, 0.4)),
               aes(x = x, y = y), fill = "#477ca8") +
  annotate("text", x = max(bin_df$index)/2, y = -0.15, 
           label = "Prob of false positive (FP) predicted by m6ACali", 
           size = 3.6, colour = "black", hjust = 0.5) +
  xlab(expression(paste("Ranked high confidence / false positive m"^6, "A sites"))) +
  scale_x_continuous(breaks = c(min(bin_df$index), max(bin_df$index)), 
                     labels = c("2903280", "1"), 
                     position = "top",
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0))) +
  theme(
    panel.background = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 10),
    axis.text.x.top = element_text(colour = "black", size = 9),
    axis.ticks.x.top = element_line(size = 0.5)
  )

ggsave("~/plots/code.pdf", width = 8, height = 1.5)


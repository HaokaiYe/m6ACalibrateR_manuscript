
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

##############################################################
####                        anti free                     ####
##############################################################

free <- readRDS("./rds/m6A_Atlas2_high_resolution_hg38.rds")
free = free[free$Technique_Num >= 2]


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
  xlab("Total width of peaks called (million bases)") +
  ylab("Overlapped high-resolution sites") +
  scale_x_continuous(limits = c(0, 2)) +
  scale_y_continuous(limits = c(250, 32768), trans = log10_trans()) +
  scale_color_manual(values = c("#0076ba", "#cb3335"),
                     name = NULL,
                     breaks = c("Before", "After"),
                     labels = c("Pre-Calibration", "Post-Calibration")) +
  theme_classic() + 
  theme(axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12),
        legend.title = element_blank(),
        legend.position = c(0.87, 0.2),
        legend.background = element_blank(),
        legend.text = element_text(colour = "black", size = 10),
        strip.text = element_text(colour = "black", size = 10),
        strip.background = element_rect(size = 0.8, color = "black")) +
  facet_wrap(facet = . ~ method, scales = "free")

ggsave("./plots/callmethods.pdf", width = 7.3, height = 3.4)


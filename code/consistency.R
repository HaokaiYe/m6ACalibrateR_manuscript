
########################################################################################################################
####                                technical independent validation via odds ratios                                ####
########################################################################################################################

library(SummarizedExperiment)
library(h2o)

library(rtracklayer)
library(AnnotationHub)

library(predictiveFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)


################################################################################
####                         get hg38 merip path                            ####
################################################################################
merip_hg38 <- list.dirs("~/merip_data_hg38")
merip_hg38_exomePeak <- grep("exomePeak2_output", merip_hg38, value=T)
merip_hg38_peaks <- paste0(merip_hg38_exomePeak, "/peaks.csv")


################################################################################
####                         construct matrix                               ####
################################################################################
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

exbtx <- exonsBy(txdb, by = "tx")
motif_all <- sampleSequence("DRACH", exbtx, bsgenome)


modification_merip <- as.data.frame(matrix(data = NA, nrow = length(motif_all), ncol = length(merip_hg38_peaks)))
for (i in 1:length(merip_hg38_peaks)) {
  peaks <- read.csv(merip_hg38_peaks[i])
  
  peaks_grg <- GRanges(seqnames = peaks$chr, 
                             ranges = IRanges(peaks$chromStart, peaks$chromEnd),
                             strand = peaks$strand)
  
  sample_exbtx <- subsetByOverlaps(exbtx, peaks_grg)
  
  # Overlap exons containing m6a in transcripts with all DRACH motifs and observe all m6a motifs on the target transcript
  motif_all_TP_index <- unique(as.data.frame(findOverlaps(motif_all, peaks_grg))$queryHits)
  motif_all_TN_index <- unique(as.data.frame(findOverlaps(motif_all, sample_exbtx))$queryHits)
  
  # Set all m6a motifs on the target transcript to 0 (TN), and then use the original 1 to cover them (TP)
  modification_merip[motif_all_TN_index, i] <- 0
  modification_merip[motif_all_TP_index, i] <- 1
}  


################################################################################
####                         import data from Zhen                          ####
################################################################################

single_based_m6A <- readRDS("~/single_based_m6A.rds")
sample <- as.data.frame(colData(single_based_m6A))
motif_hg19 <- rowRanges(single_based_m6A)
modification <- as.data.frame(assay(single_based_m6A))


################################################################################
####                         liftOver hg19 to hg38                          ####
################################################################################

chain = import.chain("~/hg19ToHg38.over.chain")

motif_hg38 <- liftOver(motif_hg19, chain)

# Obtain the coordinates that failed to liftOver
discard_index <- which(is.na(as.numeric(start(motif_hg38))))

# Remove the rows with failed liftOver coordinates
motif_hg38_discarded <- unlist(motif_hg38[-discard_index])
modification_discarded <- modification[-discard_index,]


################################################################################
####                        generate false positives                        ####
################################################################################

# Remove positive positions that are no longer DRACH motifs after liftOver
remove_index <- as.data.frame(findOverlaps(motif_hg38_discarded, motif_all))$queryHits
motif_hg38_removed <- motif_hg38_discarded[remove_index]
modification_removed <- modification_discarded[remove_index,]

# Define a new matrix
modification_new <- as.data.frame(matrix(data = NA, nrow = length(motif_all), ncol = 40))

for (i in 1:ncol(modification_new)) {
  # Obtain the coordinates of m6a TP detected in the samples
  sample1_TP_index <- which(modification_removed[, i] == 1)
  sample1_TP <- motif_hg38_removed[sample1_TP_index]
  
  # Map the coordinates of m6a TP detected in the samples to the transcript, and obtain exons containing m6a in the transcript with m6a TP
  sample1_exbtx <- subsetByOverlaps(exbtx, sample1_TP)
  
  # Overlap exons containing m6a in transcripts with all DRACH motifs and observe all m6a motifs on the target transcript
  motif_all_TP_index <- unique(as.data.frame(findOverlaps(motif_all, sample1_TP))$queryHits)
  motif_all_TN_index <- unique(as.data.frame(findOverlaps(motif_all, sample1_exbtx))$queryHits)
  
  # Set all m6a motifs on the target transcript to 0 (TN), and then use the original 1 to cover them (TP)
  modification_new[motif_all_TN_index, i] <- 0
  modification_new[motif_all_TP_index, i] <- 1
}


# Modify the DRACA motif in mazter-seq and change the 0 that is not DRACA to NA
mazter_index <- which(sample$Technique == "MAZTER-Seq")

for (m in mazter_index) {
  mazter_0_index = which(modification_new[, m] == 0)
  mazter_0_motif = motif_all[mazter_0_index]
  
  mazter_0_motif_seq = DNAStringSet(Views(bsgenome, mazter_0_motif))
  change_index = mazter_0_index[which(vcountPattern("DRACA", mazter_0_motif_seq, fixed = FALSE) == 0)]
  modification_new[change_index, m] = NA
}


# Merge samples
modification_merip_new <- cbind(modification_merip, modification_new[, c(1:4)], modification_new[, c(20:40)], modification_new[, c(5:19)])

# Remove rows with all samples that are all NA
row_na_sum <- apply(modification_merip_new, 1, function(x){sum(is.na(x))})
na_index <- which(row_na_sum == ncol(modification_merip_new))

modification_na_rm <- modification_merip_new[-na_index,]
motif_all_na_rm <- motif_all[-na_index]


################################################################################
####                         calculate odds ratios                          ####
################################################################################

# Define a new matrix to store odds ratio
odds_ratio_table <- matrix(data = NA, nrow = ncol(modification_merip_new), ncol = ncol(modification_merip_new))

for (i in 1:ncol(modification_merip_new)) {
  for (j in 1:ncol(modification_merip_new)) {
    odds_ratio_0_0 <- which((modification_na_rm[,i] + modification_na_rm[,j]) == 0)
    odds_ratio_1_1 <- which((modification_na_rm[,i] + modification_na_rm[,j]) == 2)
    odds_ratio_0_1 <- which((modification_na_rm[,i] - modification_na_rm[,j]) == -1)
    odds_ratio_1_0 <- which((modification_na_rm[,i] - modification_na_rm[,j]) == 1)
    
    odds_ratio_matrix <- matrix(c(length(odds_ratio_1_1), length(odds_ratio_1_0), length(odds_ratio_0_1), length(odds_ratio_0_0)), nrow = 2, ncol = 2)
    
    odds_ratio_ftest <- fisher.test(odds_ratio_matrix)
    odds_ratio_table[i, j] = odds_ratio_ftest$estimate
  }
}


################################################################################
####                         extract genome feature                         ####
################################################################################

### Extract features (only genome features)
gfeatures <- genomeDerivedFeatures(x = motif_all_na_rm,
                                   transcriptdb = txdb)


################################################################################
####                         perform h2o rf models                          ####
################################################################################

### h2o model
h2o.init()

h2o_gfeatures <- as.h2o(gfeatures)

h2o_rf_Abcam <- h2o.loadModel("~/h2o_rf_Abcam")
h2o_rf_Abcam_pred <- as.data.frame(h2o.predict(h2o_rf_Abcam, newdata = h2o_gfeatures))


h2o_rf_NEB <- h2o.loadModel("~/h2o_rf_NEB")
h2o_rf_NEB_pred <- as.data.frame(h2o.predict(h2o_rf_NEB, newdata = h2o_gfeatures))


h2o_rf_SYSY <- h2o.loadModel("~/h2o_rf_SYSY")
h2o_rf_SYSY_pred <- as.data.frame(h2o.predict(h2o_rf_SYSY, newdata = h2o_gfeatures))


pred = data.frame(abcam = h2o_rf_Abcam_pred$predict, neb = h2o_rf_NEB_pred$predict, sysy = h2o_rf_SYSY_pred$predict)
predp0 = data.frame(abcam = h2o_rf_Abcam_pred$p0, neb = h2o_rf_NEB_pred$p0, sysy = h2o_rf_SYSY_pred$p0)
pred$mean = rowMeans(predp0) < 0.5


filter_index <- which(pred$mean == FALSE)

modification_FP_rm <- modification_na_rm[-filter_index,]


################################################################################
####                      calculate odds ratios again                       ####
################################################################################

# Define a new matrix to store odds ratio
odds_ratio_table_FP_rm <- matrix(data = NA, nrow = ncol(modification_merip_new), ncol = ncol(modification_merip_new))

for (i in 1:ncol(modification_merip_new)) {
  for (j in 1:ncol(modification_merip_new)) {
    odds_ratio_0_0 <- which((modification_FP_rm[,i] + modification_FP_rm[,j]) == 0)
    odds_ratio_1_1 <- which((modification_FP_rm[,i] + modification_FP_rm[,j]) == 2)
    odds_ratio_0_1 <- which((modification_FP_rm[,i] - modification_FP_rm[,j]) == -1)
    odds_ratio_1_0 <- which((modification_FP_rm[,i] - modification_FP_rm[,j]) == 1)
    
    odds_ratio_matrix <- matrix(c(length(odds_ratio_1_1), length(odds_ratio_1_0), length(odds_ratio_0_1), length(odds_ratio_0_0)), nrow = 2, ncol = 2)
    
    odds_ratio_ftest <- fisher.test(odds_ratio_matrix)
    odds_ratio_table_FP_rm[i, j] = odds_ratio_ftest$estimate
  }
}


################################################################################
####                          plot boxplot  odds ratio                      ####
################################################################################

# antibody merip: 1:24
# antibody single base: 25:49
# antibody-free: 50:64

# before filter
merip_vs_anti_free_before <- c(as.numeric(odds_ratio_table[c(1:24), c(50:64)]))
single_vs_anti_free_before <- c(as.numeric(odds_ratio_table[c(25:49), c(50:64)]))
anti_vs_anti_free_before <- c(as.numeric(odds_ratio_table[c(1:49), c(50:64)]))

odds_ratio_before <- data.frame(odds_ratio = c(merip_vs_anti_free_before, single_vs_anti_free_before, anti_vs_anti_free_before),
                                class = c(rep("merip vs anti-free", length(merip_vs_anti_free_before)),
                                          rep("single vs anti-free", length(single_vs_anti_free_before)),
                                          rep("anti vs anti-free", length(anti_vs_anti_free_before))))


# after filter
merip_vs_anti_free_after <- c(as.numeric(odds_ratio_table_FP_rm[c(1:24), c(50:64)]))
single_vs_anti_free_after <- c(as.numeric(odds_ratio_table_FP_rm[c(25:49), c(50:64)]))
anti_vs_anti_free_after <- c(as.numeric(odds_ratio_table_FP_rm[c(1:49), c(50:64)]))


odds_ratio_after <- data.frame(odds_ratio = c(merip_vs_anti_free_after, single_vs_anti_free_after, anti_vs_anti_free_after),
                               class = c(rep("merip vs anti-free", length(merip_vs_anti_free_after)),
                                         rep("single vs anti-free", length(single_vs_anti_free_after)),
                                         rep("anti vs anti-free", length(anti_vs_anti_free_after))))


# combine before and after
odds_ratio_data <- rbind(odds_ratio_before, odds_ratio_after)
odds_ratio_data$group <- c(rep("Before", nrow(odds_ratio_before)), rep("After", nrow(odds_ratio_after)))

saveRDS(odds_ratio_data, "./rds/tech_indep_valid.rds")



################################################################################
####                                   plot                                 ####
################################################################################

library(ggplot2)

odds_ratio_data <- readRDS("./rds/tech_indep_valid.rds")
odds_ratio_data$group <- factor(odds_ratio_data$group, levels = c("Before", "After"))


merip_anti_free_before <- which(odds_ratio_data$class == "merip vs anti-free" & odds_ratio_data$group == "Before") # 1:360
merip_anti_free_after <- which(odds_ratio_data$class == "merip vs anti-free" & odds_ratio_data$group == "After") # 1471:1830
set.seed(123)
merip_anti_free_index <- sample(1:360, 100)
merip_anti_free <- c(merip_anti_free_before[merip_anti_free_index], merip_anti_free_after[merip_anti_free_index])


single_anti_free_before <- which(odds_ratio_data$class == "single vs anti-free" & odds_ratio_data$group == "Before") # 361:735
single_anti_free_after <- which(odds_ratio_data$class == "single vs anti-free" & odds_ratio_data$group == "After") # 1831:2205
set.seed(123)
single_anti_free_index <- sample(1:375, 100)
single_anti_free <- c(single_anti_free_before[single_anti_free_index], single_anti_free_after[single_anti_free_index])


anti_anti_free_before <- which(odds_ratio_data$class == "anti vs anti-free" & odds_ratio_data$group == "Before") # 736:1470
anti_anti_free_after <- which(odds_ratio_data$class == "anti vs anti-free" & odds_ratio_data$group == "After") # 2206:2940
set.seed(123)
anti_anti_free_index <- sample(1:735, 100)
anti_anti_free <- c(anti_anti_free_before[anti_anti_free_index], anti_anti_free_after[anti_anti_free_index])


sub_index <- unique(c(merip_anti_free, single_anti_free, anti_anti_free))
sub_data <- odds_ratio_data[sub_index,]


ggplot(odds_ratio_data, aes(x = group, y = odds_ratio, color = group)) + 
  geom_boxplot(outlier.colour = NA, size = 0.75) +
  geom_jitter(data = sub_data, alpha = 0.5, size = 0.8, width = 0.2) +
  scale_color_manual(values = c("#0076ba", "#f27200")) +
  ylab("Consistency (Odds ratio)") +
  theme_bw() +
  theme(panel.border = element_rect(size = 0.8),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.text.y = element_text(colour= "black", size = 10),
        axis.text.x = element_text(colour = "black", size = 10),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour= "black", size = 12),
        strip.text = element_text(colour = "black", size = 10),
        strip.background = element_rect(size = 0.8),
        panel.spacing = unit(0.5, "cm")) +
  facet_wrap(facet = . ~ class, scales = "free_y")

ggsave("~/plots/odds_ratio.pdf", width = 8, height = 2.6)





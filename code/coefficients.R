
########################################################################################################################
####                          Consistency of sequence features between true & false positive                        ####
########################################################################################################################

library(predictiveFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

################################################################################
####            Using IVT and mRNA data to determine two groups:            ####
####                    true positive and false positive                    ####
################################################################################

peaks_mRNA <- readRDS("./rds/peaks_mRNA.rds")
peaks_IVT <- readRDS("./rds/peaks_IVT.rds")

# Modification sites identified exclusively in the mRNA sample were considered true positives,
# while all sites identified in the IVT sample were deemed false positives
TP = peaks_mRNA[-queryHits(findOverlaps(peaks_mRNA, peaks_IVT))]
FP = peaks_IVT

################################################################################
####  Select negative sites that are equal in quantity and non-overlapping  ####
################################################################################

# Extract all DRACH motifs
exbtx <- exonsBy(txdb, by = "tx")
motif_all <- sort(sampleSequence("DRACH", exbtx, bsgenome) - 2)

position_coef <- function(positive_grg) {

  # Select exons containing m6a
  exbtx_contain_m6a <- subsetByOverlaps(exbtx, positive_grg)

  # Define positive and negative (1bp center of DRACH)
  positive = subsetByOverlaps(motif_all, positive_grg)

  neg = subsetByOverlaps(motif_all, exbtx_contain_m6a)
  negative_all = neg[-unique(queryHits(findOverlaps(neg, positive)))]


  # Ensure that the number of negative and positive samples is equal and the central A is non-overlapping
  set.seed(110)
  ind_negative = sample(1:length(negative_all), length(positive))
  negative = negative_all[ind_negative]

################################################################################
####                 Build LR model to calculate coefficients               ####
################################################################################

  # Obtain an 11-bp one-hot encoding
  seq = c(positive + 5, negative + 5)
  seq_onehot <- sequenceDerivedFeatures(seq,
                                        bsgenome,
                                        encoding = "onehot")

  # "1" is "positive", "0" is "negative"
  seq_onehot$class <- factor(c(rep(1, length(positive)), rep(0, length(negative))))


  lr = glm(class ~ ., data = seq_onehot, family = binomial(link = "logit"))
  coef = coefficients(lr)[-1]

  ################################################################################
  ####                          Scale the coefficients                        ####
  ################################################################################

  # Exclude NA values from scaling
  coef_split = split(coef, factor(rep(-5:5, each = 4)))
  coef_scaled = unname(unlist(lapply(coef_split, function(x) {x - mean(x, na.rm = T)})))

  coef_scaled_0_1 = unname(unlist(lapply(coef_split, function(x) {(x - mean(x, na.rm = T))/sd(x, na.rm = T)})))

  # Include NA values in scaling
  coef2 = coef
  coef2[is.na(coef)] = 0
  coef_NA_split = split(coef2, factor(rep(-5:5, each = 4)))
  coef_NA_scaled = c()
  for (i in seq_along(coef_NA_split)) {
    four_coef = unname(coef_NA_split[[i]])
    if(i == 4) {
      four_coef_scaled = four_coef[-3] - mean(four_coef[-3])
      coef_NA_scaled = c(coef_NA_scaled, append(four_coef_scaled, NA, 2))
    } else if(i == 5) {
      four_coef_scaled = four_coef[-c(2,3)] - mean(four_coef[-c(2,3)])
      coef_NA_scaled = c(coef_NA_scaled, append(four_coef_scaled, c(NA, NA), 1))
    } else if(i == 6 || i == 7) {
      coef_NA_scaled = c(coef_NA_scaled, rep(NA, 4))
    } else if(i == 8) {
      four_coef_scaled = four_coef[-4] - mean(four_coef[-4])
      coef_NA_scaled = c(coef_NA_scaled, append(four_coef_scaled, NA, 3))
    } else {
      four_coef_scaled = four_coef - mean(four_coef)
      coef_NA_scaled = c(coef_NA_scaled, four_coef_scaled)
    }
  }

  coef_NA_scaled_0_1 = c()
  for (i in seq_along(coef_NA_split)) {
    four_coef = unname(coef_NA_split[[i]])
    if(i == 4) {
      four_coef_scaled = (four_coef[-3] - mean(four_coef[-3]))/sd(four_coef[-3])
      coef_NA_scaled_0_1 = c(coef_NA_scaled_0_1, append(four_coef_scaled, NA, 2))
    } else if(i == 5) {
      four_coef_scaled = (four_coef[-c(2,3)] - mean(four_coef[-c(2,3)]))/sd(four_coef[-c(2,3)])
      coef_NA_scaled_0_1 = c(coef_NA_scaled_0_1, append(four_coef_scaled, c(NA, NA), 1))
    } else if(i == 6 || i == 7) {
      coef_NA_scaled_0_1 = c(coef_NA_scaled_0_1, rep(NA, 4))
    } else if(i == 8) {
      four_coef_scaled = (four_coef[-4] - mean(four_coef[-4]))/sd(four_coef[-4])
      coef_NA_scaled_0_1 = c(coef_NA_scaled_0_1, append(four_coef_scaled, NA, 3))
    } else {
      four_coef_scaled = (four_coef - mean(four_coef))/sd(four_coef)
      coef_NA_scaled_0_1 = c(coef_NA_scaled_0_1, four_coef_scaled)
    }
  }

  # Normalize all ATCG without considering DRACH
  coef_NA_scaled_all = unname(unlist(lapply(coef_NA_split, function(x) {x - mean(x)})))
  coef_NA_scaled_all_0_1 = unname(unlist(lapply(coef_NA_split, function(x) {(x - mean(x))/sd(x)})))


  result = list(coef = coef, coef_scaled = coef_scaled, coef_scaled_0_1 = coef_scaled_0_1,
                coef_NA_scaled = coef_NA_scaled, coef_NA_scaled_0_1 = coef_NA_scaled_0_1,
                coef_NA_scaled_all = coef_NA_scaled_all, coef_NA_scaled_all_0_1 = coef_NA_scaled_all_0_1)
  return(result)
}

################################################################################
####                                 plot                                   ####
################################################################################
library(ggplot2)

### TP

TP_coef = position_coef(TP)

coef_df = data.frame(coef = TP_coef$coef_NA_scaled_0_1,
                     position = factor(rep(-5:5, each = 4)),
                     base = factor(rep(c("A", "T", "C", "G"), 11), levels = c("A", "T", "C", "G")))

coef_df$coef[coef_df$coef == 0] = NA
coef_df$coef[coef_df$position == 0] = 0
coef_df$base[coef_df$position == 0] = "A"
coef_df$coef[coef_df$position == 1] = 0
coef_df$base[coef_df$position == 1] = "C"

ggplot(coef_df, aes(x = position, y = coef, color = base)) +
  #geom_point() +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_text(aes(label = base), size = 5, fontface = "bold") +
  xlab("Positions") +
  ylab("Coefficient") +
  ggtitle(expression(paste("High confidence m"^6, "A"))) +
  scale_color_manual(values = c("#e22127",
                                "#98519f",
                                "#337fba",
                                "#4caf47")) +
  theme_classic() +
  theme(axis.text = element_text(colour= "black", size = 9),
        legend.position = 'none',
        axis.title = element_text(colour= "black", size = 12),
        plot.title = element_text(colour = "black", size = 12, hjust = 0.5))

ggsave("~/plots/seq_SYSY_TP.pdf", width = 3.8, height = 3.0)

### FP

FP_coef = position_coef(FP)

FP_coef_df = data.frame(coef = FP_coef$coef_NA_scaled_0_1,
                     position = factor(rep(-5:5, each = 4)),
                     base = factor(rep(c("A", "T", "C", "G"), 11), levels = c("A", "T", "C", "G")))

FP_coef_df$coef[FP_coef_df$coef == 0] = NA
FP_coef_df$coef[FP_coef_df$position == 0] = 0
FP_coef_df$base[FP_coef_df$position == 0] = "A"
FP_coef_df$coef[FP_coef_df$position == 1] = 0
FP_coef_df$base[FP_coef_df$position == 1] = "C"

ggplot(FP_coef_df, aes(x = position, y = coef, color = base)) +
  #geom_point() +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_text(aes(label = base), size = 5, fontface = "bold") +
  xlab("Positions") +
  ylab("Coefficient") +
  ggtitle(expression(paste("False positive m"^6, "A"))) +
  scale_color_manual(values = c("#e22127",
                                "#98519f",
                                "#337fba",
                                "#4caf47")) +
  theme_classic() +
  theme(axis.text = element_text(colour= "black", size = 9),
        legend.position = 'none',
        axis.title = element_text(colour= "black", size = 12),
        plot.title = element_text(colour = "black", size = 12, hjust = 0.5))

ggsave("~/plots/seq_SYSY_FP.pdf", width = 3.8, height = 3.0)


################################################################################
####                           Coefficient plot                             ####
################################################################################
cor.test(TP_coef$coef_NA_scaled_0_1, FP_coef$coef_NA_scaled_0_1)
cor.test(TP_coef$coef_NA_scaled_0_1, FP_coef$coef_NA_scaled_0_1, method = "spearman")

all_coef_df = data.frame(TP_coef = TP_coef$coef, TP_coef_NA_scaled = TP_coef$coef_NA_scaled, TP_coef_NA_scaled_0_1 = TP_coef$coef_NA_scaled_0_1,
                         FP_coef = FP_coef$coef, FP_coef_NA_scaled = FP_coef$coef_NA_scaled, FP_coef_NA_scaled_0_1 = FP_coef$coef_NA_scaled_0_1)

ggplot(all_coef_df, aes(x = FP_coef_NA_scaled_0_1, y = TP_coef_NA_scaled_0_1)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm",size = 1.3, color = "#4a68b1", se = FALSE) +
  ylab("Coefficient TP") +
  xlab("Coefficient FP") +
  #ggtitle("Correlation of coefficient TP and FP") +
  annotate('text', x = -1.48, y = 1, label = "R = 0.7472\np = 8.991e-07\nrho = 0.7947\np = 8.302e-07", hjust = 0, size = 3.6) +
  theme_classic() +
  theme(axis.title = element_text(colour= "black", size = 12),
        plot.title = element_text(colour = "black", size = 12, hjust = 0.5))

ggsave("~/plots/seq_SYSY_corr.pdf", width = 3.5, height = 3.0)


################################################################################
####                          Contingency table                             ####
################################################################################

coef_rank = function(coef) {

  coef_split = split(coef, factor(rep(-5:5, each = 4)))

  coef_ranks = c()
  for (i in seq_along(coef_split)) {
    four_coef = coef_split[[i]]

    if(i < 4 || i > 8) {
      coef_rank = paste0(c("A", "T", "C", "G"), 5-rank(four_coef))
      coef_ranks = c(coef_ranks, coef_rank)
    }
  }

  return(coef_ranks)
}

TP_rank = coef_rank(TP_coef$coef_NA_scaled)
FP_rank = coef_rank(FP_coef$coef_NA_scaled)

rank_table = table(TP_rank == FP_rank)


coef_rank_matrix = matrix(c(rank_table["TRUE"], rank_table["FALSE"], 1/4, 3/4), nrow = 2, byrow = T,
                          dimnames = list(c("Observed", "Expected"), c("same", "diff")))

chisq.test(x = c(rank_table["TRUE"], rank_table["FALSE"]), p = c(1/4, 3/4))





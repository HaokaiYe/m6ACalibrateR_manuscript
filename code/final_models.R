
########################################################################################################################
####                                               Build up the final models                                        ####
########################################################################################################################

library(pROC)
library(h2o)
library(randomForest)

library(predictiveFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

################################################################################
####                     Load mRNA peaks and IVT peaks                      ####
################################################################################

peaks_mRNA <- readRDS("./rds/peaks_mRNA.rds")
peaks_IVT <- readRDS("./rds/peaks_IVT.rds")

################################################################################
####           Define true positives (TP) and false positives (FP)          ####
################################################################################

# Modification sites identified exclusively in the mRNA sample were considered true positives,
# while all sites identified in the IVT sample were deemed false positives
TP = peaks_mRNA[-queryHits(findOverlaps(peaks_mRNA, peaks_IVT))]
FP = peaks_IVT

################################################################################
####                            Extract features                            ####
################################################################################

# all positives (AP) = TP + FP
AP_grg <- c(TP, FP)

# AP genome features
AP_gfeatures <- genomeDerivedFeatures(x = AP_grg,
                                      transcriptdb = txdb,
                                      sequence = bsgenome)

# "1" is "TP", "0" is "FP"
AP_gfeatures$class <- as.factor(c(rep(1, nrow(TP)), rep(0, nrow(FP))))


selection_rm_gc_dist <- readRDS("./rds/selection_SYSY_rm_gc_dist.rds")

imp_features <- selection_rm_gc_dist$least_imp_arr[c((which.max(selection_rm_gc_dist$h2o_pred_minus_auc_arr) + 1):length(selection_rm_gc_dist$least_imp_arr))]
final_features <- rev(c(imp_features, "log2_length_exons"))

AP_gfeatures_final <- AP_gfeatures[, c(final_features, "class")]

################################################################################
####                        Fit the randomForest model                      ####
################################################################################

### Random select training set(80%) and testing set(20%)
set.seed(110)
ind_TP <- sample(1:nrow(TP), round(nrow(TP)*0.8))
ind_FP <- sample((nrow(TP)+1):(nrow(TP)+nrow(FP)), round(nrow(FP)*0.8))
ind <- c(ind_TP, ind_FP)

training_set_geo <- AP_gfeatures_final[ind,]
testing_set_geo <- AP_gfeatures_final[-ind,]

set.seed(1234)
rf_model <- randomForest(class ~ ., data = training_set_geo)
rf_pred <- predict(rf_model, newdata = testing_set_geo, type = "prob")
rf_auc = roc(testing_set_geo$class, rf_pred[,1])

# Calculate AUPRC and AUROC using the PRROC package
fg <- rf_pred[,2][testing_set_geo$class == 1]
bg <- rf_pred[,2][testing_set_geo$class == 0]
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)

### save
results <- list(rf_auc, roc, pr)
names(results) <- c("pROC_roc", "PRROC_roc", "PRROC_rf")

saveRDS(results, "./rds/final_models.rds")

################################################################################
####                            AUROC & AUPRC plot                          ####
################################################################################

library(ggplot2)
library(pROC)
library(PRROC)

final_models <- readRDS("./rds/final_models.rds")

ggplot(final_models$ROC, aes(x = X1, y = X2, color = class)) + 
  geom_line() + 
  xlab("False Positive Rate") + 
  ylab("True Positive Rate") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(colour = "black", size = 12),
        legend.position = c(0.61, 0.28),
        legend.title=element_blank(),
        legend.text=element_text(colour = "black", size = 10)) +
  scale_color_manual(values = c("dodgerblue4",
                                "darkorchid3",
                                "#d75427",
                                "darkolivegreen4"),
                     breaks = c("Abcam", "NEB", "SYSY", "Mouse"), 
                     labels = c("Abcam (AUC = 0.9708)", "NEB (AUC = 0.9876)", "SYSY (AUC = 0.9852)", "Mouse (AUC = 0.9862)"))

ggsave("~/plots/AUROC.pdf", width = 3.2, height = 2.8)


ggplot(final_models$PR, aes(x = X1, y = X2, color = class)) + 
  geom_line() + 
  xlab("Recall") + 
  ylab("Precision") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(colour = "black", size = 12),
        legend.position = c(0.38, 0.28),
        legend.title=element_blank(),
        legend.text=element_text(colour = "black", size = 10)) +
  scale_color_manual(values = c("dodgerblue4",
                                "darkorchid3",
                                "#d75427",
                                "darkolivegreen4"),
                     breaks = c("Abcam", "NEB", "SYSY", "Mouse"), 
                     labels = c("Abcam (AP = 0.9419)", "NEB (AP = 0.9959)", "SYSY (AP = 0.9900)", "Mouse (AP = 0.9941)"))

ggsave("~/plots/AUPRC.pdf", width = 3.2, height = 2.8)




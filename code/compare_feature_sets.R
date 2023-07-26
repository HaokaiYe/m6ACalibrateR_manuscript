
########################################################################################################################
####                                  Comparative analysis of machine learning models                               ####
########################################################################################################################

library(h2o)
library(pROC)
library(PRROC)
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

# AP sequence features
AP_101_grg <- AP_grg + 50

AP_seq_onehot <- sequenceDerivedFeatures(AP_101_grg, 
                                         bsgenome, 
                                         encoding = "onehot")

# combine features
AP_seq <- AP_seq_onehot
AP_geo <- AP_gfeatures
AP_seq_geo <- cbind(AP_gfeatures, AP_seq_onehot)


# "1" is "TP", "0" is "FP"
AP_seq$class <- as.factor(c(rep(1, nrow(TP)), rep(0, nrow(FP))))
AP_geo$class <- as.factor(c(rep(1, nrow(TP)), rep(0, nrow(FP))))
AP_seq_geo$class <- as.factor(c(rep(1, nrow(TP)), rep(0, nrow(FP))))


################################################################################
####                              Fit the model                             ####
################################################################################

# Random select training set(80%) and testing set(20%)
set.seed(110)
ind_TP <- sample(1:nrow(TP), round(nrow(TP)*0.8))
ind_FP <- sample((nrow(TP)+1):(nrow(TP)+nrow(FP)), round(nrow(FP)*0.8))
ind <- c(ind_TP, ind_FP)

# seq
training_set_seq <- AP_seq[ind,]
testing_set_seq <- AP_seq[-ind,]  
# geo
training_set_geo <- AP_geo[ind,]
testing_set_geo <- AP_geo[-ind,]
# seq + geo
training_set_seq_geo <- AP_seq_geo[ind,]
testing_set_seq_geo <- AP_seq_geo[-ind,]


### h2o model
h2o.init()

# seq
h2o_training_set_seq <- as.h2o(training_set_seq)
h2o_testing_set_seq <- as.h2o(testing_set_seq)  
# geo
h2o_training_set_geo <- as.h2o(training_set_geo)
h2o_testing_set_geo <- as.h2o(testing_set_geo)
# seq + geo
h2o_training_set_seq_geo <- as.h2o(training_set_seq_geo)
h2o_testing_set_seq_geo <- as.h2o(testing_set_seq_geo)


### randomForest
##############################################################
####                          seq                         ####
##############################################################
h2o_fit_seq <- h2o.randomForest(y = "class", training_frame = h2o_training_set_seq)
h2o_pred_seq <- as.data.frame(h2o.predict(h2o_fit_seq, newdata = h2o_testing_set_seq))
h2o_pred_seq_auc = roc(testing_set_seq$class, h2o_pred_seq$p0)

# Calculate ROC and PR curves
seq_fg <- h2o_pred_seq$p1[testing_set_seq$class == 1]
seq_bg <- h2o_pred_seq$p1[testing_set_seq$class == 0]
seq_roc <- roc.curve(scores.class0 = seq_fg, scores.class1 = seq_bg, curve = T)
seq_pr <- pr.curve(scores.class0 = seq_fg, scores.class1 = seq_bg, curve = T)

##############################################################
####                          geo                         ####
##############################################################
h2o_fit_geo <- h2o.randomForest(y = "class", training_frame = h2o_training_set_geo)
h2o_pred_geo <- as.data.frame(h2o.predict(h2o_fit_geo, newdata = h2o_testing_set_geo))
h2o_pred_geo_auc = roc(testing_set_geo$class, h2o_pred_geo$p0)

# Calculate ROC and PR curves
geo_fg <- h2o_pred_geo$p1[testing_set_geo$class == 1]
geo_bg <- h2o_pred_geo$p1[testing_set_geo$class == 0]
geo_roc <- roc.curve(scores.class0 = geo_fg, scores.class1 = geo_bg, curve = T)
geo_pr <- pr.curve(scores.class0 = geo_fg, scores.class1 = geo_bg, curve = T)

##############################################################
####                       seq + geo                      ####
##############################################################
h2o_fit_seq_geo <- h2o.randomForest(y = "class", training_frame = h2o_training_set_seq_geo)
h2o_pred_seq_geo <- as.data.frame(h2o.predict(h2o_fit_seq_geo, newdata = h2o_testing_set_seq_geo))
h2o_pred_seq_geo_auc = roc(testing_set_seq_geo$class, h2o_pred_seq_geo$p0)

# Calculate ROC and PR curves
seq_geo_fg <- h2o_pred_seq_geo$p1[testing_set_seq_geo$class == 1]
seq_geo_bg <- h2o_pred_seq_geo$p1[testing_set_seq_geo$class == 0]
seq_geo_roc <- roc.curve(scores.class0 = seq_geo_fg, scores.class1 = seq_geo_bg, curve = T)
seq_geo_pr <- pr.curve(scores.class0 = seq_geo_fg, scores.class1 = seq_geo_bg, curve = T)


##############################################################
####                          save                        ####
##############################################################
results <- list(h2o_pred_seq_roc, seq_roc, seq_pr,
                h2o_pred_geo_roc, geo_roc, geo_pr,
                h2o_pred_seq_geo_roc, seq_geo_roc, seq_geo_pr)
names(results) <- c("seq_pROC_roc", "seq_PRROC_roc", "seq_PRROC_rf",
                    "geo_pROC_roc", "geo_PRROC_roc", "geo_PRROC_rf",
                    "seq_geo_pROC_roc", "seq_geo_PRROC_roc", "seq_geo_PRROC_rf")

saveRDS(results, "./rds/compare_feature_sets.rds")

################################################################################
####                                Plot curves                             ####
################################################################################

library(pROC)

# Read the model selection results from an RDS file
compare_feature_sets <- readRDS("./rds/compare_feature_sets.rds")

roc_df = rbind(data.frame(compare_feature_sets$seq_PRROC_roc$curve, class = rep("Seq", nrow(compare_feature_sets$seq_PRROC_roc$curve))),
               data.frame(compare_feature_sets$geo_PRROC_roc$curve, class = rep("Geo", nrow(compare_feature_sets$geo_PRROC_roc$curve))),
               data.frame(compare_feature_sets$seq_geo_PRROC_roc$curve, class = rep("Seq+Geo", nrow(compare_feature_sets$seq_geo_PRROC_roc$curve))))


pr_df = rbind(data.frame(compare_feature_sets$seq_PRROC_rf$curve, class = rep("Seq", nrow(compare_feature_sets$seq_PRROC_rf$curve))),
              data.frame(compare_feature_sets$geo_PRROC_rf$curve, class = rep("Geo", nrow(compare_feature_sets$geo_PRROC_rf$curve))),
              data.frame(compare_feature_sets$seq_geo_PRROC_rf$curve, class = rep("Seq+Geo", nrow(compare_feature_sets$seq_geo_PRROC_rf$curve))))


# ROC curve
ggplot(roc_df, aes(x = X1, y = X2, color = class)) + 
        geom_line() + 
        xlab("False Positive Rate") + 
        ylab("True Positive Rate") +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(colour = "black", size = 12),
              axis.text = element_text(colour= "black", size = 10),
              legend.position = c(0.63, 0.2),
              legend.title=element_blank(),
              legend.background = element_blank(),
              legend.text=element_text(colour = "black", size = 10)) +
        scale_color_manual(values = c("#85B0A9",
                                      "#E39CB1",
                                      "#EFDA8E"),
                           breaks = c("Geo", "Seq", "Seq+Geo"), 
                           labels = c("Geo (AUC = 0.97)", "Seq (AUC = 0.60)", "Seq+Geo (AUC = 0.92)"))

ggsave("~/plots/cs_fea_roc.pdf", width = 3.2, height = 3)

# PR curve
ggplot(pr_df, aes(x = X1, y = X2, color = class)) + 
        geom_line() + 
        xlab("Recall") + 
        ylab("Precision") +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(colour = "black", size = 12),
              axis.text = element_text(colour= "black", size = 10),
              legend.position = c(0.4, 0.2),
              legend.title=element_blank(),
              legend.background = element_blank(),
              legend.text=element_text(colour = "black", size = 10)) +
        scale_color_manual(values = c("#85B0A9",
                                      "#E39CB1",
                                      "#EFDA8E"),
                           breaks = c("Geo", "Seq", "Seq+Geo"), 
                           labels = c("Geo (AP = 0.98)", "Seq (AP = 0.69)", "Seq+Geo (AP = 0.94)"))

ggsave("~/plots/cs_fea_prc.pdf", width = 3.2, height = 3)


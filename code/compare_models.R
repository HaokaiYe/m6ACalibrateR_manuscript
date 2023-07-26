
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
AP_all_features <- cbind(AP_gfeatures, AP_seq_onehot)

# "1" is "TP", "0" is "FP"
AP_all_features$class <- as.factor(c(rep(1, nrow(TP)), rep(0, nrow(FP))))

################################################################################
####                              Fit the model                             ####
################################################################################

# Random select training set(80%) and testing set(20%)
set.seed(110)
ind_TP <- sample(1:nrow(TP), round(nrow(TP)*0.8))
ind_FP <- sample((nrow(TP)+1):(nrow(TP)+nrow(FP)), round(nrow(FP)*0.8))
ind <- c(ind_TP, ind_FP)

training_set <- AP_all_features[ind,]
testing_set <- AP_all_features[-ind,]


# h2o model
h2o.init()

h2o_training_set <- as.h2o(training_set)
h2o_testing_set <- as.h2o(testing_set)

##############################################################
####                          glm                         ####
##############################################################
h2o_fit_glm <- h2o.glm(y = "class", training_frame = h2o_training_set, family = "binomial", link = "logit")
h2o_pred_glm <- as.data.frame(h2o.predict(h2o_fit_glm, newdata = h2o_testing_set))
h2o_pred_glm_roc = roc(testing_set$class, h2o_pred_glm$p0)

# Calculate ROC and PR curves
glm_fg <- h2o_pred_glm$p1[testing_set$class == 1]
glm_bg <- h2o_pred_glm$p1[testing_set$class == 0]
glm_roc <- roc.curve(scores.class0 = glm_fg, scores.class1 = glm_bg, curve = T)
glm_pr <- pr.curve(scores.class0 = glm_fg, scores.class1 = glm_bg, curve = T)

##############################################################
####                        xgboost                       ####
##############################################################
h2o_fit_xgboost <- h2o.xgboost(y = "class", training_frame = h2o_training_set)
h2o_pred_xgboost <- as.data.frame(h2o.predict(h2o_fit_xgboost, newdata = h2o_testing_set))
h2o_pred_xgboost_auc = roc(testing_set$class, h2o_pred_xgboost$p0)

# Calculate ROC and PR curves
xgboost_fg <- h2o_pred_xgboost$p1[testing_set$class == 1]
xgboost_bg <- h2o_pred_xgboost$p1[testing_set$class == 0]
xgboost_roc <- roc.curve(scores.class0 = xgboost_fg, scores.class1 = xgboost_bg, curve = T)
xgboost_pr <- pr.curve(scores.class0 = xgboost_fg, scores.class1 = xgboost_bg, curve = T)

##############################################################
####                      randomForest                    ####
##############################################################
h2o_fit_rf <- h2o.randomForest(y = "class", training_frame = h2o_training_set)
h2o_pred_rf <- as.data.frame(h2o.predict(h2o_fit_rf, newdata = h2o_testing_set))
h2o_pred_rf_auc = roc(testing_set$class, h2o_pred_rf$p0)

# Calculate ROC and PR curves
rf_fg <- h2o_pred_rf$p1[testing_set$class == 1]
rf_bg <- h2o_pred_rf$p1[testing_set$class == 0]
rf_roc <- roc.curve(scores.class0 = rf_fg, scores.class1 = rf_bg, curve = T)
rf_pr <- pr.curve(scores.class0 = rf_fg, scores.class1 = rf_bg, curve = T)

##############################################################
####                      deep learning                   ####
##############################################################
h2o_fit_dp <- h2o.deeplearning(y = "class", training_frame = h2o_training_set)
h2o_pred_dp <- as.data.frame(h2o.predict(h2o_fit_dp, newdata = h2o_testing_set))
h2o_pred_dp_roc = roc(testing_set$class, h2o_pred_dp$p0)

# Calculate ROC and PR curves
dp_fg <- h2o_pred_dp$p1[testing_set$class == 1]
dp_bg <- h2o_pred_dp$p1[testing_set$class == 0]
dp_roc <- roc.curve(scores.class0 = dp_fg, scores.class1 = dp_bg, curve = T)
dp_pr <- pr.curve(scores.class0 = dp_fg, scores.class1 = dp_bg, curve = T)


##############################################################
####                          save                        ####
##############################################################
results <- list(h2o_pred_glm_roc, glm_roc, glm_pr,
                h2o_pred_xgboost_roc, xgboost_roc, xgboost_pr,
                h2o_pred_rf_roc, rf_roc, rf_pr,
                h2o_pred_dp_roc, dp_roc, dp_pr)
names(results) <- c("glm_pROC_roc", "glm_PRROC_roc", "glm_PRROC_rf",
                    "xgboost_pROC_roc", "xgboost_PRROC_roc", "xgboost_PRROC_rf",
                    "rf_pROC_roc", "rf_PRROC_roc", "rf_PRROC_rf",
                    "dp_pROC_roc", "dp_PRROC_roc", "dp_PRROC_rf")

saveRDS(results, "./rds/compare_models.rds")


################################################################################
####                                Plot curves                             ####
################################################################################

library(pROC)

# Read the model selection results from an RDS file
compare_models <- readRDS("./rds/compare_models.rds")

roc_df = rbind(data.frame(compare_models$glm_PRROC_roc$curve, class = rep("GLM", nrow(compare_models$glm_PRROC_roc$curve))),
               data.frame(compare_models$xgboost_PRROC_roc$curve, class = rep("XGBoost", nrow(compare_models$xgboost_PRROC_roc$curve))),
               data.frame(compare_models$rf_PRROC_roc$curve, class = rep("RF", nrow(compare_models$rf_PRROC_roc$curve))),
               data.frame(compare_models$dp_PRROC_roc$curve, class = rep("DL", nrow(compare_models$dp_PRROC_roc$curve))))


pr_df = rbind(data.frame(compare_models$glm_PRROC_rf$curve, class = rep("GLM", nrow(compare_models$glm_PRROC_rf$curve))),
              data.frame(compare_models$xgboost_PRROC_rf$curve, class = rep("XGBoost", nrow(compare_models$xgboost_PRROC_rf$curve))),
              data.frame(compare_models$rf_PRROC_rf$curve, class = rep("RF", nrow(compare_models$rf_PRROC_rf$curve))),
              data.frame(compare_models$dp_PRROC_rf$curve, class = rep("DL", nrow(compare_models$dp_PRROC_rf$curve))))


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
              legend.position = c(0.63, 0.25),
              legend.title=element_blank(),
              legend.background = element_blank(),
              legend.text=element_text(colour = "black", size = 10)) +
        scale_color_manual(values = c("#5E519B",
                                      "#76C6A5",
                                      "#F47C42",
                                      "#900943"),
                           breaks = c("DL", "GLM", "RF", "XGBoost"), 
                           labels = c("DL (AUC = 0.83)", "GLM (AUC = 0.82)", "RF (AUC = 0.92)", "XGBoost (AUC = 0.87)"))

ggsave("~//plots/cs_mod_roc.pdf", width = 3.2, height = 3)

# PR curve
pr_df_clean = pr_df[-c(56500:56743),]
ggplot(pr_df_clean, aes(x = X1, y = X2, color = class)) + 
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
              legend.position = c(0.4, 0.25),
              legend.title=element_blank(),
              legend.background = element_blank(),
              legend.text=element_text(colour = "black", size = 10)) +
        scale_color_manual(values = c("#5E519B",
                                      "#76C6A5",
                                      "#F47C42",
                                      "#900943"),
                           breaks = c("DL", "GLM", "RF", "XGBoost"), 
                           labels = c("DL (AP = 0.87)", "GLM (AP = 0.85)", "RF (AP = 0.94)", "XGBoost (AP = 0.91)"))

ggsave("~/plots/cs_mod_prc.pdf", width = 3.2, height = 3)


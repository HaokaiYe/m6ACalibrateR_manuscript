
########################################################################################################################
####                                      Cross validation for three models                                         ####
########################################################################################################################

library(pROC)
library(randomForest)
library(predictiveFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 
bsgenome <- BSgenome.Hsapiens.UCSC.hg38


################################################################################
####                             Load rf models                             ####
################################################################################

rf_Abcam <- readRDS("./rds/Abcam.rds")
rf_NEB <- readRDS("./rds/NEB.rds")
rf_SYSY <- readRDS("./rds/SYSY.rds")
rf_mouse <- readRDS("./rds/mouse.rds")

##### human data sets
antibody <- c("Abcam", "NEB", "SYSY")

crossValidAUC_rf <- data.frame()

for (t in 1:length(antibody)) {
  
  
  ################################################################################
  ####                     Load mRNA peaks and IVT peaks                      ####
  ################################################################################
  
  peaks_mRNA <- readRDS("./rds/", antibody[t],"_peaks_mRNA.rds")
  peaks_IVT <- readRDS("./rds/", antibody[t],"_peaks_IVT.rds")
  
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
                                        transcriptdb = txdb)
  
  # "1" is "TP", "0" is "FP"
  AP_gfeatures$class <- as.factor(c(rep(1, nrow(TP)), rep(0, nrow(FP))))

  ################################################################################
  ####                        random forest predicting                        ####
  ################################################################################
  
  # Random select training set(80%) and testing set(20%)
  set.seed(110)
  ind_TP <- sample(1:nrow(TP), round(nrow(TP)*0.8))
  ind_FP <- sample((nrow(TP)+1):(nrow(TP)+nrow(FP)), round(nrow(FP)*0.8))
  ind <- c(ind_TP, ind_FP)
  
  # training & testing sets
  training_set_geo <- AP_gfeatures[ind,]
  testing_set_geo <- AP_gfeatures[-ind,]
  
  # predicting
  rf_Abcam_pred <- predict(rf_Abcam, newdata = testing_set_geo, type = "prob")
  rf_NEB_pred <- predict(rf_NEB, newdata = testing_set_geo, type = "prob")
  rf_SYSY_pred <- predict(rf_SYSY, newdata = testing_set_geo, type = "prob")
  rf_mouse_pred <- predict(rf_mouse, newdata = testing_set_geo, type = "prob")

  # AUC
  rf_Abcam_auc = roc(testing_set_geo$class, rf_Abcam_pred[,1])$auc[1]
  rf_NEB_auc = roc(testing_set_geo$class, rf_NEB_pred[,1])$auc[1]
  rf_SYSY_auc = roc(testing_set_geo$class, rf_SYSY_pred[,1])$auc[1]
  rf_mouse_auc = roc(testing_set_geo$class, rf_mouse_pred[,1])$auc[1]
  
  
  print(paste0("random forest Abcam predict ", antibody[t], ": ", rf_Abcam_auc))
  print(paste0("random forest NEB predict ", antibody[t], ": ", rf_NEB_auc))
  print(paste0("random forest SYSY predict ", antibody[t], ": ", rf_SYSY_auc))
  print(paste0("random forest mouse predict ", antibody[t], ": ", rf_mouse_auc))
  
  
  rf_auc_vec <- c(rf_Abcam_auc, rf_NEB_auc, rf_SYSY_auc, rf_mouse_auc)
  crossValidAUC_rf <- rbind(crossValidAUC_rf, rf_auc_vec)
}

colnames(crossValidAUC_rf) <- c(antibody, "mouse")
rownames(crossValidAUC_rf) <- c(antibody)

saveRDS(crossValidAUC_rf, "./rds/cross_validation.rds")


################################################################################
####                                Plot heatmap                            ####
################################################################################

library(aplot)
library(ggtree)
library(ggplot2)
library(reshape2)

cross_validation = readRDS("./rds/cross_validation.rds")

data = melt(as.matrix(cross_validation))
p = ggplot(data = data, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = NA, width = 0.9, height = 0.9) + 
  #geom_text(aes(label = sprintf("%.4f", round(value, 4))), color = "#ffffff", size = 2.8) +
  scale_fill_gradient(name="AUC", low = "#F8D7C9", high = "#860D24", limits = c(0.7, 1), breaks = c(0.7, 0.8, 0.9, 1)) +
  xlab("Testing dataset") +
  ylab("Training dataset") +
  scale_x_discrete(limits = c("Mouse", "Abcam", "NEB", "SYSY"), expand = c(0,0)) +
  scale_y_discrete(limits = c("Mouse", "Abcam", "NEB", "SYSY"), expand = c(0,0)) +
  theme(axis.ticks = element_blank(),
        axis.title.x = element_text(colour= "black", size = 12),
        axis.title.y = element_text(colour= "black", size = 12),
        axis.text.x = element_text(colour= "black", size = 10, angle = 90, vjust = 0.5, margin = margin(0, 0, 3, 0)),
        axis.text.y = element_text(colour= "black", size = 10),
        legend.title = element_text(vjust = 1.3, hjust = -1.6),
        #legend.direction = "horizontal",
        #legend.position = "top",
        #legend.margin = margin(t = 0, r = 0, b = -10, l = 0),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(1.0, "cm"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))


hc <- hclust(dist(cross_validation))
phr <- ggtree(hc, layout="rectangular",branch.length="none") + layout_dendrogram()
insert_top(p, phr, height =.1)

ggsave("~/plots/cross_validation.pdf", width = 3.6, height = 3.1)



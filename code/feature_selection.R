
########################################################################################################################
####                                                 Feature selection                                              ####
########################################################################################################################

library(pROC)
library(h2o)

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
####                      Extract genomic features only                     ####
################################################################################

# all positives (AP) = TP + FP
AP_grg <- c(TP, FP)

# AP genome features
AP_gfeatures <- genomeDerivedFeatures(x = AP_grg,
                                      transcriptdb = txdb,
                                      sequence = bsgenome)
  
# "1" is "TP", "0" is "FP"
AP_gfeatures$class <- as.factor(c(rep(1, nrow(TP)), rep(0, nrow(FP))))

AP_gfeatures_colnames <- colnames(AP_gfeatures)
AP_gfeatures_rm_gc_dist <- AP_gfeatures[, which(!grepl("GC|dist", AP_gfeatures_colnames))]


################################################################################
####                              Fit the model                             ####
################################################################################

### Random select training set(80%) and testing set(20%)
set.seed(110)
ind_TP <- sample(1:nrow(TP), round(nrow(TP)*0.8))
ind_FP <- sample((nrow(TP)+1):(nrow(TP)+nrow(FP)), round(nrow(FP)*0.8))
ind <- c(ind_TP, ind_FP)

training_set_geo <- AP_gfeatures_rm_gc_dist[ind,]
testing_set_geo <- AP_gfeatures_rm_gc_dist[-ind,]


### h2o model
h2o.init()

h2o_training_set_geo <- as.h2o(training_set_geo)
h2o_testing_set_geo <- as.h2o(testing_set_geo)


### randomForest
h2o_fit_geo <- h2o.randomForest(y = "class", training_frame = h2o_training_set_geo, seed = 1234)
h2o_pred_geo <- as.data.frame(h2o.predict(h2o_fit_geo, newdata = h2o_testing_set_geo))
h2o_pred_geo_auc = roc(testing_set_geo$class, h2o_pred_geo$p0)$auc[1]

################################################################################
####                        Perform feature selection                       ####
################################################################################

### variable importance
variable <- h2o.varimp(h2o_fit_geo)$variable
len <- length(variable)

variable_plus <- variable
variable_minus <- variable


most_imp_arr <- c()
least_imp_arr <- c()
h2o_pred_plus_auc_arr <- c()
h2o_pred_minus_auc_arr <- c()

for (i in 1:len) {
  print(date())
  h2o.removeAll()
  
  
  # Add the most important feature one by one and each time select the first one
  most_imp <- variable_plus[1]
  most_imp_arr <- c(most_imp_arr, most_imp)
  
  print(paste0("It's adding ", i, " variables. ", most_imp, " is doing!"))
  
  # Regenerate the training and testing sets
  training_set_plus <- dplyr::select(training_set_geo, c(most_imp_arr, "class"))
  testing_set_plus <- dplyr::select(testing_set_geo, c(most_imp_arr, "class"))
  
  h2o_training_set_plus <- as.h2o(training_set_plus)
  h2o_testing_set_plus <- as.h2o(testing_set_plus)
  
  
  # Model the data and calculate the AUC
  h2o_fit_plus <- h2o.randomForest(y = "class", training_frame = h2o_training_set_plus, seed = 1234)
  h2o_pred_plus <- as.data.frame(h2o.predict(h2o_fit_plus, newdata = h2o_testing_set_plus))
  h2o_pred_plus_auc = roc(testing_set_plus$class, h2o_pred_plus$p0)$auc[1]
  
  
  # Store the results
  h2o_pred_plus_auc_arr <- c(h2o_pred_plus_auc_arr, h2o_pred_plus_auc)
  
  
  # Recalculate the variable importance
  if(length(most_imp_arr) < len) {
    
    # Remove the most important feature
    training_set_plus_re <- dplyr::select(training_set_geo, -c(most_imp_arr))
    
    h2o_training_set_plus_re <- as.h2o(training_set_plus_re)
    
    # Reassign the variable importance
    h2o_fit_plus_re <- h2o.randomForest(y = "class", training_frame = h2o_training_set_plus_re)
    variable_plus <- h2o.varimp(h2o_fit_plus_re)$variable
  }      
  
  
  
  # Remove the least important feature one by one and each time select the last one
  if(length(most_imp_arr) < len) {
    least_imp <- tail(variable_minus, 1)
    least_imp_arr <- c(least_imp_arr, least_imp)
    
    print(paste0("It's deleting ", i, " variables. ", least_imp, " is doing!"))
    
    # Regenerate the training and testing sets
    training_set_minus <- dplyr::select(training_set_geo, -c(least_imp_arr))
    testing_set_minus <- dplyr::select(testing_set_geo, -c(least_imp_arr))
    
    h2o_training_set_minus <- as.h2o(training_set_minus)
    h2o_testing_set_minus <- as.h2o(testing_set_minus)
    
    
    # Model the data and calculate the AUC
    h2o_fit_minus <- h2o.randomForest(y = "class", training_frame = h2o_training_set_minus, seed = 1234)
    h2o_pred_minus <- as.data.frame(h2o.predict(h2o_fit_minus, newdata = h2o_testing_set_minus))
    h2o_pred_minus_auc = roc(testing_set_minus$class, h2o_pred_minus$p0)$auc[1]
    
    
    # Store the results
    h2o_pred_minus_auc_arr <- c(h2o_pred_minus_auc_arr, h2o_pred_minus_auc)
    
    
    # Reassign the variable importance
    variable_minus <- h2o.varimp(h2o_fit_minus)$variable
    
  } else {
    least_imp_arr <- c("All features", least_imp_arr)
    
    h2o_pred_minus_auc_arr <- c(h2o_pred_geo_auc, h2o_pred_minus_auc_arr)
  }
  
}

### save
result = data.frame(variable, most_imp_arr, h2o_pred_plus_auc_arr, least_imp_arr, h2o_pred_minus_auc_arr)
saveRDS(result, "./rds/feature_selection.rds")


########################################################################################################################
####                                            Plot feature selection                                              ####
########################################################################################################################

library(ggplot2)
library(ggrepel)

feature_selection <- readRDS("./rds/feature_selection.rds")

ggplot(feature_selection, aes(x = features_num, y = auc, color = antibody)) + 
  geom_point(size = 0.4) +
  geom_point(data = feature_selection[which(feature_selection$max == TRUE),], size = 1.5) +
  geom_point(data = feature_selection[which(feature_selection$features_num == 1),], shape = 1, size = 2.5) +
  geom_point(data = feature_selection[which(feature_selection$features_num == 2),], shape = 1, size = 2) +
  geom_line(aes(color=antibody), size = 0.5) + 
  ylim(c(0.85, 1)) +
  ylab("Performence of random forest model (AUC)") +
  xlab("Number of features") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_rect(size = 0.8),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.text.y = element_text(colour= "black", size = 10),
        axis.text.x = element_text(colour = "black", size = 10),
        legend.position = 'none',
        axis.title.y  =element_text(colour= "black", size=12, vjust = 2),
        axis.title.x = element_text(colour= "black", size=12),
        strip.text = element_text(colour = "black", size = 10),
        strip.background = element_rect(size = 0.8)) +
  geom_label_repel(aes(label = label), 
                   size = 3.6, 
                   box.padding = 0.5, 
                   min.segment.length = 0,
                   segment.curvature = 0) +
  geom_text_repel(aes(label = text), 
                  size = 3.6, 
                  min.segment.length = 0,
                  nudge_x = 1,
                  nudge_y = 0.1) +
  scale_color_manual(values = c("dodgerblue4",
                                "darkorchid3",
                                "#d75427",
                                "darkolivegreen4")) +
  facet_wrap(facet = antibody ~ ., ncol = 1) +
  geom_rect(aes(xmin = 20, xmax = 40, ymin = 0.86, ymax = 0.9), fill = NA, size = 0.05, linetype = "longdash") +
  geom_text(data = feature_selection[which(feature_selection$max == TRUE),], aes(label = paste0("Best: [ ", features_num, ", ", sprintf("%.4f", round(auc, 4)), " ]")), x = 30, y = 0.88, size = 3.6)


ggsave("~/plots/feature_select.pdf", width = 3.6, height = 7.6)




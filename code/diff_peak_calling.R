#! /usr/bin/env Rscript

library(exomePeak2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg38 <- BSgenome.Hsapiens.UCSC.hg38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


data_dir <- "~/merip/Data/"
writer_list <- c("KIAA1429", "METTL3", "METTL14", "WTAP", "METTL16", "M3M4")
peak_dir <- grep("PEAK", list.files(data_dir), value = TRUE)
peak_dir <- unlist(lapply(writer_list, function(writer) {
  grep(writer, peak_dir, value = TRUE)
}))

library(SummarizedExperiment)
bamDir = "/data/heterogenous/"
dir = "~/merip/diffPeakAll/"
for (i in 1:10) {
  peak_se <- readRDS(file.path(data_dir, peak_dir[i]))
  
  indxCtrl = peak_se$Treatment == "Ctrl"
  indxTreat = peak_se$Treatment != "Ctrl"
  
  IP_CTRL = peak_se$SRR_IP[indxCtrl]
  INPUT_CTRL = peak_se$SRR_input[indxCtrl]
  IP_TREAT = peak_se$SRR_IP[indxTreat]
  INPUT_TREAT = peak_se$SRR_input[indxTreat]
  
  IP_BAM_CTRL = paste0(bamDir, IP_CTRL, ".bam")
  INPUT_BAM_CTRL = paste0(bamDir, INPUT_CTRL, ".bam")
  IP_BAM_TREAT = paste0(bamDir, IP_TREAT, ".bam")
  INPUT_BAM_TREAT = paste0(bamDir, INPUT_TREAT, ".bam")
  
  # diff peak calling
  res_diff <- exomePeak2(bam_ip = IP_BAM_CTRL,
                         bam_input = INPUT_BAM_CTRL,
                         bam_ip_treated = IP_BAM_TREAT,
                         bam_input_treated = INPUT_BAM_TREAT,
                         txdb = txdb,
                         genome = hg38,
                         save_dir = dir,
                         experiment_name = paste0("diffPeakAll_", i))
  date()
  print(i)
}


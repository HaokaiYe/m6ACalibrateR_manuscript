
########################################################################################################################
####                                           Peak calling using exomePeak2                                        ####
########################################################################################################################

library(exomePeak2)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg38 <- BSgenome.Hsapiens.UCSC.hg38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

################################################################################
####                              IVT samples                               ####
################################################################################

IP_BAM_IVT = c("~/bam/SRR14765585.bam", "~/bam/SRR14765587.bam")
INPUT_BAM_IVT = c("~/bam/SRR14765584.bam", "~/bam/SRR14765586.bam")
dir_IVT = "~/peak/IVT/"

res_IVT <- exomePeak2(bam_ip = IP_BAM_IVT,
                      bam_input = INPUT_BAM_IVT,
                      txdb = txdb,
                      genome = hg38,
                      p_cutoff = 1e-10,
                      motif_based = TRUE,
                      save_dir = dir_IVT)

################################################################################
####                              mRNA samples                              ####
################################################################################

IP_BAM_mRNA = c("~/bam/SRR14765589.bam", "~/bam/SRR14765591.bam")
INPUT_BAM_mRNA = c("~/bam/SRR14765588.bam", "~/bam/SRR14765590.bam")
dir_mRNA = "~/peak/mRNA/"

res_mRNA <- exomePeak2(bam_ip = IP_BAM_mRNA,
                       bam_input = INPUT_BAM_mRNA,
                       txdb = txdb,
                       genome = hg38,
                       p_cutoff = 1e-10,
                       motif_based = TRUE,
                       save_dir = dir_mRNA)


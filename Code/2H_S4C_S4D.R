library(openxlsx)
library(org.Hs.eg.db)
library(goseq)
library(pheatmap)
library(viridis)
library(SummarizedExperiment)
library(tidyr)
library(ggplot2)

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################





############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

# PARAMETER
dataDir <- file.path("../Data")
outputDir <- file.path("../CibersortX_heatmaps")
dir.create(outputDir)

# INPUT FILES
signature.2H <- "CIBERSORTx_Job61_dataset3_V2_TPM_reference_sample_inferred_phenoclasses.CIBERSORTx_Job61_dataset3_V2_TPM_reference_sample_inferred_refsample.bm.K999.txt"
signature.S4C <- "CIBERSORTx_Job64_Richards_Dev_TPM_reference_sample_inferred_phenoclasses.CIBERSORTx_Job64_Richards_Dev_TPM_reference_sample_inferred_refsample.bm.K999.txt"
signature.S4D <- "CIBERSORTx_Job66_Richards_Inj_TPM_reference_sample_inferred_phenoclasses.CIBERSORTx_Job66_Richards_Inj_TPM_reference_sample_inferred_refsample.bm.K999.txt"


#####################################################################################
# SIGNATURE MATRIX HEATMAPS

###########
# FIGURE 2H
setwd(dataDir)
signMat <- read.delim(signature.2H)
rownames(signMat) <- signMat$NAME
signMat <- signMat[, -1]
signMat <- t(scale(t(log2(signMat))))

orderMat <- lapply(1:ncol(signMat), function(i) signMat[, i] - rowMeans(signMat[, -i]))
orderMat <- do.call(cbind, orderMat)
rownames(orderMat) <- rownames(signMat)
colnames(orderMat) <- colnames(signMat)
orderMat[orderMat < 0] <- 0
newOrder <- order(-orderMat[,1],
                  -orderMat[,2],
                  -orderMat[,3],
                  -orderMat[,4],
                  -orderMat[,5],
                  -orderMat[,6])

paletteLength <- 50
myColor <- magma(paletteLength)
myBreaks <- c(seq(0, max(orderMat), length.out=ceiling(paletteLength)))

setwd(outputDir)
pheatmap(orderMat[newOrder, ], color = myColor, breaks = myBreaks, filename = "ds3_V2_signature_heatmap.pdf",
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE, clustering_distance_cols = "correlation",
         cellwidth = 10, fontsize_col = 8,
         height = 4,
         border_color = NA
)

##############
# FIGURE S4C
setwd(dataDir)
signMat <- read.delim(signature.S4C)
rownames(signMat) <- signMat$NAME
signMat <- signMat[, -1]
signMat <- t(scale(t(log2(signMat))))

orderMat <- lapply(1:ncol(signMat), function(i) signMat[, i] - rowMeans(signMat[, -i]))
orderMat <- do.call(cbind, orderMat)
rownames(orderMat) <- rownames(signMat)
colnames(orderMat) <- colnames(signMat)
orderMat[orderMat < 0] <- 0
newOrder <- order(-orderMat[,1],
                  -orderMat[,2],
                  -orderMat[,3])

paletteLength <- 50
myColor <- magma(paletteLength)
myBreaks <- c(seq(0, max(orderMat), length.out=ceiling(paletteLength)))

setwd(outputDir)
pheatmap(orderMat[newOrder, ], color = myColor, breaks = myBreaks, filename = "Richards_dev_signature_heatmap.pdf",
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE, clustering_distance_cols = "correlation",
         cellwidth = 10, fontsize_col = 8,
         height = 4,
         border_color = NA
)

##############
# FIGURE S4D
setwd(dataDir)
signMat <- read.delim(signature.S4C)
rownames(signMat) <- signMat$NAME
signMat <- signMat[, -1]
signMat <- t(scale(t(log2(signMat))))

orderMat <- lapply(1:ncol(signMat), function(i) signMat[, i] - rowMeans(signMat[, -i]))
orderMat <- do.call(cbind, orderMat)
rownames(orderMat) <- rownames(signMat)
colnames(orderMat) <- colnames(signMat)
orderMat[orderMat < 0] <- 0
newOrder <- order(-orderMat[,1],
                  -orderMat[,2],
                  -orderMat[,3])

paletteLength <- 50
myColor <- magma(paletteLength)
myBreaks <- c(seq(0, max(orderMat), length.out=ceiling(paletteLength)))

setwd(outputDir)
pheatmap(orderMat[newOrder, ], color = myColor, breaks = myBreaks, filename = "Richards_inj_signature_heatmap.pdf",
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE, clustering_distance_cols = "correlation",
         cellwidth = 10, fontsize_col = 8,
         height = 4,
         border_color = NA
)




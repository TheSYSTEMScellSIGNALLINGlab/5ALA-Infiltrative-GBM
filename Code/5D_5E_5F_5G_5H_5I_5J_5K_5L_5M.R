library(R.matlab)
library(dplyr)
library(Seurat)
library(patchwork)

x = readMat("gbm_sparse.mat")
# write.csv(colnames(df),"genes.csv")
df = x$RAW.COUNTS

df = as.data.frame(as.matrix(df))

#barcodes
# bars = x$barcodes
# bars = as.matrix(bars)
# first = bars[[1]][[1]]
# second = bars[[2]][[1]]
# third = bars[[3]][[1]]
# fourth = bars[[4]][[1]]
# fifth = bars[[5]][[1]]
# sixth = bars[[6]][[1]]
# seventh = bars[[7]][[1]]
# br = rbind(first, second, third, fourth, fifth, sixth, seventh)
# write.csv(br, "barcodes.csv")

codes = read.csv("barcodes.csv")

gene_name = read.csv("genes.csv")

df = t(df)
dim(df)
colnames(df) = codes$Barcodes
rownames(df) = gene_name$Genes
df = as.matrix(df)
# df = cbind(x$sample, df)



exp.mat <- Matrix::Matrix(as.matrix(df), sparse = T)
# if(.Platform$OS.type == "windows") withAutoprint({
#   memory.size()
#   memory.size(TRUE)
#   memory.limit()
# })
# memory.limit(size=56000)
# write.csv(df,"matrix.csv")
# remove(sp)

# annot = readMat("annotated_cancer_data.mat")
# remotes::install_github("carmonalab/UCell")
library(Seurat)
library(UCell)
library(data.table)
set.seed(123)

library(readxl)
gene_sig = read_excel("ALA pos gene list_final.xlsx")


signatures <- list(ALA_sig = c(gene_sig$ALA_pos_signature_251),
                   GP = c(gene_sig$GP),
                   OLC = c(gene_sig$OLC),
                   Others = c(gene_sig$Others))
u.scores <- ScoreSignatures_UCell(exp.mat, features = signatures)
write.csv(u.scores, "U_scores_Courtier.csv")
u.scores[1:8, 1]
dim(u.scores)

library(reshape2)
library(ggplot2)
melted <- melt(u.scores)
colnames(melted) <- c("Cell", "Signature", "Uscore")
p <- ggplot(melted, aes(x = Signature, y = Uscore)) + geom_violin(aes(fill = Signature), 
                                                                  scale = "width") + geom_boxplot(width = 0.1) + theme_bw() + theme(axis.text.x = element_blank())
p

set.seed(123)
ranks <- StoreRankings_UCell(exp.mat)

ranks[1:5, 1:5]

seurat.object <- CreateSeuratObject(counts = exp.mat, project = "Courtier")

seurat.object <- AddMetaData(seurat.object, metadata = as.data.frame(u.scores))
head(seurat.object@meta.data)

seurat.object <- NormalizeData(seurat.object, verbose = FALSE)
seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 1000, 
                                      verbose = FALSE)

seurat.object <- ScaleData(seurat.object)
seurat.object <- RunPCA(seurat.object, features = seurat.object@assays$RNA@var.features, 
                        verbose = FALSE)
seurat.object <- RunTSNE(seurat.object, reduction = "pca", dims = 1:20, seed.use = 123, 
                         verbose = FALSE)
signature.names <- paste0(names(signatures), "_UCell")


# Figure 5D, 5E, 5F, 5G

FeaturePlot(seurat.object, reduction = "tsne",
            features = signature.names, ncol = 3,
            order = T,
            cols = c("darkred", "orange", "yellow","white"))




##########################################
############# Seurat Analysis ############
##########################################

seurat.object <- CreateSeuratObject(counts = exp.mat, project = "Courtier")
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)


seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat.object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat.object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(seurat.object)
seurat.object <- ScaleData(seurat.object, features = all.genes)


seurat.object <- RunPCA(seurat.object, features = VariableFeatures(object = seurat.object))

# Examine and visualize PCA results a few different ways
print(seurat.object[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(seurat.object, dims = 1:2, reduction = "pca")

DimPlot(seurat.object, reduction = "pca")

DimHeatmap(seurat.object, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat.object, dims = 1:15, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
seurat.object <- JackStraw(seurat.object, num.replicate = 100)
seurat.object <- ScoreJackStraw(seurat.object, dims = 1:20)

JackStrawPlot(seurat.object, dims = 1:15)
ElbowPlot(seurat.object)


seurat.object <- FindNeighbors(seurat.object, dims = 1:10)
seurat.object <- FindClusters(seurat.object, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat.object), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
seurat.object <- RunTSNE(seurat.object, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(seurat.object, reduction = "tsne")

# Figure 5H, 5I, 5J, 5K, 5L, 5M

FeaturePlot(seurat.object, reduction = "tsne",
            features = c("CD44", "AQP4", "FAM107A", "SOX9", "GLI3", "TIMP1"))

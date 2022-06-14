# remotes::install_github("carmonalab/UCell")

library(Seurat)
library(openxlsx)
library(UCell)
library(data.table)
set.seed(123)


## Load expression matrix and metadata
exp.mat <- fread("Neftel.csv")
exp.mat <- unique(exp.mat, by = "GENE")

Rw <- exp.mat[,1]
Rw <- as.matrix(Rw)

exp.mat <- as.matrix(exp.mat[,1:4917])
row.names(exp.mat) <- Rw
exp.mat <- exp.mat[,2:4917]
exp.mat <- as.matrix(exp.mat)


# Make sparse
exp.mat <- Matrix::Matrix(as.matrix(exp.mat), sparse = T)
# saveRDS(exp.mat, "Neftel.rds")


set.seed(123)
gene_sig <- read.xlsx("ALA pos gene list_final.xlsx")
signatures <- lapply(1:ncol(gene_sig), function(i) {as.character(gene_sig[, i])} )
names(signatures) <- colnames(gene_sig)
signatures <- lapply(signatures, function(i) unique(i[!is.na(i)]))


u.scores <- ScoreSignatures_UCell(exp.mat, features = signatures)
write.csv(u.scores, "U_scores_neftel data.csv")
u.scores[1:8, 1]
dim(u.scores)

library(reshape2)
library(ggplot2)
melted <- melt(u.scores)
colnames(melted) <- c("Cell", "Signature", "Uscore")

p <- ggplot(melted, aes(x = Signature, y = Uscore)) + 
     geom_violin(aes(fill = Signature), scale = "width") + 
     geom_boxplot(width = 0.1) + theme_bw() + theme(axis.text.x = element_blank())
p

set.seed(123)
ranks <- StoreRankings_UCell(exp.mat)
ranks[1:5, 1:5]
# write.csv(ranks, "Ranks_neftel.csv")

seurat.object <- CreateSeuratObject(counts = exp.mat, project = "Neftel")

seurat.object <- AddMetaData(seurat.object, metadata = as.data.frame(u.scores))
head(seurat.object@meta.data)

seurat.object <- NormalizeData(seurat.object, verbose = FALSE)
seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 1000, 
                                      verbose = FALSE)

seurat.object <- ScaleData(seurat.object)
seurat.object <- RunPCA(seurat.object, features = seurat.object@assays$RNA@var.features, 
                        verbose = FALSE)
seurat.object <- RunUMAP(seurat.object, reduction = "pca", dims = 1:20, seed.use = 123, 
                         verbose = FALSE)
signature.names <- paste0(names(signatures), "_UCell")
FeaturePlot(seurat.object, reduction = "umap", 
            features = signature.names, ncol = 3, 
            order = T,
            cols = c("darkred", "orange", "yellow"))

q(save= "no")
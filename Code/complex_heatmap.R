library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(digest)
library(cluster)
library(readxl)
library(dplyr)

gene_list <- read_excel("Hallmarks gene.xlsx", sheet = 1)

mat <- read.csv('log2_tpm_t.csv')


### Hypoxia geneset selection
Hypoxia <- mat %>%filter(Genes %in% gene_list$Hypoxia)
row.names(Hypoxia) <- Hypoxia$Genes
Hypoxia <- Hypoxia[2:51]
Hypoxia <- t(scale(t(Hypoxia)))
x <- Hypoxia # Matrix for Hypoxia


### Glycolytic geneset selection
Glycolysis <- mat %>%filter(Genes %in% gene_list$Glycolysis)
row.names(Glycolysis) <- Glycolysis$Genes
Glycolysis <- Glycolysis[2:51]
Glycolysis <- t(scale(t(Glycolysis)))
x <- Glycolysis # Matrix for Glycolysis


### Mixed (Inflammatory response + TNFA signaling geneset selection)
Mixed <- mat %>% filter(Genes %in% gene_list$Mixed)
row.names(Mixed) <- Mixed$Genes
Mixed <- Mixed[2:51]
Mixed <- t(scale(t(Mixed)))
x <- Mixed # Matrix for Inflammatory response + TNFA signaling

annotation_core <- read.csv('Annotation_core.csv',row.names <- 1)
annotation_rim <- read.csv('Annotation_rim.csv',row.names <- 1)
annotation_inv <- read.csv('Annotation_inv.csv',row.names <- 1)
annotation_neg <- read.csv('Annotation_neg.csv',row.names <- 1)
annotation_pos <- read.csv('Annotation_pos.csv',row.names <- 1)

ann_colors <- list(
  Region <- c(tumor_core <- "#f4f4f4", 
             tumor_rim <- "#f9ab8c",
             invasive_region<-'#ffd66b',
             neg<-'#71cac7',
             pos<-'#e26240'))

core <- as.matrix(x[,1:10])
invasive <- as.matrix(x[,11:20])
negative <- as.matrix(x[,21:30])
positive <- as.matrix(x[,31:40])
rim <- as.matrix(x[,41:50])
Breaks <- c(-3,0,3)


h1 <- ComplexHeatmap::pheatmap(as.matrix(core), 
                              annotation_col = annotation_core,
                              annotation_colors = ann_colors,
                              annotation_names_row = F,
                              border_color = NA, breaks = Breaks,
                              border = TRUE, show_colnames = F)

h2 <- ComplexHeatmap::pheatmap(as.matrix(rim),
                              annotation_col = annotation_rim,
                              annotation_colors = ann_colors,
                              annotation_names_row = F,
                              border_color = NA, breaks = Breaks,
                              border = TRUE, show_colnames = F)

h3 <- ComplexHeatmap::pheatmap(as.matrix(invasive), 
                              annotation_col = annotation_inv,
                              annotation_colors = ann_colors,
                              annotation_names_row = F,
                              border_color = NA, breaks = Breaks,
                              border = TRUE, show_colnames = F)

h4 <- ComplexHeatmap::pheatmap(as.matrix(negative), 
                              annotation_col = annotation_neg,
                              annotation_colors = ann_colors,
                              annotation_names_row = F,
                              border_color = NA, breaks = Breaks,
                              border = TRUE, show_colnames = F)

h5 <- ComplexHeatmap::pheatmap(as.matrix(positive), 
                              annotation_col = annotation_pos,
                              annotation_colors = ann_colors,
                              annotation_names_col = F,
                              border_color = NA, breaks = Breaks,
                              border = TRUE, show_colnames = F)


### This will merge the 5 heatmaps in a single page
h1 + h2 + h3 + h4 + h5

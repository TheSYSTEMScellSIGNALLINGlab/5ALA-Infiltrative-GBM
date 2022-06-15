library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(digest)
library(cluster)

annotation_GPM = read.csv('S9I_Annotation_row_GPM.csv',row.names = 1)
annotation_Inflammatory = read.csv('S9I_Annotation_row_Inflammatory.csv',row.names = 1)
annotation_Injury = read.csv('S9I_Annotation_row_Injury.csv',row.names = 1)
annotation_Mesenchymal = read.csv('S9I_Annotation_row_Mesenchymal.csv',row.names = 1)
annotation_NFKB = read.csv('S9I_Annotation_row_NFKB.csv',row.names = 1)

column_ann = read.csv("S9I_Column_annotation.csv",row.names = 1)

ann_colors = list(
  Gene_set = c(GPM = "#F14668", 
             Inflammatory = "#f9ab8c",
             Inf_wound='#ffd66b',
             Mesenchymal='#71cac7',
             TNFA="#6ECB63"),
  Tumor_regions = c(Cellular_tumor = "#0E918C", 
                      Hyperplastic_blood_vessels = "#FF95C5",
                      Infiltrating_tumor='#A20A0A',
                      Leading_edge='#F0C929',
                      Microvascular_proliferation='#EFB7B7',
                      Perinecrotic_zone="#8FD6E1",
                      Pseudopalisading_cells="#EF8D32"))



mat = read.csv('S9I_Heatmap_matrix.csv',row.names = 1)
mat = t(scale(t(mat)))
colnames(mat)
Mesenchymal = as.matrix(mat[17:21,])
Injury = as.matrix(mat[5:16,])
GPM = as.matrix(t(mat[1:1,]))
Inflammatory = as.matrix(mat[2:4,])
NFKB = as.matrix(mat[22:29,])
row.names(GPM) = 'IL10RA'

# Breaks <- seq(min(mat), max(mat))
Breaks <- c(-4,4)

h1 = ComplexHeatmap::pheatmap(as.matrix(Mesenchymal), 
                              annotation_colors = ann_colors,
                              annotation_row = annotation_Mesenchymal,
                              annotation_col = column_ann,
                              cluster_cols = F,
                              show_colnames = F,
                              gaps_col = c(111,133,157,176,204,230),
                              border_color = NA, breaks = Breaks, border = TRUE)

h2 = ComplexHeatmap::pheatmap(as.matrix(Injury),
                              annotation_colors = ann_colors,
                              annotation_row = annotation_Injury,
                              cluster_cols = F,
                              # annotation_col = column_ann,
                              show_colnames = F,
                              annotation_names_row = F,
                              border_color = NA, breaks = Breaks, border = TRUE)

h3 = ComplexHeatmap::pheatmap(as.matrix(GPM), 
                              annotation_colors = ann_colors,
                              annotation_row = annotation_GPM,
                              cluster_cols = F,
                              # annotation_col = column_ann,
                              show_colnames = F,
                              annotation_names_row = F,
                              border_color = NA, breaks = Breaks, border = TRUE)

h4 = ComplexHeatmap::pheatmap(as.matrix(Inflammatory), 
                              annotation_colors = ann_colors,
                              annotation_row = annotation_Inflammatory,
                              cluster_cols = F,
                              # annotation_col = column_ann,
                              show_colnames = F,
                              annotation_names_row = F,
                              border_color = NA, breaks = Breaks, border = TRUE)

h5 = ComplexHeatmap::pheatmap(as.matrix(NFKB), 
                              annotation_colors = ann_colors,
                              annotation_row = annotation_NFKB,
                              cluster_cols = F,
                              # annotation_col = column_ann,
                              show_colnames = F,
                              annotation_names_col = F,
                              annotation_names_row = F,
                              border_color = NA, breaks = Breaks, border = TRUE)


h1 %v% h2 %v% h3 %v% h4 %v% h5
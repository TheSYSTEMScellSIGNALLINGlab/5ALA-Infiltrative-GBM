library(readxl)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(openxlsx)

# import matrix
mat = read.xlsx("1F.xlsx", rowNames = T)

mat[is.na(mat)] = 0
mat = t(mat)
pheatmap(mat,  show_colnames = T, 
         fontsize_row = 10, fontsize_col = 10,
         angle_col = 90,
         cluster_cols = F,
         cluster_rows = F,
         clustering_method = "complete",
         color = inferno(100),
         cellwidth = 10, cellheight = 10, border_color = 'grey')


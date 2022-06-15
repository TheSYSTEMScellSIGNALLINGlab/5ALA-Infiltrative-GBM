library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(digest)
library(cluster)

annotation_ala = read.csv('Annotation_ala_5N.csv',row.names = 1)
annotation_gp = read.csv('Annotation_gp_5N.csv',row.names = 1)
annotation_olc = read.csv('Annotation_olc_5N.csv',row.names = 1)
annotation_others = read.csv('Annotation_others_5N.csv',row.names = 1)

ann_colors = list(
  Region = c(ALA = "#f4f4f4", 
             GP = "#f9ab8c",
             OLC ='#ffd66b',
             Others ='#71cac7'))

mat = read.csv('5N.csv',row.names = 1)
mat = log2(mat+1)
mat = t(scale(t(mat)))
dim(mat)

ala = as.matrix(mat[,1:104]) #104
gp = as.matrix(mat[,105:296]) #192
olc = as.matrix(mat[,297:431]) #135
others = as.matrix(mat[,432:600]) #169


Breaks <- c(-3,3)

h1 = ComplexHeatmap::pheatmap(as.matrix(ala), 
                              annotation_col = annotation_ala,
                              annotation_colors = ann_colors,
                              annotation_names_row = F,
                              show_colnames = F,
                              border_color = NA, breaks = Breaks,
                              border = TRUE)

h2 = ComplexHeatmap::pheatmap(as.matrix(gp),
                              annotation_col = annotation_gp,
                              annotation_colors = ann_colors,
                              annotation_names_row = F,
                              show_colnames = F,
                              border_color = NA, breaks = Breaks,
                              border = TRUE)

h3 = ComplexHeatmap::pheatmap(as.matrix(olc), 
                              annotation_col = annotation_olc,
                              annotation_colors = ann_colors,
                              annotation_names_row = F,
                              show_colnames = F,
                              border_color = NA, breaks = Breaks,
                              border = TRUE)

h4 = ComplexHeatmap::pheatmap(as.matrix(others),
                              annotation_col = annotation_others,
                              annotation_colors = ann_colors,
                              annotation_names_row = F,
                              show_colnames = F,
                              border_color = NA, breaks = Breaks,
                              border = TRUE)


h1 + h2 + h3 + h4 




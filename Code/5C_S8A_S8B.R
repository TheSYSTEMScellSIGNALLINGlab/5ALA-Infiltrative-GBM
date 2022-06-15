library(readxl)
library(pheatmap)
library(RColorBrewer)

#load matrix for heatmap
mat = read.csv(file.choose(), row.names = 1) # 5C / S8A / S8B

#load annotation column for heatmap
annotation_col = read.csv(file.choose(), row.names = 1) # 5C annotation column / S8A annotation column / S8B annotation column

#set annotation colors according to the annotation files
ann_colors = list(
  Marker = c(Yamanaka_factor = "#c6dc6a", 
             Cancer_sc = "#48c3cf",
             Normal_sc='#f47541')
)

#generate heatmap
pheatmap(mat,  show_colnames = T, 
         annotation_col = annotation_col,
         annotation_colors = ann_colors, 
         fontsize_row = 10, fontsize_col = 10,
         angle_col = 90,
         cluster_cols = F,
         cluster_rows = F,
         cellwidth = 10, cellheight = 10, border_color = 'black',
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
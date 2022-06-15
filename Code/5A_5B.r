library(pheatmap)
library(readxl)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(openxlsx)


# Figure 5A
mat = read.xlsx("5A.xlsx",rowNames = T)

ComplexHeatmap::pheatmap(mat,  show_colnames = T,cluster_rows = T,
         fontsize_row = 10, fontsize_col = 10,
         clustering_distance_cols = 'spearman', 
         clustering_distance_rows = 'spearman',
         cellwidth = 10, cellheight = 10, border_color = 'black',
         color = colorRampPalette(c("#f2ea0f","#ee2829"))(100))


# Figure 5B
jitter = read.csv('5B.csv')
attach(jitter)
names(jitter)

level_order = c('Core', 'Rim', 'Inv',	'5ALA-', '5ALA+', 'Classical',
                'Mesenchymal', 'Neural', 'Proneural')

ggplot(jitter, aes(x = factor(Condition, level = level_order),
                   y = StemnessScore, fill = Condition))+
  stat_summary(fun.y='mean', geom='bar', colour="black") +
  scale_fill_manual(values = c("#71cac7","#e46141","#d3942b","#f6f8f7",
                               "#fed66a","#93ab3c",
                               "#31b34a","#ea0d8d","#f9aa8c"))+
  geom_jitter()+
  theme_bw() + facet_grid(~Source, scales = "free", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90))+
  xlab("") + theme(legend.position = "none")

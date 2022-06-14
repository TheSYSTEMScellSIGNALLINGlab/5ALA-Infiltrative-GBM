library(readxl)
library(pheatmap)
library(openxlsx)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(scales)

# Figure 4A
data = read_excel("4A.xlsx")

attach(data)


ggplot(data, aes(X,Y)) +                     
  geom_point(aes(color = U_cell), size = 1) +
  theme_bw()+ 
  scale_colour_gradientn(colours = c("darkred", "orange", "yellow", "white"))+
  ggtitle("Neftel et al. dataset") +
  xlab("Relative meta-module score") + ylab("Relative meta-module score")

  
 
# Figure 4B
data = read.csv("4B.csv")

ggplot(data, aes(X,Y)) +                     
  geom_point(aes(color = ALA_sig_UCell), size = 1) +
  theme_bw()+ scale_colour_gradientn(colours = c("darkred", "orange", "yellow", "white")) + 
  xlab('pc1') + ylab("pc2")+
  ggtitle("Richards et al. dataset")
setwd("H:/lab data/Glioblastoma/Figures/Supplementary figure 9")

library(readxl)
data = read_excel("Sup 9.xlsx")
library(ggplot2)
library(ggthemes)
attach(data)
library(dplyr)
library(hrbrthemes)
library(viridis)
data$pathway = tolower(pathway)

ggplot(data, aes(x = -log10(padj), y = reorder(pathway,-log10(padj)), 
                 size = size,
                 fill = NES)) +
  geom_point(shape=21, color = "Black") + theme_bw() + 
  scale_color_viridis(option = "D")+
  facet_grid(.~Category, scales = "fixed") +
  ylab("Hallmarks") +
  xlab("-log10 (adj. p-value)")

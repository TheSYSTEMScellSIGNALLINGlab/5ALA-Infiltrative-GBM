library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(corrr)


### Reading the annotation file
Annot <- read_excel("TCGA annotation wo normal.xlsx")


### Reading the cibersortx dataframes
Neftel <- read_excel("Neftel_cibersortx.xlsx")
Garofano <- read_excel("Garofano_cibersortx.xlsx")
Richards <- read_excel("Richards_cibersortx.xlsx")


######################################## Part 1 #######################################


### Filtering out only the paired GBM-TCGA samples
Neftel <- Neftel %>% filter(Mixture %in% Annot$Samples) 
Garofano <- Garofano %>% filter(Mixture %in% Annot$Samples) 
# Injury <- Injury %>% filter(Mixture %in% Annot$Samples) 
# Dev <- Dev %>% filter(Mixture %in% Annot$Samples) 
Richards <- Richards %>% filter(Mixture %in% Annot$Samples) 


### Shaping the dataframe

Neftel_re <- Neftel[1:7]
Neftel_re <- Neftel_re %>% distinct(Mixture, .keep_all = TRUE)
Neftel_re <- melt(Neftel_re, id.vars=1)
Neftel_re <- as.matrix(Neftel_re)
Neftel_re <- as_tibble(Neftel_re)
Neftel_re
Neftel_re$value <- as.double(Neftel_re$value)
Neftel_merged <- merge(Neftel_re, Annot, by.x = "Mixture", by.y = "Samples",  all.x = T)


Garofano_re <- Garofano[1:5]
Garofano_re <- Garofano_re %>% distinct(Mixture, .keep_all = TRUE)
Garofano_re <- melt(Garofano_re, id.vars=1)
Garofano_re <- as.matrix(Garofano_re)
Garofano_re <- as_tibble(Garofano_re)
Garofano_re
Garofano_re$value <- as.double(Garofano_re$value)
Garofano_merged <- merge(Garofano_re, Annot, by.x = "Mixture", by.y = "Samples",  all.x = T)


Richards_re <- Richards[1:3]
Richards_re <- Richards_re %>% distinct(Mixture, .keep_all = TRUE)
Richards_re <- melt(Richards_re, id.vars=1)
Richards_re <- as.matrix(Richards_re)
Richards_re <- as_tibble(Richards_re)
Richards_re
Richards_re$value <- as.double(Richards_re$value)
Richards_merged <- merge(Richards_re, Annot, by.x = "Mixture", by.y = "Samples",  all.x = T)

### ggplot creation
library(ggplot2)

neftel_plot <-  ggplot(Neftel_merged, aes(x = reorder(Mixture, value), y = value, fill = variable)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("#17A19E","#FCB247","#E86D76","#22577A",
                               "#97D3BE","#4590CE"))+
  theme_bw() + facet_grid(~Group, scales = "free", space = "free_x") +
  # theme(axis.text.x = element_text(angle = 90))+
  xlab("") + ylab("Fraction")+ ggtitle("Neftel et al")

garofano_plot <-  ggplot(Garofano_merged, aes(x = reorder(Mixture, value), y = value, fill = variable)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("#CCCE79","#2EB086","#337D3D","#21325E"))+
  theme_bw() + facet_grid(~Group, scales = "free",  space = "free_x") +
  # theme(axis.text.x = element_text(angle = 90))+
  xlab("") + ylab("Fraction") + ggtitle("Garofano et al")

richards_plot <-  ggplot(Richards_merged, aes(x = reorder(Mixture, value), y = value, fill = variable)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("#FAB8A1","#B97A95"))+
  theme_bw() + facet_grid(~Group, scales = "free",  space = "free_x") +
  # theme(axis.text.x = element_text(angle = 90))+
  xlab("") + ylab("Fraction") + ggtitle("Richards et al")

require(gridExtra)

# Figure S12A, S12B, S12C

Plot <- grid.arrange(neftel_plot, garofano_plot, richards_plot, ncol=1)


ggsave(Plot, file = file.path("H:/lab data/Glioblastoma/cibersortx/TCGA/TCGA_Cibersotx.pdf"),
       width = 8.27, height = 6)


#################################### End of Part 1 ####################################


####################################### Part 2 ######################################## 

### Filtering out only the paired GBM-TCGA samples
Neftel <- Neftel %>% filter(Mixture %in% Annot$Samples) 
Garofano <- Garofano %>% filter(Mixture %in% Annot$Samples) 
Richards <- Richards %>% filter(Mixture %in% Annot$Samples) 


### Shaping the dataframe

Neftel_re <- Neftel[1:7]
Neftel_re <- melt(Neftel_re, id.vars=1)

Garofano_re <- Garofano[1:5]
Garofano_re <- melt(Garofano_re, id.vars=1)

Richards_re <- Richards[1:3]
Richards_re <- melt(Richards_re, id.vars=1)


TCGA <- rbind(Neftel_re, Garofano_re, Richards_re)


### Adding dataset annotation by dplyr-mutate
TCGA <- TCGA %>% mutate(Datasets = ifelse(variable == "AClike", "Neftel", 
                                            ifelse(variable == "MESlike1", "Neftel",
                                                   ifelse(variable == "MESlike2", "Neftel",
                                                          ifelse(variable == "NPClike1", "Neftel",
                                                                 ifelse(variable == "NPClike2", "Neftel",
                                                                        ifelse(variable == "OPClike", "Neftel",
                                                                               ifelse(variable == "GPM", "Garofano",
                                                                                      ifelse(variable == "MTC", "Garofano",
                                                                                             ifelse(variable == "PPR", "Garofano",
                                                                                                    ifelse(variable == "NEU", "Garofano",
                                                                                                           ifelse(variable == "Injury", "Richards",
                                                                                                                  ifelse(variable == "Dev", "Richards", "")))))))))))))


unique(TCGA$variable)

TCGA_merged <- merge(TCGA, Annot, by.x = "Mixture", by.y = "Samples",  all.x = T)
attach(TCGA_merged)
write.csv(TCGA_merged, "TCGA_ggplot.csv")

### ggplot creation

x <- read.csv("TCGA_ggplot.csv")
x <- TCGA_merged
p <- ggplot(x, aes(x = Group, y = value, fill = variable)) +
  geom_bar(position="fill", stat="identity") +
  # scale_fill_manual(values = c("#71cac7","#e46141","#d3942b","#f6f8f7",
  #                              "#fed66a","#93ab3c"))+
  theme_bw() + facet_grid(~Datasets, scales = "free") +
  theme(axis.text.x = element_text(angle = 90))+
  xlab("") + ylab("Fraction")

# Figure 7A

p

ggsave(p, file = file.path("H:/lab data/Glioblastoma/cibersortx/TCGA/TCGA_Cibersotx_mean.pdf"),
       width = 4, height = 4)


#################################### End of Part 2 ####################################
library(psych)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggthemes)
library(dplyr)


##boxplot 1A
log_tpm <- read.delim('Sup 1A.txt', sep = ",", row.names = 1)

box <- log_tpm %>% select(names(log_tpm)) %>%
  pivot_longer(., cols = c(colnames(log_tpm)), 
               names_to = "Var", values_to = "Val")

x <- box %>% mutate(Category= ifelse(grepl("Core", Var), "a_Core",
                        ifelse(grepl("Rim", Var), "b_Rim",
                               ifelse(grepl("Inv", Var), "c_Inv",
                                      ifelse(grepl("Neg", Var), "Neg",
                                             ifelse(grepl("Pos", Var), "Pos","Other"))))))

ggplot(x, aes(x = Var, y = Val, fill = Category)) +
geom_boxplot(size = .1) + theme_light()+ facet_grid(~Category, scales = "free")+
  theme(axis.text.x = element_text(angle = 90))



##correlation plot 1C
cor_data <- read.delim('H:\\lab data\\Glioblastoma\\Github\\5ALA-Infiltrative-GBM\\Data\\Sup 1C.txt', sep=",")
attach(cor_data)
pairs(cor_data[2:6], pch = 19, lower.panel = NULL)

pairs.panels(cor_data[,2:6], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE, # show correlation ellipses
             col="#69b3a2"
)
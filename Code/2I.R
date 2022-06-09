library(readxl)
data = read_excel("2I.xlsx")
attach(data)

level_order = c('Core', 'Rim', 'Inv',	'5ALA neg', '5ALA pos')

library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggplot2)
library(readxl)

# Stacked bar plot

ggplot(data, aes(fill=Categories, y=Values, x=factor(Group, level = level_order))) + 
  geom_bar(position="fill", stat="identity") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Fraction") + xlab("")

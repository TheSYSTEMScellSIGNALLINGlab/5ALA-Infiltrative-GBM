library(readxl)
library(ggthemes)
library(ggpubr)
library(ggplot2)
library(readxl)

# Stacked + percent
# Figure 3F
data = read_excel("3F.xlsx")
attach(data)

level_order = c('Core', 'Rim', 'Inv',	'5ALA-', '5ALA+')

ggplot(data, aes(fill=Category, y=Values, x=factor(Group, level = level_order))) + 
  geom_bar(position="fill", stat="identity") + 
  facet_grid(~Pathway, scales = "free", space = "free_x")+
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("fraction") + xlab("")


# All patient bar jitter
# Figure 3G and 3H
data = read_excel("3G_3H.xlsx", sheet = 1) # 3G
data = read_excel("3G_3H.xlsx", sheet = 2) # 3H

ggplot(data, aes(x = Category, y = Value, fill = Category))+
  stat_summary(fun.y='mean', geom='bar', colour="black") +
  scale_fill_manual(values = c("#71cac7","#e46141","#d3942b"))+
  geom_jitter()+
  theme_bw() + facet_grid(~Group, scales = "free", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90))+
  xlab("") + ylab("Fraction")



# Injury_high/Dev_high
# Figure 3I
Injury_dev = read_excel("3I.xlsx")


ggplot(Injury_dev, aes(Mixture, Logratio))+
  geom_line(aes(group = Pair, size =2))+
  geom_point(aes(color = Group, size =4, fill = Group))+
  theme_bw()+
  # scale_shape_identity()+ facet_grid(~Group, scales = "free", space = "free_x")+
  theme(axis.text.x = element_text(angle = 90)) + xlab("") + 
  ylab("Log 10 (Injury-High/ Developmental-High)")



library(dplyr)
library(tidyr)
library(tidyverse)

# Figure S3L
library(readxl)
data = read_excel("3L.xlsx")
attach(data)

ggplot(data, aes(fill=Categories, y=Values, x=Group)) + 
  geom_bar(position="fill", stat="identity") + 
  facet_grid(~Patients, scales = "free", space = "free_x")+
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("fraction") + xlab("")


# Figure S4E and S4F
library(readxl)
data = read_excel("S4E.xlsx") # S4E
data = read_excel("S4F.xlsx") # S4F
attach(data)

level_order = c('Core', 'Rim', 'Inv',	'5ALA-', '5ALA+')

library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggplot2)
library(readxl)

# Stacked + percent

ggplot(data, aes(fill=Group, y=Value, x=factor(Category, levels = level_order))) + 
  geom_bar(position="fill", stat="identity") + 
  facet_grid(~Mixture, scales = "free", space = "free_x")+
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("fraction") + xlab("")


# Dumbell
# Figure S4G
dum = read_excel("S4G.xlsx")

ggplot(dum, aes(Mixture, Value))+
  geom_line(aes(group = Pair, size =2))+
  geom_point(aes(color = Category, size =4, fill = Group))+
  theme_bw()+
  scale_shape_identity()+ facet_grid(~Group, scales = "free", space = "free_x")+
  theme(axis.text.x = element_text(angle = 90)) + xlab("") + 
  ylab("Fraction")

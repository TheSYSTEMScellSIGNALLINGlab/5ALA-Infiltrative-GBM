library(readxl)
library(ggpubr)
library(ggplot2)
library(ggthemes)

data = read_excel('Sup 6.xlsx')   #data loading
attach(data)
uniq_genes = unique(data$Genes)
level_order <- c('5ALA+ - Core', '5ALA+ - Rim', '5ALA+ - Inv',	'5ALA+ - 5ALA-')

for(i in uniq_genes) {                              
  p = ggplot(data = subset(data, Genes ==i))+
    geom_bar(stat='identity', color ='black',
             position=position_dodge(),
             aes(x = factor(Region, level = level_order), 
                 y = Delta, fill = Category))+ggtitle(i)+
    scale_fill_manual(values = c('#f9c838','#e77e39','#79aa40'))+
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    xlab("")
  
  ggsave(p, file=paste0(i,".pdf"), width = 4, height = 3.5, units = "in")
}




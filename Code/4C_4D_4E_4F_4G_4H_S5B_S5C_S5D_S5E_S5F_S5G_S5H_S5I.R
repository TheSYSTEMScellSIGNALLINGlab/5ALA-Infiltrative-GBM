mat = read.csv("4C_4D_4E_4F_4G_4H_S5B_S5C_S5D_S5E_S5F_S5G_S5H_S5I.csv")


#################################################
#################################################
# Figure 4C, 4D, 4E, 4F, 4G, 4H, S5D, S5E, S5F, S5G, S5H, S5I

getDensity <- function(x, y)
{
  dc <- densCols(x, y, colramp=colorRampPalette(c("black", "white")))
  dens <- col2rgb(dc)[1,] + 1L
  dens <- dens/max(dens)
  return(dens)
}

A = getDensity(mat$ALA_sig_UCell, mat$GPM_UCell)
B = getDensity(mat$ALA_sig_UCell, mat$MTC_UCell)
C = getDensity(mat$ALA_sig_UCell, mat$NEU_UCell)
D = getDensity(mat$ALA_sig_UCell, mat$PPR_UCell)
E = getDensity(mat$ALA_sig_UCell, mat$MES2_UCell)
F = getDensity(mat$ALA_sig_UCell, mat$MES1_UCell)
G = getDensity(mat$ALA_sig_UCell, mat$AC_UCell)
H = getDensity(mat$ALA_sig_UCell, mat$OPC_UCell)
I = getDensity(mat$ALA_sig_UCell, mat$NPC1_UCell)
J = getDensity(mat$ALA_sig_UCell, mat$NPC2_UCell)
K = getDensity(mat$ALA_sig_UCell, mat$Injury_UCell)
L = getDensity(mat$ALA_sig_UCell, mat$Developmental_UCell)
M = getDensity(mat$MES2_UCell, mat$Injury_UCell)
N = getDensity(mat$MES1_UCell, mat$Injury_UCell)


library(ggplot2)
library(ggthemes)

myCol <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k)) 
 
c1 = cor(mat$ALA_sig_UCell, mat$GPM_UCell, method = "pearson") # Calculation of pearson correlation
c1 = specify_decimal(c1, 3)
c1 = as.character(c1)
c1 = paste("p = ", c1)

cor.test(mat$ALA_sig_UCell, mat$Injury_UCell, method = "pearson")

p = ggplot(mat, aes(ALA_sig_UCell, GPM_UCell, col = A)) +
    geom_point(size = 1, alpha=0.5) +
    scale_colour_gradientn(limits=c(0, 1), colours = myCol) +
    theme_bw() + 
    annotate("text", x = 0.2, y=0.01, label = c1)

#################################################
#################################################

library(reshape2)

x = mat
dim(x)

Neftel = x[, c(1:2,4, 9:14)]
Garofano = x[, c(1,3:8)]

Neftel_GG <- melt(Neftel, id.vars=1:2)
Garofano_GG <- melt(Garofano, id.vars = 1:2)

Neftel_factor = c("ALA_sig_UCell","AC_UCell", "MES1_UCell", "MES2_UCell", "NPC1_UCell", "NPC2_UCell", "OPC_UCell")


library(ggplot2)
attach(x)

# Figure S5B
ggplot(Neftel_GG, aes(x = factor(variable, level = Neftel_factor), y =value, fill = variable)) +
  geom_violin() +
  theme_bw() + facet_wrap(~Neftel,  ncol = 1) +
  theme(axis.text.x = element_text(angle = 90))+
  xlab("") + ylab("Ucell score")

# Figure S5C
ggplot(Garofano_GG, aes(x = variable, y =value, fill = variable)) +
  geom_violin() +
  theme_bw() + facet_wrap(~Garofano,  ncol = 1) +
  theme(axis.text.x = element_text(angle = 90))+
  xlab("") + ylab("Ucell score")

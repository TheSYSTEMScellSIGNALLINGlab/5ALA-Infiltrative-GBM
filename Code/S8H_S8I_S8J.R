mat = read.csv("S8H_S8I_S8J.csv")

getDensity <- function(x, y)
{
  dc <- densCols(x, y, colramp=colorRampPalette(c("black", "white")))
  dens <- col2rgb(dc)[1,] + 1L
  dens <- dens/max(dens)
  return(dens)
}

A = getDensity(mat$ALA_sig_UCell, mat$GP_UCell)
B = getDensity(mat$ALA_sig_UCell, mat$OLC_UCell)
C = getDensity(mat$ALA_sig_UCell, mat$Others_UCell)


library(ggplot2)
library(ggthemes)

myCol <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k)) 
 
c1 = cor(mat$ALA_sig_UCell, mat$GP_UCell, method = "pearson")
c1 = specify_decimal(c1, 3)
c1 = as.character(c1)
c1 = paste("p = ", c1)
cor.test(mat$ALA_sig_UCell, mat$Others_UCell, method = "pearson")
p1 = ggplot(mat, aes(ALA_sig_UCell, GP_UCell, col = A)) +
    geom_point(size = 1, alpha=0.5) +
    scale_colour_gradientn(limits=c(0, 1), colours = myCol) +
    theme_bw() + annotate("text", x = 0.2, y=0.01, label = c1)

c1 = cor(mat$ALA_sig_UCell, mat$OLC_UCell, method = "pearson")
c1 = specify_decimal(c1, 3)
c1 = as.character(c1)
c1 = paste("p = ", c1)

p2 = ggplot(mat, aes(ALA_sig_UCell, OLC_UCell, col = B)) +
  geom_point(size = 1, alpha=0.5) +
  scale_colour_gradientn(limits=c(0, 1), colours = myCol) +
  theme_bw() + annotate("text", x = 0.19, y=0.04, label = c1)

c1 = cor(mat$ALA_sig_UCell, mat$Others_UCell, method = "pearson")
c1 = specify_decimal(c1, 3)
c1 = as.character(c1)
c1 = paste("p = ", c1)

p3 = ggplot(mat, aes(ALA_sig_UCell, Others_UCell, col = C)) +
  geom_point(size = 1, alpha=0.5) +
  scale_colour_gradientn(limits=c(0, 1), colours = myCol) +
  theme_bw() + annotate("text", x = 0.19, y=0.01, label = c1)


require(gridExtra)
grid.arrange(p1, p2, p3, ncol=3)

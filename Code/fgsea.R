library(fgsea)
library(data.table)


#Rank file defining
rnk.file = system.file("extdata", "Rank_file.rnk", package="fgsea")

#Gene set file defining
gmt.file <- system.file("extdata", "Gene_sets.gmt", package="fgsea")

ranks <- read.table(rnk.file,
                    header=TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$t, ranks$ID)
str(ranks)
pathways <- gmtPathways(gmt.file)
str(head(pathways))
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=1000)
head(fgseaRes)

collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.05], 
                                      pathways, ranks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathways[mainPathways], ranks, fgseaRes, 
              gseaParam = 0.5)


#plot the enriched pathways by defining the pathway names from the Gene set file
plotEnrichment(pathways[["VERHAAK_GLIOBLASTOMA_CLASSICAL"]],ranks)

#write the fgsea result in a table
fwrite(fgseaRes, file="pos_rim.txt", sep="\t", sep2=c("", " ", ""))

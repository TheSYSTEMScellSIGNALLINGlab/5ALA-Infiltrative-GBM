
library(limma)
library(edgeR)
library(openxlsx)
library(UpSetR)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################


entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Hs.egSYMBOL, ifnotfound=NA)
  symbol <- unlist(lapply(symbol, function(i) return(i[1])))
  return(symbol)
}

ensembl2entrez <- function(ensembl)
{
  entrez <- mget(as.character(ensembl), org.Hs.egENSEMBL2EG, ifnotfound=NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

symbol2entrez <- function(symbol)
{
  entrez <- mget(as.character(symbol), org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrez <- unlist(lapply(entrez, function(i) return(i[1])))
  entrez <- unique(entrez[!is.na(entrez)])
  return(entrez)
}

entrez2name <- function(entrez)
{
  gn <- mget(as.character(entrez), org.Hs.egGENENAME, ifnotfound=NA)
  gn <- unlist(lapply(gn, function(i) return(i[1])))
  return(gn)
}


getDEG <- function(limma, pv, fc, dsign = "both", adjusted = TRUE)
{
  limmaFC <- limma$logFC
  if(class(limmaFC)=="factor") limmaFC <- as.numeric(levels(limmaFC))[limmaFC]
  ifelse(adjusted, limmaPV <- limma$adj.P.Val, limmaPV <- limma$P.Value)
  if(class(limmaPV)=="factor") limmaPV <- as.numeric(levels(limmaPV))[limmaPV]
  
  if(dsign == "up") return(as.character(limma$entrez)[limmaFC > fc & limmaPV < pv])
  else if(dsign == "down") return(as.character(limma$entrez)[limmaFC < -fc & limmaPV < pv])
  else if(dsign == "both") return(as.character(limma$entrez)[abs(limmaFC) > fc & limmaPV < pv])
}

getNbDEG <- function(limma, pv, fc, adjusted = FALSE)
{
  limmaFC <- limma$logFC
  if(class(limmaFC)=="factor") limmaFC <- as.numeric(levels(limmaFC))[limmaFC]
  ifelse(adjusted, limmaPV <- limma$adj.P.Val, limmaPV <- limma$P.Value)
  if(class(limmaPV)=="factor") limmaPV <- as.numeric(levels(limmaPV))[limmaPV]
  
  return(sum(abs(limmaFC)>= fc & limmaPV <= pv, na.rm = TRUE))
}



ggVolcanoV2 <- function(limma, pvalue = 0, genes = NULL, outFile){
  ggmat <- data.frame(X = limma$logFC, Y = -log10(limma$adj.P.Val), GENE = limma$symbol, AVG = limma$AveExpr)
  ggmat$isSignif <- "no"
  ggmat$isSignif[limma$adj.P.Val < pvalue] <- "yes"
  ggmat$toShow <- "no"
  if(!is.null(genes)){
    ggmat$toShow <- "no"
    ggmat$toShow[match(genes, ggmat$GENE)] <- "yes"
    
    ggmat$isSignif <- "no"
    ggmat$isSignif[match(genes, ggmat$GENE)] <- "yes"
  }
  
  p <- ggplot(ggmat, aes(x=X, y=Y, color = AVG, label = GENE))
  p <- p + geom_point(size = 2, alpha = 0.5) + scale_color_viridis_c()
  p <- p + geom_point(size = 2, alpha = 0.5, shape = 21, color = "black", data = subset(ggmat, toShow == "yes"))
  if(pvalue != 0) p <- p + geom_hline(yintercept=-log10(pvalue), linetype="dashed", color = "black")
  p <- p + theme_bw(base_size = 20)
  p <- p + theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18))
  p <- p + xlab("log2 Fold Change")+ ylab("-log10 Pvalue")
  #p <- p + theme(legend.position="none")
  p <- p + geom_text_repel(data = subset(ggmat, toShow == "yes"), color = "black", box.padding = 0.5, max.overlaps = Inf)
  
  pdf(outFile)
  plot(p)
  dev.off()
}

getNbMat <- function(limma, pvs, fcs, sign="ALL"){
  nbMat <- c()
  for(i in pvs){
    for(j in fcs){
      if(sign == "ALL"){
        nb <- sum(limma$adj.P.Val < i & abs(limma$logFC) > j)
      }else if(sign == "UP"){
        nb <- sum(limma$adj.P.Val < i & limma$logFC > j)
      }else if(sign == "DOWN"){
        nb <- sum(limma$adj.P.Val < i & limma$logFC < -j)
      }else(return(NA))
      nbMat <- rbind(nbMat, c(i, j, nb))
    }
  }
  colnames(nbMat) <- c("Pvalue", "logFC", "NB")
  
  nbMat <- as.data.frame(nbMat)
  nbMat$Pvalue <- as.character(nbMat$Pvalue)
  nbMat$Sign <- sign
  return(nbMat)
}





############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

# PARAMETERS
countDir <- file.path("../Data/")

limmaDir <- file.path("../limma_voom/")
dir.create(limmaDir, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES
exonCount <- "Final_count_exon.csv"
sampleAnnotation <- "Annotation_col.csv"

##############
# LOAD SAMPLES

setwd(countDir)
cntEx <- read.csv(exonCount)
rownames(cntEx) <- cntEx[,1]
cntEx <- cntEx[, -1]
Rex <- cntEx[, colnames(cntEx) != "width"]
countMat <- Rex


# Create geneMat
entrezID <- rownames(countMat)
symbol <- entrez2symbol(entrezID)
geneName <- entrez2name(entrezID)

geneMat <- data.frame(entrez = entrezID, symbol = symbol, gene.name = geneName)

#rownames(countMat) <- ensID

# Min row CPM >= 2
#logcpmMat <- cpm(countMat, log = TRUE)
#keep.exprs <- rowSums(logcpmMat > 0) >= 2

# count != 0 in at least 2 samples
nonZero <- apply(countMat, 1, function(i) sum(i != 0))
keep.exprs <- nonZero >= 5

# SAMPLE ANNOTATION
ann.sample <- read.csv(sampleAnnotation)
ann.sample$SAMPLE <- gsub(" ", ".", ann.sample$id)
ann.sample <- ann.sample[match(colnames(countMat), ann.sample$SAMPLE), ]
ann.sample$PAIR <- sapply(strsplit(ann.sample$id, split = "_"), function(i) i[1])
ann.sample$PAIR <- gsub(" ", "", ann.sample$PAIR)

# BUILD DGE
dge <- DGEList(count = countMat, genes = geneMat)
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge, method = "TMM") # normalization

m <- dge$counts
rownames(m) <- dge$genes$ensembl


# SAVE COUNT
#setwd(countDir)
#write.table(m, "count.txt", sep = "\t", quote = FALSE)
#write.table(cpm(m, log = TRUE), "logCPM.txt", sep = "\t", quote = FALSE)



###################################################################################
# LIMMA 

# PLEASE CHECK "groups", "design" and "contrast.matrix" variables

#######################
# Differential analysis
groups <- ann.sample$Annotation
pairs <- ann.sample$PAIR

design <- model.matrix(~0+groups+pairs) # paired
colnames(design) <- gsub("groups", "", colnames(design))
v <- voom(dge, design)

# UNPAIRED
fit <- lmFit(v, design, method = "robust")# robust or not 

fit <- eBayes(fit)
plotSA(fit)

tt <- topTable(fit, adjust="BH", number=nrow(dge), p.value = 1)

contrast.matrix <- makeContrasts(
  Pos-Core,
  Pos-Rim,
  Pos-Inv,
  Pos-Neg,
  levels=design
)
contrast.name <- gsub(" ", "", colnames(contrast.matrix))

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)	

resList <- lapply(1:ncol(contrast.matrix), function(i)
  topTable(fit2, coef=i, adjust="BH", number=nrow(dge), sort.by="P", p.value = 1)
)
names(resList) <- contrast.name

# SAVE results
setwd(limmaDir)
lapply(1:length(resList), function(i)
  write.xlsx(resList[[i]], paste(names(resList)[i], "_limma.xlsx", sep = ""),
             firstRow = T, headerStyle = createStyle(textDecoration = 'bold'))
)



###################################################################################
# VOLCANO PLOTS


# Load limma outputs
setwd(limmaDir)
limmaFiles <- list.files(pattern = "_limma.xlsx")
names(limmaFiles) <- gsub("_limma.xlsx", "", limmaFiles)
limmaList <- lapply(limmaFiles, read.xlsx, sheet = 1)

limmaList <- lapply(limmaList, function(i) i[!is.na(i$symbol), ])

lapply(1:length(limmaList), function(i)
  ggVolcanoV2(limmaList[[i]], pvalue = 0.05, genes = head(limmaList[[i]]$symbol, 25),
              outFile = paste0(names(limmaList)[i], "_volcano.pdf"))
)



###################################################################################
# NUMBER OF DEG PLOT

setwd(limmaDir)
limmaFiles <- list.files(pattern = "_limma.xlsx")
names(limmaFiles) <- gsub("_limma.xlsx", "", limmaFiles)
limmaList <- lapply(limmaFiles, read.xlsx, sheet = 1)

pvs <- c(0.1, 0.05, 0.01, 0.005, 0.001)
fcs <- seq(0,3, by = 0.5)

for(i in 1:length(limmaList)){
  limmaMat <- limmaList[[i]]
  
  # plot
  ggmat.ALL <- getNbMat(limmaMat, pvs, fcs, sign="ALL")
  ggmat.UP <- getNbMat(limmaMat, pvs, fcs, sign="UP")
  ggmat.DOWN <- getNbMat(limmaMat, pvs, fcs, sign="DOWN")
  
  ggmat <- do.call(rbind,list(ggmat.ALL, ggmat.UP, ggmat.DOWN))
  
  p <- ggplot(data=ggmat, aes(x=logFC, y=NB, group=Pvalue, colour=Pvalue))
  p <- p + theme_bw()
  p <- p + geom_line()
  p <- p + geom_point()
  p <- p + facet_grid(. ~ Sign)
  #p <- p + scale_y_log10()
  p <- p + ylab("Number of DEG")
  pdf(paste0(names(limmaList)[i], "_nb_diff.pdf"), width = 10, height = 5)
  plot(p)
  dev.off()
  
}



library(QuasR)
library(eisaR)
library(edgeR)
library(openxlsx)
library(org.Hs.eg.db)
library(pheatmap)
library(ggplot2)
library(corrplot)
library(RColorBrewer)
library(tidyr)
library(fgsea)
library(viridis)

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################


my.read_xlsx <- function(inFile)
{
  mysheets <- getSheetNames(inFile)
  mList <- lapply(mysheets, read.xlsx, xlsxFile = inFile)
  names(mList) <- mysheets
  return(mList)
}


entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Hs.egSYMBOL, ifnotfound=NA)
  symbol <- unlist(lapply(symbol, function(i) return(i[1])))
  return(symbol)
}

std <- function(x) sd(x)/sqrt(length(x))

getEmpiricalPV <- function(x, y, type)
{
  x.ecdf <- stats::ecdf(x)
  if(type == "greater") return(1-x.ecdf(y))
  else if(type == "less") return(x.ecdf(y))
  return(NA)    
}

getHeatmapMatrix <- function(gsea_list, scoreName)
{
  trms <- sort(unique(unlist(lapply(gsea_list, function(i) i[,1]))))
  m <- matrix(NA, nrow = length(trms), ncol = length(gsea_list))
  
  for(i in 1:length(gsea_list))
  {
    gs <- gsea_list[[i]]
    idx <- match(unlist(gs[,1]), trms)
    m[idx, i] <- as.numeric(gs[, scoreName])
  }
  m[is.na(m)] <- 0	
  rownames(m) <- trms
  colnames(m) <- names(gsea_list)
  return(m)
}

getHeatmapMatrix.pv <- function(gsea_list, adjusted = FALSE)
{
  trms <- sort(unique(unlist(lapply(gsea_list, function(i) i[,1]))))
  m <- matrix(NA, nrow = length(trms), ncol = length(gsea_list))
  
  for(i in 1:length(gsea_list))
  {
    gs <- gsea_list[[i]]
    idx <- match(unlist(gs[,1]), trms)
    if(adjusted)m[idx, i] <- as.numeric(gs$"pval")		
    else(m[idx, i] <- as.numeric(gs$"padj"))
  }
  m[is.na(m)] <- 1	
  rownames(m) <- trms
  colnames(m) <- names(gsea_list)
  return(m)
}

fgseaHeatmap <- function(fh_list, outFile, nb, adjusted = FALSE, kw = NULL, gs = NULL, group = NULL)
{
  gsMat <- getHeatmapMatrix.pv(fh_list, adjusted)
  colnames(gsMat) <- names(fh_list)	
  
  # re-order gsMat columns
  if(is.null(group)){
    idx.up <- grep("UP$", colnames(gsMat))
    idx.down <- grep("DOWN$", colnames(gsMat))
    gsMat <- gsMat[, c(idx.up, idx.down)]
    
    # add annotation
    ann.col <- data.frame(Sign = c(rep("UP", length(idx.up)), rep("DOWN", length(idx.down))))
    rownames(ann.col) <- colnames(gsMat)
    myColor.ann <- list(Sign = c(UP = "red", DOWN = "blue"))
  }else{
    upANDdown <- ifelse(grepl("UP$", colnames(gsMat)), "UP", "DOWN")
    upANDdown <- factor(upANDdown, levels = c("UP", "DOWN"))
    newOrder <- order(upANDdown, group)
    
    gsMat <- gsMat[, newOrder]
    upANDdown <- upANDdown[newOrder]
    group <- group[newOrder]
    
    # add annotation
    ann.col <- data.frame(Sign = upANDdown,
                          Group = group)
    rownames(ann.col) <- colnames(gsMat)
    myColor.ann <- list(Sign = c(UP = "red", DOWN = "blue"),
                        Group = getColor(group))
  }
  
  if(!is.null(kw)) gsMat <- gsMat[grepl(kw, rownames(gsMat), ignore.case = TRUE), ]
  if(!is.null(gs)) gsMat <- gsMat[match(intersect(rownames(gsMat), gs), rownames(gsMat)), ]
  
  # select top X gene-sets per column
  if(nrow(gsMat) > nb){
    idxMat <- apply(gsMat, 2, order)
    idxMat <- idxMat[1:nb, ]
    gsMat <- gsMat[unique(as.numeric(idxMat)), ]	
  }
  gsMat <- -log10(gsMat)	
  gsMat[gsMat > 5] <- 5 # set the limit to 5
  
  paletteLength <- 10
  myColor <- magma(paletteLength)
  myBreaks <- c(seq(0, max(gsMat), length.out=ceiling(paletteLength)))
  
  gsMat <- gsMat[, c("Dex.Pos-Core.UP", "Dex.Pos-Rim.UP", "Dex.Pos-Inv.UP", "Dex.Pos-Neg.UP",
                     "Din.Pos-Core.UP", "Din.Pos-Rim.UP", "Din.Pos-Inv.UP", "Din.Pos-Neg.UP",
                     "Dex.Pos-Core.DOWN", "Dex.Pos-Rim.DOWN", "Dex.Pos-Inv.DOWN", "Dex.Pos-Neg.DOWN",
                     "Din.Pos-Core.DOWN", "Din.Pos-Rim.DOWN", "Din.Pos-Inv.DOWN", "Din.Pos-Neg.DOWN")]
  
  
  #paletteLength <- 25
  #myMax <- ceiling(max(abs(gsMat)))
  
  #myBreaks <- seq(-myMax , myMax, length.out=paletteLength)
  #myBreaks <- myBreaks[myBreaks != 0]
  #myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength-2)
  
  doClust <- ifelse(nrow(gsMat) > 1, TRUE, FALSE)
  
  mysize <- nrow(gsMat) * 15 / 40
  mysize <- max(c(mysize, 10))
  if(mysize > 30) mysize <- 30
  
  pheatmap::pheatmap(gsMat, color = myColor, breaks = myBreaks, filename = outFile,
                     cluster_cols = FALSE, cluster_rows = doClust,
                     annotation_col = ann.col, annotation_colors = myColor.ann,
                     cellwidth = 10, cellheight = 10, fontsize_row = 8, fontsize_col = 8
                     #width = mysize, height = mysize
  )
}






############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################



# PARAMETERS
countDir <- file.path("../Data/")

outDir <- file.path("../eisaR/")
dir.create(outDir)

# INPUT FILES
exonCount <- "Final_count_exon.csv"
intronCount <- "Final_count_intron.csv"
sampleAnnotation <- "Annotation_col.csv"


############################
# COUNT
setwd(countDir)
cntEx <- read.csv(exonCount)
rownames(cntEx) <- cntEx[,1]
cntEx <- cntEx[, -1]

cntIn <- read.csv(intronCount)
rownames(cntIn) <- cntIn[,1]
cntIn <- cntIn[, -1]

# RUN EISA
setwd(outDir)

# remove "width" column
Rex <- cntEx[, colnames(cntEx) != "width"]
Rin <- cntIn[, colnames(cntIn) != "width"]

# SAMPLE ANNOTATION
ann.sample <- read.csv(sampleAnnotation)
ann.sample$SAMPLE <- gsub(" ", ".", ann.sample$id)
ann.sample <- ann.sample[match(colnames(Rex), ann.sample$SAMPLE), ]


###################
# EISA STEP BY STEP

# Normalization
Rex <- cntEx[,colnames(cntEx) != "width"]
Rin <- cntIn[,colnames(cntIn) != "width"]
Rall <- Rex + Rin
fracIn <- colSums(Rin)/colSums(Rall)
summary(fracIn)

# scale counts to the mean library size,
# separately for exons and introns
Nex <- t(t(Rex) / colSums(Rex) * mean(colSums(Rex)))
Nin <- t(t(Rin) / colSums(Rin) * mean(colSums(Rin)))

# log transform (add a pseudocount of 8)
NLex <- log2(Nex + 8)
NLin <- log2(Nin + 8)

# identify quantifiable genes
quantGenes <- rownames(Rex)[ rowMeans(NLex) > 5.0 & rowMeans(NLin) > 5.0 ]
length(quantGenes)


# RUN EDGER
myContrast <- matrix(c("Pos", "Core",
                       "Pos", "Rim",
                       "Pos", "Inv",
                       "Pos", "Neg"),
                     ncol = 2, byrow = TRUE)

lapply(1:nrow(myContrast), function(i){
  
  idx.test <- which(ann.sample$Annotation == myContrast[i, 1])
  idx.ref <- which(ann.sample$Annotation == myContrast[i, 2])
  
  # Calculate delta
  Dex <- NLex[,idx.test] - NLex[,idx.ref]
  Din <- NLin[,idx.test] - NLin[,idx.ref]
  Dex.Din <- Dex - Din
  
  # Save delta
  write.xlsx(list(Dex = Dex, Din = Din), file.path(outDir, paste0(myContrast[i, 1], "-", myContrast[i, 2], "_delta_matrices.xlsx")), row.names = TRUE)
  write.xlsx(list(Dex = cor(Dex), Din = cor(Din)), file.path(outDir, paste0(myContrast[i, 1], "-", myContrast[i, 2], "_delta_matrices.cor.xlsx")), row.names = TRUE)
  
  
  # DGE object
  cnt <- data.frame(Ex = Rex[, c(idx.test, idx.ref)], In = Rin[, c(idx.test, idx.ref)])
  y <- DGEList(counts = cnt, genes = data.frame(entrez = rownames(cnt),
                                                symbol = entrez2symbol(rownames(cnt))))
  
  # select quantifiable genes and normalize
  y <- y[quantGenes, ]
  y <- calcNormFactors(y)
  
  # design matrix with interaction term
  region <- factor(c(rep("ex", length(c(idx.test, idx.ref))), rep("in", length(c(idx.test, idx.ref)))), levels = c("in", "ex"))
  cond <- rep(factor(ann.sample$Annotation[c(idx.test, idx.ref)], levels = c(myContrast[i, 2], myContrast[i, 1])), 2)
  design <- model.matrix(~0+ region * cond)
  rownames(design) <- colnames(cnt)
  #colnames(design) <- gsub("cond", "", colnames(design))
  
  # estimate model parameters
  y <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  
  # calculate likelihood-ratio between full and reduced models
  lrt <- glmLRT(fit)
  
  # create results table
  tt <- topTags(lrt, n = nrow(y), sort.by = "none")
  toxlsx <- tt$table[order(tt$table$PValue, decreasing = FALSE), ]
  
  # ADD Dex / Din / Dex.Din
  toxlsx$Dex <- rowMeans(Dex[match(toxlsx$entrez, rownames(Dex)), ])
  toxlsx$Din <- rowMeans(Din[match(toxlsx$entrez, rownames(Din)), ])
  toxlsx$Dex.Din <- toxlsx$Dex - toxlsx$Din
  
  # SAVE
  outName <- file.path(outDir, paste0(myContrast[i, 1], "-", myContrast[i, 2], "_edgeR.xlsx"))
  write.xlsx(toxlsx, outName)
})

# DELTA SCATTERPLOT (GGPLOT)
setwd(outDir)
degFiles <- list.files(pattern = "_edgeR.xlsx$")
names(degFiles) <- gsub("_edgeR.xlsx", "", degFiles)
degList <- lapply(degFiles, read.xlsx, sheet = 1)

pvCutoff <- 0.05

addGenes <- c("IER2", "PLS1", "MOV10L1", "CCL2", "MMP25", "ADAMTSL5",
              "MMP19", "DMD", "PLCXD3", "MGLL", "BTBD11",
              "NDUFA7", "NDUDB7", "NDUFB1", "ALKBH7", "MRPS15",
              "CCL4", "IER2", "CCL2", "PANX1",
              "SIGLEC9", "CCL2")

lapply(1:length(degList), function(i){
  
  ggmat <- degList[[i]][, c("symbol", "PValue", "FDR", "Dex", "Din", "Dex.Din")]
  ggmat$GROUP <- "NS"
  ggmat$GROUP[ggmat$PValue < pvCutoff & ggmat$Dex.Din > 0 ] <- "Dex_UP"
  ggmat$GROUP[ggmat$PValue < pvCutoff & ggmat$Dex.Din < 0] <- "Dex_DOWN"
  
  ggmat$GROUP2 <- "NS"
  ggmat$GROUP2[ggmat$FDR < pvCutoff & ggmat$Dex > 0 & ggmat$Din > 0] <- "DexPOS_DinPOS"
  ggmat$GROUP2[ggmat$FDR < pvCutoff & ggmat$Dex > 0 & ggmat$Din < 0] <- "DexPOS_DinNEG"
  ggmat$GROUP2[ggmat$FDR < pvCutoff & ggmat$Dex < 0 & ggmat$Din > 0] <- "DexNEG_DinPOS"
  ggmat$GROUP2[ggmat$FDR < pvCutoff & ggmat$Dex < 0 & ggmat$Din < 0] <- "DexNEG_DinNEG"
  
  ggmat$toShow <- "no"
  ggmat$toShow[ggmat$symbol %in% addGenes & ggmat$FDR < pvCutoff] <- "yes"
  
  # PLOT
  p <- ggplot(ggmat, aes(x=Din, y=Dex, label = symbol))
  p <- p + geom_abline(intercept = 0, linetype= 2)
  p <- p + geom_point(aes(colour = GROUP), alpha = 0.75, shape=16, size = -log10(ggmat$PValue)[ggmat$GROUP == "NS"], data = subset(ggmat, GROUP == "NS"))
  p <- p + geom_point(aes(colour = GROUP), alpha = 0.5, shape=16, size = -log10(ggmat$PValue)[ggmat$GROUP != "NS"], data = subset(ggmat, GROUP != "NS"))
  p <- p + geom_point(size = -log10(ggmat$PValue)[ggmat$toShow == "yes"], shape = 21, color = "black", data = subset(ggmat, toShow == "yes"))
  p <- p + scale_color_manual(values = c(NS = "grey",
                                         Dex_UP = rgb(196, 49, 42, maxColorValue = 255),
                                         Dex_DOWN = rgb(70, 107, 161, maxColorValue = 255)))
  p <- p + theme_bw(base_size = 16) + theme(legend.position="none") 
  #p <- p + geom_text_repel(data=geneWiseMat[geneWiseMat$score > -log10(0.05),], aes(label=symbol), cex = 2, segment.size = 0.1)
  p <- p + geom_text_repel(data = subset(ggmat, toShow == "yes"), color = "black", box.padding = 0.5, max.overlaps = Inf)
  
  ggsave(plot = p, filename = paste0(names(degList)[i], "_delta_scatterplot.pdf"), width = 7, height = 7)
  
  
})



##############################################################
# CHECK GENES OF INTEREST

# LOAD EDGER
setwd(outDir)
degFiles <- list.files(pattern = "_edgeR.xlsx$")
names(degFiles) <- gsub("_edgeR.xlsx", "", degFiles)
degList <- lapply(degFiles, read.xlsx, sheet = 1)
deg <- unique(unlist(lapply(degList, function(i) i$symbol[i$FDR < 0.05])))


# LOAD GENES OF INTEREST
geneMat <- read.xlsx(file.path(dataDir, "Six cat_Leading edge genes.xlsx"))
geneMat <- geneMat[geneMat$Genes %in% deg, ]

geneList <- lapply(unique(geneMat$Cat), function(i) geneMat$Genes[geneMat$Cat == i])
names(geneList) <- unique(geneMat$Cat)

lapply(1:length(geneList), function(i){
  degList.sub <- lapply(degList, function(j){
    j[match(geneList[[i]], j$symbol), ]
  })
  
  # DELTA
  deltaMat <- lapply(degList.sub, function(j) j$Dex.Din)
  deltaMat <- do.call(cbind, deltaMat)
  rownames(deltaMat) <- geneList[[i]]
  
  # CORPLOT LIKE
  pvMat <- lapply(degList.sub, function(j) j$FDR)
  pvMat <- do.call(cbind, pvMat)
  rownames(pvMat) <- geneList[[i]]
  
  # RE-ORDER COLUMNS
  deltaMat <- deltaMat[, c("Pos-Core", "Pos-Rim", "Pos-Inv", "Pos-Neg"), drop = FALSE]
  pvMat <- pvMat[, c("Pos-Core", "Pos-Rim", "Pos-Inv", "Pos-Neg"), drop = FALSE]
  
  # RE-ORDER ROWS
  if(nrow(deltaMat)>1){
    newOrder <- order(-rowMeans(deltaMat))
    deltaMat <- deltaMat[newOrder, ]
    pvMat <- pvMat[newOrder, ]
  }
  
  mymax <- ceiling(max(abs(deltaMat)))
  pdf(file.path(outDir,
                paste0(names(geneList)[i], "_Dex.Din.delta_heatmap_V2.pdf")))
  corrplot(deltaMat, is.corr = FALSE,
           p.mat = pvMat, sig.level = 0.05, insig = "label_sig",
           cl.lim = c(-mymax, mymax),
           tl.col = "black", tl.cex = 1.5,
           cl.cex = 1.5, cl.ratio = 0.3,
           col = rev(brewer.pal(n = 8, name = "RdYlBu")))
  dev.off()
  
  # SAVE
  degList.sub <- lapply(degList.sub, function(j) j[order(j$PValue), ])
  write.xlsx(degList.sub,
             file.path(outDir,
                       paste0(names(geneList)[i], "_edgerR.xlsx")))
  
})



##############################################################
# SIGNIFICANCE OF GENES OF INTEREST

# LOAD DELTA
deltaFiles <- list.files(pattern = "_delta_matrices.xlsx")
names(deltaFiles) <- gsub("_delta_matrices.xlsx", "", deltaFiles)
DexList <- lapply(deltaFiles, read.xlsx, sheet = "Dex", rowNames = TRUE)
DinList <- lapply(deltaFiles, read.xlsx, sheet = "Din", rowNames = TRUE)
DexDinList <- lapply(1:length(DexList), function(i) DexList[[i]] - DinList[[i]])
names(DexDinList) <- names(DexList)

DexMat <- lapply(DexList, rowMeans)
DexMat <- do.call(cbind, DexMat)

DinMat <- lapply(DinList, rowMeans)
DinMat <- do.call(cbind, DinMat)

# LOAD GENES OF INTEREST
geneMat <- read.xlsx(file.path("../Data/ALA_signature.xlsx"))
geneList <- lapply(1:ncol(geneMat), function(i) as.character(geneMat[, i]))
geneList <- lapply(geneList, function(i) i[!is.na(i)])
geneList <- lapply(geneList, unique)
names(geneList) <- colnames(geneMat)
geneList <- lapply(geneList, intersect, degList[[1]]$symbol)
geneList.entrez <- lapply(geneList, function(i) degList[[1]]$entrez[match(i, degList[[1]]$symbol)])


# AVG DELTA
geneList.dex <- lapply(geneList.entrez, function(i) colMeans(DexMat[i, ]))
geneList.din <- lapply(geneList.entrez, function(i) colMeans(DinMat[i, ]))

nb <- 10000
geneList.dex.rd <- lapply(geneList.entrez, function(i){
  m <- lapply(1:nb, function(j){
    genes.rd <- sample(rownames(DexMat), length(i))
    colMeans(DexMat[genes.rd, ])
  })
  return(do.call(rbind, m))
})

geneList.dex.pv <- lapply(1:length(geneList.dex), function(i){
  
  pv.greater <- sapply(1:length(geneList.dex[[i]]), function(j){
    getEmpiricalPV(x = geneList.dex.rd[[i]][, j], y = geneList.dex[[i]][j], type = "greater")
  })
  
  pv.less <- sapply(1:length(geneList.dex[[i]]), function(j){
    getEmpiricalPV(x = geneList.dex.rd[[i]][, j], y = geneList.dex[[i]][j], type = "less")
  })
  
  return(data.frame(COMP = names(geneList.dex[[i]]),
                    PV.greater = pv.greater,
                    PV.less = pv.less))
  
})


# FGSEA
DexMat.rank <- lapply(1:ncol(DexMat), function(i){
  x <- as.numeric(DexMat[, i])
  names(x) <- rownames(DexMat)
  return(sort(x))
})
names(DexMat.rank) <- paste0("Dex.", colnames(DexMat))

DinMat.rank <- lapply(1:ncol(DinMat), function(i){
  x <- as.numeric(DinMat[, i])
  names(x) <- rownames(DinMat)
  return(sort(x))
})
names(DinMat.rank) <- paste0("Din.", colnames(DinMat))

fgseaDir <- file.path(outDir, "fgsea")
dir.create(fgseaDir)
setwd(fgseaDir)

dbList <- list(Signature = geneList.entrez)

rankedGenesList <- c(DexMat.rank, DinMat.rank)

lapply(1:length(dbList), function(i){
  mydb <- dbList[[i]]
  mydb.name <- names(dbList)[i]
  
  plotDir <- fgseaDir   
  
  dir.create(file.path(plotDir, mydb.name), showWarnings = FALSE)
  setwd(file.path(plotDir, mydb.name))   
  
  #########################
  # Perform fgsea analysis
  
  fgseaResList <- lapply(rankedGenesList,
                         fgsea, pathways = mydb, minSize=5, maxSize=1000)
  fgseaResList <- lapply(fgseaResList, function(i) i[order(i$pval),])	
  
  # Add symbols
  fgseaResList <- lapply(fgseaResList, as.data.frame)
  
  fgseaResList <- lapply(fgseaResList, function(k){
    k.symbol <- lapply(k$leadingEdge, entrez2symbol)
    k.symbol <- lapply(k.symbol, function(j) unique(j[!is.na(j)]))
    k.symbol <- lapply(k.symbol, toString)
    k$leadingEdge.symbol <- unlist(k.symbol)
    k
  })
  
  fgseaResList <- lapply(fgseaResList, function(j){
    j$NES[is.na(j$NES)] <- 0
    return(j)
  })
  
  # Divide UP and DOWN
  fgseaResList.UP <- lapply(fgseaResList, function(j) j[j$NES > 0, ])
  names(fgseaResList.UP) <- paste0(names(fgseaResList), ".UP")
  fgseaResList.DOWN <- lapply(fgseaResList, function(j) j[j$NES < 0, ])
  names(fgseaResList.DOWN) <- paste0(names(fgseaResList), ".DOWN")
  fgseaResList <- c(fgseaResList.UP, fgseaResList.DOWN)
  
  # Save fgsea output   
  lapply(1:length(rankedGenesList), function(j)
    write.xlsx(list(UP =  fgseaResList.UP[[j]], DOWN = fgseaResList.DOWN[[j]]),
               paste(names(rankedGenesList)[j], "_", mydb.name, "_fgsea.xlsx", sep = ""),
               row.names = FALSE, firstRow = T, headerStyle = createStyle(textDecoration = 'bold'))
  )	
  
})

##########
# HEATMAPS

#comp <- c("A-B",
#          "B-C")
comp <- names(rankedGenesList)

setwd(fgseaDir)
lapply(1:length(dbList), function(i){
  mydb.name <- names(dbList)[i]
  
  # load fisher results
  fhFiles <- file.path(mydb.name, paste0(comp, "_", mydb.name, "_fgsea.xlsx"))
  
  fhList <- lapply(fhFiles, my.read_xlsx)
  names(fhList) <- gsub(paste0("_", mydb.name, "_fgsea.xlsx"), "", fhFiles)
  names(fhList) <- gsub(paste0(mydb.name, "\\/"), "", names(fhList))
  fhList <- do.call(c, fhList)
  
  fgseaHeatmap(fhList, file.path(mydb.name, paste0(mydb.name, "_fgsea_heatmap.pdf")), nb = 5, adjusted = FALSE,
               group = NULL)
})







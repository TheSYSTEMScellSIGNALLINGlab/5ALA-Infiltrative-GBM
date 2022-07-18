
library(org.Hs.eg.db)
library(openxlsx)
library(pheatmap)
library(GO.db)
library(fgsea)
library(GeneAnswers)
library(viridis)
library(gridExtra)
library(msigdbr)
library(ggplot2)
library(ungeviz)
library(limma)
library(tidyr)


############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

getKbLength <- function(entrez, organism = "hsa"){
  require(EDASeq)
  lgth <- getGeneLengthAndGCContent(entrez, organism, mode = "org.db")
  return(lgth[,1] / 1000)
}

getTPM <- function(countMat, organism = "hsa"){
  
  geneLength <- getKbLength(rownames(countMat), organism)
  
  # remove NA gene length
  toRemove <- is.na(geneLength)
  geneLength <- geneLength[!toRemove]
  countMat <-countMat[!toRemove, ]
  
  # get TPM
  x <- countMat / geneLength
  tpmMat <- t( t(x) * 1e6 / colSums(x) )
  
  return(tpmMat)
}


my.read_xlsx <- function(inFile)
{
  mysheets <- getSheetNames(inFile)
  mList <- lapply(mysheets, read.xlsx, xlsxFile = inFile)
  names(mList) <- mysheets
  return(mList)
}


fgsea_pipeline <- function(rank, pathways, minSize=15, maxSize=500, nperm=10000){
  fg <- fgsea(rank, pathways = pathways, minSize=minSize, maxSize=maxSize, nperm=nperm)
  fg <- as.data.frame(fg)
  fg.entrez <- fg$leadingEdge
  fg.symbol <- lapply(fg.entrez, entrez2symbol)
  fg.symbol <- unlist(lapply(fg.symbol, function(i) paste(i, collapse = ",")))
  avg.delta <- unlist(lapply(fg.entrez, function(i) mean(rank[i])))
  fg <- data.frame(fg, leadingEdge.symbol = fg.symbol, avgRank = avg.delta)
  fg <- fg[order(fg$pval), ]
  return(fg)
}


symbol2entrez <- function(symbol)
{
  entrez <- mget(as.character(symbol), org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrez <- unique(unlist(lapply(entrez, function(i) return(i[1]))))
  entrez <- entrez[!is.na(entrez)]
  return(entrez)
}

entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Hs.egSYMBOL, ifnotfound=NA)
  symbol <- unique(unlist(lapply(symbol, function(i) return(i[1]))))
  return(symbol)
}

ensembl2entrez <- function(ensembl)
{
  entrez <- mget(as.character(ensembl), org.Hs.egENSEMBL2EG, ifnotfound=NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

refseq2entrez <- function(refseq){
  entrez <- mget(as.character(refseq), org.Hs.egREFSEQ2EG, ifnotfound=NA)
  entrez <- unique(unlist(lapply(entrez, function(i) return(i[1]))))
  return(entrez)
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

fgseaHeatmapSingle <- function(fh_list, outFile, nb, adjusted = FALSE, kw = NULL, gs = NULL)
{
  gsMat <- getHeatmapMatrix.pv(fh_list, adjusted)
  colnames(gsMat) <- names(fh_list)	
  
  if(!is.null(kw)) gsMat <- gsMat[grepl(kw, rownames(gsMat), ignore.case = TRUE), ]
  if(!is.null(gs)) gsMat <- gsMat[match(intersect(rownames(gsMat), gs), rownames(gsMat)), ]
  
  # select top X gene-sets per column
  idxMat <- apply(gsMat, 2, order)
  if(nrow(idxMat) > nb) idxMat <- idxMat[1:nb, ]
  gsMat <- gsMat[unique(as.numeric(idxMat)), ]	
  
  # Merge UP and DOWN columns
  idx.up <- grep("UP$", colnames(gsMat))
  idx.down <- grep("DOWN$", colnames(gsMat))
  gsMat.up <- gsMat[, idx.up]
  gsMat.down <- gsMat[, idx.down]
  
  gsMat <- matrix(NA, nrow = nrow(gsMat.up), ncol = ncol(gsMat.up))
  for(rowIdx in 1:nrow(gsMat)){
    for(colIdx in 1:ncol(gsMat)){
      pv.up <- gsMat.up[rowIdx, colIdx] 
      pv.down <- gsMat.down[rowIdx, colIdx]
      if(pv.up <= pv.down) gsMat[rowIdx, colIdx] <- -log10(pv.up)
      if(pv.up > pv.down) gsMat[rowIdx, colIdx] <- log10(pv.down)
    }
  }
  rownames(gsMat) <- rownames(gsMat.up)
  colnames(gsMat) <- gsub("\\.UP", "", colnames(gsMat.up))
  
  gsMat[gsMat > 6] <- 6 # set the limit to 5
  gsMat[gsMat < -6] <- -6
  
  paletteLength <- 25
  myMax <- ceiling(max(abs(gsMat)))
  
  myBreaks <- seq(-myMax , myMax, length.out=paletteLength)
  myBreaks <- myBreaks[myBreaks != 0]
  myColor <- colorRampPalette(c("deepskyblue", "snow", "orangered"))(paletteLength-2)
  
  doClust <- ifelse(nrow(gsMat) > 1, TRUE, FALSE)
  
  mysize <- nrow(gsMat) * 15 / 40
  mysize <- max(c(mysize, 10))
  if(mysize > 30) mysize <- 30
  
  pheatmap::pheatmap(gsMat, color = myColor, breaks = myBreaks, filename = outFile,
                     annotation_row = NULL, annotation_col = NULL,
                     cluster_cols = TRUE, cluster_rows = doClust, show_rownames = TRUE,
                     cellwidth = 10, cellheight = 10, fontsize_row = 8, fontsize_col = 8,
                     width = mysize, height = mysize)
}


scaledAnnoHeatmap <- function(mat, ann.row, ann.col, outFile)
{
  paletteLength <- 25
  myMax <- ceiling(max(abs(mat)))
  
  myBreaks <- seq(-myMax , myMax, length.out=paletteLength)
  myBreaks <- myBreaks[myBreaks != 0]
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength-2)
  
  #myColor.ann <- list(Group = c(Dc = "black", Tc = "grey", WEHI = "red"))
  
  if(nrow(mat)<2) doCluster <- FALSE
  else(doCluster <- TRUE)
  
  pheatmap(mat, color = myColor, breaks = myBreaks, filename = outFile,
           annotation_row = ann.row, annotation_col = ann.col,
           #annotation_colors = myColor.ann,
           cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE,
           cellwidth = 12, cellheight = 12
           #height = 8, width = 8
  )
  
}	

############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


# PARAMETERS
dataDir <- file.path("../Data")

fgseaDir <- file.path("../gsea/fgsea_single")
dir.create(fgseaDir, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES
exonCount <- "Final_count_exon.csv"
sampleAnnotation <- "Annotation_col.csv"
signature5ALA <- "ALA_signature.xlsx"

# BUILD GENE-SETS FROM FILE
mygs.raw <- read.xlsx(file.path(dataDir, signature5ALA), sheet = 1)
mygs <- lapply(1:ncol(mygs.raw), function(i) as.character(mygs.raw[, i]))
mygs <- lapply(mygs, function(i) i[!is.na(i)])
mygs <- lapply(mygs, symbol2entrez)
names(mygs) <- colnames(mygs.raw)

dbList <- list("leadingEdge" = mygs)


##################################
# Rank genes based on TPM

setwd(dataDir)
countMat <- read.csv(exonCount)
rownames(countMat) <- countMat[,1]
countMat <- countMat[, -c(1:2)]

# COUNT TO TPM
tpmMat <- log2(getTPM(countMat+3))

rankedGenesList <- lapply(1:ncol(tpmMat), function(i){
  tpm.current <- as.numeric(tpmMat[,i])
  names(tpm.current) <- rownames(tpmMat)
  return(sort(tpm.current))
})
names(rankedGenesList) <- colnames(tpmMat)

# SAMPLE ANNOTATION
setwd(dataDir)
ann.sample <- read.csv(sampleAnnotation)
ann.sample$id <- gsub(" ", ".", ann.sample$id)
ann.sample <- ann.sample[match(names(rankedGenesList), ann.sample$id), ]
ann.sample$patient <- sapply(strsplit(ann.sample$id, split = "\\."), function(i) i[2])
ann.sample$patient <- sapply(strsplit(ann.sample$patient, split = "_"), function(i) i[1])


#########################
# Perform fgsea analysis

# FGSEA
setwd(fgseaDir)
lapply(1:length(dbList), function(i){
  mydb <- dbList[[i]]
  mydb.name <- names(dbList)[i]
  
  plotDir <- fgseaDir   
  
  dir.create(file.path(plotDir, mydb.name), showWarnings = FALSE)
  setwd(file.path(plotDir, mydb.name))   
  
  #########################
  # Perform fgsea analysis
  
  fgseaResList <- lapply(rankedGenesList,
                         fgsea, pathways = mydb, minSize=5, maxSize=1000
  )
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












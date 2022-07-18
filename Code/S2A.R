
library(gage)
library(GO.db)
library(org.Hs.eg.db)
library(openxlsx)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)
library(ggplot2)
library(gridExtra)
library(msigdbr)
library(limma)
library(edgeR)

# BUILD DB
m_df <- msigdbr(species = "Homo sapiens")
m_df$gs_subcat <- gsub("\\:", "_", m_df$gs_subcat)
m_df$gs_subcat <- paste0(m_df$gs_cat, "_", m_df$gs_subcat)
m_df$gs_subcat <- gsub("_$", "", m_df$gs_subcat)

subcat.unique <- sort(unique(m_df$gs_subcat))

dbList <- lapply(subcat.unique, function(i){
  m_df.sub <- m_df[m_df$gs_subcat == i, ]
  terms <- unique(m_df.sub$gs_name)
  db <- mclapply(terms, function(j){
    return(m_df.sub$entrez_gene[m_df.sub$gs_name == j])
  }, mc.cores = 8)
  names(db) <- terms
  return(db)
})
names(dbList) <- subcat.unique


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


gagePipeline <- function(test, ref, db, up = NULL, down = NULL)
{
  gageMat <- cbind(ref, test)
  gos <- gage(gageMat, gsets = db, ref = seq(1:ncol(ref)),
              same.dir = T, compare = "paired", rank = T,
              saaTest = gs.tTest, set.size = c(5, 500))
  significant.groups <- sigGeneSet(gos, cutoff = 1, qpval = c("q.val"))
  
  # Greater
  greater <- significant.groups$greater
  greater.entrez <- lapply(rownames(greater), function(i) return(intersect(db[[i]], up)))
  greater.symbol <- lapply(greater.entrez, entrez2symbol)
  greater.nb <- lapply(greater.entrez, length)
  
  greater.entrez <- unlist(lapply(greater.entrez, toString))
  greater.symbol <- unlist(lapply(greater.symbol, toString))
  greater.nb <- unlist(greater.nb)
  greater <- cbind(greater, nb = greater.nb, entrez = greater.entrez, symbol = greater.symbol)
  
  # Less
  less <- significant.groups$less
  less.entrez <- lapply(rownames(less), function(i) return(intersect(db[[i]], down)))
  less.symbol <- lapply(less.entrez, entrez2symbol)
  less.nb <- lapply(less.entrez, length)
  
  less.entrez <- unlist(lapply(less.entrez, toString))
  less.symbol <- unlist(lapply(less.symbol, toString))
  less.nb <- unlist(less.nb)
  less <- cbind(less, nb = less.nb, entrez = less.entrez, symbol = less.symbol)
  
  return(list(greater = greater, less = less))
}


my.read_xlsx <- function(inFile)
{
  mysheets <- getSheetNames(inFile)
  mList <- lapply(mysheets, read.xlsx, xlsxFile = inFile)
  names(mList) <- mysheets
  return(mList)
}

my.read_xlsxNoHeader <- function(inFile)
{
  mysheets <- getSheetNames(inFile)
  mList <- lapply(mysheets, read.xlsx, xlsxFile = inFile, colNames = FALSE)
  names(mList) <- mysheets
  return(mList)
}

getHeatmapMatrix <- function(gsea_list, adjusted = FALSE)
{
  trms <- sort(unique(unlist(lapply(gsea_list, function(i) i[,1]))))
  m <- matrix(NA, nrow = length(trms), ncol = length(gsea_list))
  
  for(i in 1:length(gsea_list))
  {
    gs <- gsea_list[[i]]
    idx <- match(gs[,1], trms)
    if(adjusted)m[idx, i] <- as.numeric(gs$"q.val")		
    else(m[idx, i] <- as.numeric(gs$"p.val"))
  }
  m[is.na(m)] <- 1	
  rownames(m) <- trms
  colnames(m) <- names(gsea_list)
  return(m)
}

gageHeatmap <- function(fh_list, outFile, nb, adjusted = FALSE, kw = NULL, gs = NULL)
{
  gsMat <- getHeatmapMatrix(fh_list, adjusted)	
  colnames(gsMat) <- names(fh_list)	
  
  # re-order gsMat columns
  idx.up <- grep("UP$", colnames(gsMat))
  idx.down <- grep("DOWN$", colnames(gsMat))
  gsMat <- gsMat[, c(idx.up, idx.down)]
  
  # add annotation
  ann.col <- data.frame(Sign = c(rep("UP", length(idx.up)), rep("DOWN", length(idx.down))))
  rownames(ann.col) <- colnames(gsMat)
  myColor.ann <- list(Sign = c(UP = "red", DOWN = "blue"))
  
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
  #myColor <- viridis(paletteLength)
  #myColor <- colorRampPalette(c("white", "red"))(paletteLength)
  myBreaks <- c(seq(0, max(gsMat), length.out=ceiling(paletteLength)))
  
  doClust <- ifelse(nrow(gsMat) > 1, TRUE, FALSE)
  
  mysize <- nrow(gsMat) * 15 / 40
  mysize <- max(c(mysize, 10))
  if(mysize > 30) mysize <- 30
  
  pheatmap::pheatmap(gsMat, color = myColor, breaks = myBreaks, filename = outFile,
                     cluster_cols = FALSE, cluster_rows = doClust,
                     annotation_col = ann.col, annotation_colors = myColor.ann,
                     cellwidth = 10, cellheight = 10, fontsize_row = 8, fontsize_col = 8,
                     width = mysize, height = mysize)
}

gageHeatmapSingle <- function(fh_list, outFile, nb, adjusted = FALSE, kw = NULL, gs = NULL)
{
  gsMat <- getHeatmapMatrix(fh_list, adjusted)	
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


############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

# BUILD GENE-SETS FROM FILE
mygs.raw <- read.xlsx(file.path("~/Research/Sajib/P10_GBM/doc/Gene_Set_ALA_POS.xlsx"), sheet = 1)
mygs <- lapply(1:ncol(mygs.raw), function(i) as.character(mygs.raw[, i]))
mygs <- lapply(mygs, function(i) i[!is.na(i)])
mygs <- lapply(mygs, symbol2entrez)
names(mygs) <- colnames(mygs.raw)

dbList[["leadingEdge"]] <- mygs


# PARAMETERS
countDir <- file.path("../Data/")

limmaDir <- file.path("../limma_voom/")

gageDir <- file.path("../gsea/gage_paired/")
dir.create(gageDir, recursive = TRUE, showWarnings = FALSE)


# INPUT FILES
exonCount <- "Final_count_exon.csv"
sampleAnnotation <- "Annotation_col.csv"



##############################################################
# EXPRESSION
setwd(countDir)

setwd(countDir)
cntEx <- read.csv(exonCount)
rownames(cntEx) <- cntEx[,1]
cntEx <- cntEx[, -1]
Rex <- cntEx[, colnames(cntEx) != "width"]
countMat <- Rex

# Create geneMat
entrezID <- rownames(countMat)
symbol <- entrez2symbol(entrezID)

geneMat <- data.frame(entrez = entrezID, symbol = symbol)

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
expr <- cpm(m, log = TRUE)

# DEFINE GROUPS
Group <- ann.sample$Annotation
cond <- Group


# CONTRAST
contMat <-  matrix(c("Pos", "Core",
                     "Pos", "Rim",
                     "Pos", "Inv",
                     "Pos", "Neg"),
                   ncol = 2, byrow = TRUE)
contMat <- cbind("Group", contMat)					 
contName <- apply(contMat, 1, function(i) paste(i[2], i[3], sep = "-"))


######################
lapply(1:length(dbList), function(i){
  
  mydb <- dbList[[i]]
  mydb.name <- names(dbList)[i]
  
  lapply(1:nrow(contMat), function(j){
    testCond <- contMat[j, 2]
    refCond <- contMat[j, 3]
    
    deseqMat <- read.xlsx(file.path(limmaDir, paste0(contName[j], "_limma.xlsx")), sheet=1)
    deseqMat <- deseqMat[!is.na(deseqMat$entrez), ]
    deseqFC <- deseqMat$logFC
    deseqUP <- deseqMat$entrez[deseqFC > 0]
    deseqDOWN <- deseqMat$entrez[deseqFC < 0]
    
    testMat <- expr[,cond==testCond, drop = FALSE]
    refMat <- expr[,cond==refCond, drop = FALSE]
    gage.res <- gagePipeline(testMat, refMat, mydb, deseqUP, deseqDOWN)
    
    dir.create(file.path(gageDir, mydb.name), showWarnings = FALSE)
    
    write.xlsx(list(UP = gage.res$greater, DOWN = gage.res$less),
               file.path(gageDir, mydb.name, paste(contName[j], "_", mydb.name, "_gage.xlsx", sep = "")),
               row.names = TRUE, firstRow = T, headerStyle = createStyle(textDecoration = 'bold'), overwrite = TRUE)
    
  })
  
})

##########
# HEATMAPS

comp <- contName

setwd(gageDir)
lapply(1:length(dbList), function(i){
  mydb.name <- names(dbList)[i]
  
  # load fisher results
  fhFiles <- file.path(mydb.name, paste0(comp, "_", mydb.name, "_gage.xlsx"))
  
  fhList <- lapply(fhFiles, my.read_xlsx)
  names(fhList) <- gsub(paste0("_", mydb.name, "_gage.xlsx"), "", fhFiles)
  names(fhList) <- gsub(paste0(mydb.name, "\\/"), "", names(fhList))
  fhList <- do.call(c, fhList)
  
  gageHeatmapSingle(fhList, file.path(mydb.name, paste0(mydb.name, "_gage_heatmap.pdf")), nb = 5, adjusted = FALSE)
  
})





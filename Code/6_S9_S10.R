library(SPATA2)
library(openxlsx)
library(ggplot2)
library(devtools)
library(monocle3)
library(tidyverse)
library(RColorBrewer)
library(viridis)

# devtools::install_github(repo = "theMILOlab/SPATA2")
# devtools::install_github("theMILOlab/SPATA2", "runGSEA_debug", force = T)


############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################



############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


# DEFINE ALA5pos SIGNATURES
geneMat <- read.xlsx("spata_gene_list.xlsx")
geneList <- lapply(1:ncol(geneMat), function(i) as.character(geneMat[, i]))
names(geneList) <- colnames(geneMat)
geneList <- lapply(geneList, function(i) unique(i[!is.na(i)]))
# geneList$ALL <- unique(unlist(geneList))


# LOAD SPATA2 OBJ
dataDir <- file.path("H:/lab data/Glioblastoma/SPATA/CNV/")
setwd(dataDir)
spata_files <- list.files(pattern = ".RDS")

idx <- 1#243_T
idx <- 2#251_T
idx <- 3#259_T
idx <- 4#266_T
idx <- 5#269_T
idx <- 6#275_T
idx <- 7#313_T
idx <- 8#334_T


spata_obj <- loadSpataObject(spata_files[idx])


############### Autoencoder Denoising ############### 

# require(devtools)
# install_github("rstudio/reticulate")
# install_github("rstudio/tensorflow", force = T)
# install_github("rstudio/keras")
# tensorflow::install_tensorflow()
# tensorflow::tf_config()

# all expression matrices before denoising
getExpressionMatrixNames(object = spata_obj)

# active expression matrix before denoising
getActiveMatrixName(object = spata_obj)

# denoising your data 
spata_obj <-
  runAutoencoderDenoising(
    object = spata_obj, 
    activation = "selu", 
    bottleneck = 56, 
    epochs = 20, 
    layers = c(128, 64, 32), 
    dropout = 0.1
  )

# all expression matrices after denoising
getExpressionMatrixNames(object = spata_obj)

# active expression matrix after denoising
getActiveMatrixName(object = spata_obj)

# print a summary 
printAutoencoderSummary(object = spata_obj)

spata_obj <-
  runAutoencoderAssessment(
    object = spata_obj,
    activations = c("relu", "selu", "sigmoid"), 
    bottlenecks = c(32, 40, 48, 56, 64),
    epochs = 20, 
    layers = c(128, 64, 32), 
    dropout = 0.1
  )

plotAutoencoderAssessment(object = spata_obj)

spata_obj <- setActiveExpressionMatrix(object = spata_obj, mtr_name = "scaled")

spata_obj <- setActiveExpressionMatrix(object = spata_obj, mtr_name = "denoised")



# saveRDS(spata_obj, "334_T_AD.RDS")

############### Denoising end ##############



### Coordinate extraction
coordinate <- getCoordsDf(spata_obj)
write.csv(coordinate, "Coordinate_ti.csv")


### Expression matrix extraction
mat <- spata_obj@data
count <- mat$`334_T`$counts
count <- as.data.frame(count)
write.csv(count, "334_Ti_matrix_de.csv")

# initial features in slot @fdata
getFeatureNames(spata_obj)


outDir <- file.path("H:/lab data/Glioblastoma/SPATA/CNV/")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
setwd(outDir)

############################
# Basic extracting functions

if(FALSE){
  # the essential data.frame
  getSpataDf(spata_obj)
  
  # dimensional reduction data
  getUmapDf(spata_obj)
  
  # barcode spot coordinates
  getCoordsDf(spata_obj)
}


# GENES
spata.genes.meta <- getGeneMetaData(spata_obj)
spata.genes <- getGenes(spata_obj)
geneList <- lapply(geneList, intersect, y = spata.genes)


### Adding genesets

spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'ALA', genes = geneList$ALA)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'GPM', genes = geneList$GPM)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Hypoxia', genes = geneList$Hypoxia)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Endo', genes = geneList$Endo_marker)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'MTC', genes = geneList$MTC)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'IWR', genes = geneList$IWR)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'MES', genes = geneList$MES)


##############
# PLOT SURFACE
# open application to obtain a list of plots
plots <- plotSurfaceInteractive(object = spata_obj)

names(geneList)

Pltt <- plots$ALA +
        ggplot2::labs(title = "ALA") +
        theme(plot.title = element_text(size = 8, face = "bold"))+
        plots$GPM +
        ggplot2::labs(title = "GPM") +
        theme(plot.title = element_text(size = 8, face = "bold"))+
        plots$Hypoxia +
        ggplot2::labs(title = "Hypoxia") +
        theme(plot.title = element_text(size = 8, face = "bold"))
 
ggsave(plot = Pltt, file=file.path("H:\\lab data\\Glioblastoma\\SPATA\\CNV\\New analysis_5 sig/",
                                    "334_C_genesets.pdf"),
       width = 12, height = 10)



## Plotting chr
chr <-  plots$ch7 +
        ggplot2::labs(title = "ch7") +
        theme(plot.title = element_text(size = 8, face = "bold"))+
        plots$ch10 +
        ggplot2::labs(title = "ch10") +
        theme(plot.title = element_text(size = 8, face = "bold"))

ggsave(plot = chr, file=file.path("H:\\lab data\\Glioblastoma\\SPATA\\CNV\\New analysis_5 sig/",
                                    "334_T_Chr.pdf"),
       width = 12, height = 10)


# spata_obj <- setActiveExpressionMatrix(object = spata_obj, mtr_name = "scaled")
spata_obj <- setActiveExpressionMatrix(object = spata_obj, mtr_name = "denoised")


### surface plotting
p <- plotSurfaceAverage(object = spata_obj, 
                   color_by = geneList, 
                   smooth = TRUE, 
                   pt_size = 2,
                   smooth_span = 0.2)

p

ggsave(plot = p, file=file.path("H:\\lab data\\Glioblastoma\\SPATA\\CNV\\New analysis/",
                                "259_Ti_surface_scaled.pdf"),
       width = 12, height = 10)


# ### ssgsea
# coords_df <- getCoordsDf(spata_obj)
# 
# joined_df <-
#   joinWith(object = spata_obj,
#            spata_df = coords_df,
#            method_gs = "ssgsea",
#            gene_sets = c("mygs_ALA_pos_signature_251", "mygs_GPM_whole", "mygs_Injury_whole", "mygs_Hypoxia",
#                           "mygs_Endothelial_marker", "mygs_Mesenchymal_whole", "mygs_Inf_response_whole",
#                           "mygs_TNFA_whole", "mygs_GPM_part", "mygs_Injury_part", "mygs_Mesenchymal_part",
#                           "mygs_Inf_response_part", "mygs_TNFA_part"), # expression values of the gene sets
#            features = "seurat_clusters", # cluster belonging
#            verbose = FALSE)
# 
# write.csv(joined_df, "334_T_ssgsea.csv")

# # add to spata-object via 
# spata_obj <- addFeatures(object = spata_obj,
#                          feature_names = c("mygs_ALA_pos_signature_251", "mygs_GPM_whole", "mygs_Injury_whole", "mygs_Hypoxia",
#                                            "mygs_Endothelial_marker", "mygs_Mesenchymal_whole", "mygs_Inf_response_whole",
#                                            "mygs_TNFA_whole", "mygs_GPM_part", "mygs_Injury_part", "mygs_Mesenchymal_part",
#                                            "mygs_Inf_response_part", "mygs_TNFA_part"), 
#                          feature_df = joined_df, 
#                          overwrite = TRUE)

# # visualize ssgsea on the surface
# p1 <- plotSurface(spata_obj, color_by = "mygs_ALA_pos_signature_251", smooth = TRUE, smooth_span = 0.2, pt_size = 2.0) 
# p2 <- plotSurface(spata_obj, color_by = "mygs_GPM_whole", smooth = TRUE, smooth_span = 0.2, pt_size = 2.0) 
# p3 <- plotSurface(spata_obj, color_by = "mygs_Injury_whole", smooth = TRUE, smooth_span = 0.2, pt_size = 2.0) 
# p4 <- plotSurface(spata_obj, color_by = "mygs_Hypoxia", smooth = TRUE, smooth_span = 0.2, pt_size = 2.0) 
# p5 <- plotSurface(spata_obj, color_by = "mygs_Endothelial_marker", smooth = TRUE, smooth_span = 0.2, pt_size = 2.0) 
# p6 <- plotSurface(spata_obj, color_by = "mygs_Mesenchymal_whole", smooth = TRUE, smooth_span = 0.2, pt_size = 2.0) 
# p7 <- plotSurface(spata_obj, color_by = "mygs_Inf_response_whole", smooth = TRUE, smooth_span = 0.2, pt_size = 2.0) 
# p8 <- plotSurface(spata_obj, color_by = "mygs_TNFA_whole", smooth = TRUE, smooth_span = 0.2, pt_size = 2.0) 
# p9 <- plotSurface(spata_obj, color_by = "mygs_GPM_part", smooth = TRUE, smooth_span = 0.2, pt_size = 2.0) 
# p10 <- plotSurface(spata_obj, color_by = "mygs_Injury_part", smooth = TRUE, smooth_span = 0.2, pt_size = 2.0) 
# p11 <- plotSurface(spata_obj, color_by = "mygs_Mesenchymal_part", smooth = TRUE, smooth_span = 0.2, pt_size = 2.0) 
# p12 <- plotSurface(spata_obj, color_by = "mygs_Inf_response_part", smooth = TRUE, smooth_span = 0.2, pt_size = 2.0) 
# p13 <- plotSurface(spata_obj, color_by = "mygs_TNFA_part", smooth = TRUE, smooth_span = 0.2, pt_size = 2.0) 

# mean_13 signatures
# m <- plotSurfaceComparison(object = spata_obj, 
#                       color_by = c("mygs_ALA_pos_signature_251", "mygs_GPM_whole", "mygs_Injury_whole", "mygs_Hypoxia",
#                                    "mygs_Endothelial_marker", "mygs_Mesenchymal_whole", "mygs_Inf_response_whole",
#                                    "mygs_TNFA_whole", "mygs_GPM_part", "mygs_Injury_part", "mygs_Mesenchymal_part",
#                                    "mygs_Inf_response_part", "mygs_TNFA_part"),
#                       method_gs = "mean",
#                       smooth = TRUE, 
#                       smooth_span = 0.2, 
#                       pt_size = 2, 
#                       pt_clrsp = "inferno")
# 
# ggsave(plot = m, file=file.path("H:\\lab data\\Glioblastoma\\SPATA\\CNV\\New analysis/",
#                                 "334_Ti_surface_scaled.pdf"),
#        width = 12, height = 10)

# ssgsea_13 signatures
# s <- plotSurfaceComparison(object = spata_obj, 
#                       color_by = c("mygs_ALA_pos_signature_251", "mygs_GPM_whole", "mygs_Injury_whole", "mygs_Hypoxia",
#                                    "mygs_Endothelial_marker", "mygs_Mesenchymal_whole", "mygs_Inf_response_whole",
#                                    "mygs_TNFA_whole", "mygs_GPM_part", "mygs_Injury_part", "mygs_Mesenchymal_part",
#                                    "mygs_Inf_response_part", "mygs_TNFA_part"),
#                       method_gs = "ssgsea",
#                       smooth = TRUE, 
#                       smooth_span = 0.2, 
#                       pt_size = 2, 
#                       pt_clrsp = "inferno")

# ssgsea_5 signatures
# s <- plotSurfaceComparison(object = spata_obj, 
#                            color_by = c("mygs_ALA", "mygs_GPM", "mygs_Hypoxia",
#                                         "mygs_Endo", "mygs_MTC"),
#                            method_gs = "ssgsea",
#                            smooth = TRUE, 
#                            smooth_span = 0.2, 
#                            pt_size = 2, 
#                            pt_clrsp = "inferno")
# 
# ggsave(plot = s, file=file.path("H:\\lab data\\Glioblastoma\\SPATA\\CNV\\New analysis_5 sig/",
#                                 "243_Ti_ssgsea_denoised.pdf"),
#        width = 12, height = 10)


# z-score_5 signatures
s <- plotSurfaceComparison(object = spata_obj, 
                           color_by = c("mygs_IWR", "mygs_MES"),
                           method_gs = "zscore",
                           smooth = T,
                           pt_clrsp = "inferno",
                           smooth_span = 0.2,
                           pt_size = 2, 
                           # pt_alpha = 0.4,
                           # display_image = T
                           )

# s
# 
# validColorPalettes()
# validColorSpectra()


ggsave(plot = s, file=file.path("H:\\lab data\\Glioblastoma\\SPATA\\CNV\\New analysis_5 sig/",
                                "334_Ti_zscore_denoised_IWR_MES.pdf"),
       width = 12, height = 10)

# gsva
# v <- plotSurfaceComparison(object = spata_obj, 
#                            color_by = c("mygs_ALA_pos_signature_251", "mygs_GPM_whole", "mygs_Injury_whole", "mygs_Hypoxia",
#                                         "mygs_Endothelial_marker", "mygs_Mesenchymal_whole", "mygs_Inf_response_whole",
#                                         "mygs_TNFA_whole", "mygs_GPM_part", "mygs_Injury_part", "mygs_Mesenchymal_part",
#                                         "mygs_Inf_response_part", "mygs_TNFA_part"),
#                            method_gs = "gsva",
#                            smooth = TRUE, 
#                            smooth_span = 0.2, 
#                            pt_size = 2, 
#                            pt_clrsp = "inferno")
# 
# ggsave(plot = v, file=file.path("H:\\lab data\\Glioblastoma\\SPATA\\CNV\\New analysis/",
#                                 "334_Ti_gsva_denoised.pdf"),
#        width = 12, height = 10)


# library(gridExtra)
# grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, ncol=4)

# rm(joined_df)


# Image output
writeImage(spata_obj@images$`334_T`, "334_T_IHC.jpeg", quality = 100)
display(spata_obj@images$`334_T`)




############## ALA+ surface segment plot ##############

library(readxl)
library(ggplot2)

int_genes <- c("CD44", "PECAM1", "RBFOX3", "VEGFA", "NF1")

# compare gene expression on the surface
q <- plotSurfaceComparison(object = spata_obj, 
                      color_by = int_genes,
                      smooth = TRUE, 
                      method_gs = "ssgsea",
                      smooth_span = 0.2, 
                      pt_size = 2, 
                      pt_clrsp = "inferno")

ggsave(plot = q, file=file.path("H:\\lab data\\Glioblastoma\\SPATA\\CNV",
                                "334_Ti_genes_denoised.pdf"),
       width = 8.60, height = 5.57)

#######################################################

# open application to obtain a list of plots
plots <- plotSurfaceInteractive(object = spata_obj)

lapply(1:length(geneList), function(j){
  
  mygenes <- geneList[[j]]
  if(length(mygenes) > 25) mygenes <- mygenes[1:25]
  
  p <- plotSurfaceComparison(object = spata_obj, 
                        color_by = mygenes,
                        smooth = TRUE, 
                        smooth_span = 0.2, 
                        pt_size = 0.6, 
                        pt_clrsp = "inferno")
  ggsave(plot = p, file=file.path(outDir, paste0(names(geneList)[j], "_surface.pdf")),
         width = 15, height = 15)
  
})

########################
# IDENTIFY ALA5pos CELLS

# Use graphical interface to define segmentation
spata_obj <- createSegmentation(object = spata_obj)

# saveRDS(spata_obj, "243_Ti_AD_DE_Segment.RDS")

p <- plotSegmentation(object = spata_obj, pt_size = 2.0)
ggsave(plot = p, file=file.path("H:\\lab data\\Glioblastoma\\SPATA\\CNV\\New analysis_5 sig\\334_T_segmentation.pdf"),
       width = 6, height = 5)

ala5.features <- getFeatureVariables(spata_obj, features = "segmentation", return = "data.frame")
segment_df <- getSegmentDf(spata_obj, segment_names = c("Core_ALA", "Inf_ALA")) ### 334_T


########################################
# DIFFERENTIAL ANALYSIS ALA5pos vs. REST

seu <- plotSurface(object = spata_obj, color_by = "seurat_clusters", pt_size = 2)
ggsave(plot = seu, file=file.path("H:\\lab data\\Glioblastoma\\SPATA\\CNV", "334_T_seurat.pdf"),
       width = 7, height = 7)

# TO DO
spata_obj <- 
  runDeAnalysis(object = spata_obj,
                across = "segmentation", # across which identity groups
                method_de = c("wilcox") # with which methods
  )

# get an overview about the de-analysis results stored in your object
printDeaOverview(spata_obj)


getDeaResultsDf(object = spata_obj, 
                across = "segmentation", 
                method_de = "wilcox", 
                max_adj_pval = 0.05, # of every cluser take genes with an adj. p-value of 0.025 or lower
                n_lowest_pval = 100) # from the remaining genes take only to 'top 20')


# HEATMAP
cluster_of_interest <- c("Core_ALA", "Inf_ALA") ### 243_T

p <- plotDeaHeatmap(object = spata_obj, 
               across = "segmentation", # the grouping variable
               across_subset = cluster_of_interest, # the identity groups of interest 
               method_de = "wilcox", # the method with which the results were computed
               max_adj_pval = 0.05, # the adjusted p-value threshold
               n_lowest_pval = 20, 
               n_highest_lfc = 20, 
               clrp = "jama", 
               fontsize = 8) 
ggsave(plot = p, file=file.path(outDir, "334_Ti_segmentation_wilcoxon_heatmap.pdf"),
       width = 7, height = 7)

# OTHER PLOTS
high_threshold_markers <- 
  getDeaGenes(object = spata_obj,
              across = "segmentation", 
              across_subset = cluster_of_interest, 
              method = "wilcox", 
              n_lowest_pval = 3
  ) # return only the two genes with the lowest p-value for every identity group

# print genes
high_threshold_markers

# plot
p <- plotBoxplot(object = spata_obj, 
            variables = high_threshold_markers, 
            across = "segmentation", 
            across_subset = cluster_of_interest, 
            clrp = "jama", 
            nrow = 3
) +   legendTop()

ggsave(plot = p, file=file.path(outDir, "334_T_segmentation_wilcoxon_boxplot.pdf"),
       width = 10, height = 7)


### GSEA of segmentation

spata_obj <- 
  runGSEA(object = spata_obj, 
          across = "segmentation",
          methods_de = "wilcox",
          gene_set_names = c("mygs_ALA", "mygs_GPM", "mygs_IWR",
                             "mygs_Hypoxia", "mygs_MES", "mygs_MTC"))

spata_obj <- 
  runGSEA(object = spata_obj, 
          across = "segmentation",
          methods_de = "wilcox")


mygs_gsea <- getGseaDf(
  object = spata_obj, 
  across = "segmentation",
  method_de = "wilcox", 
  n_gsets = 100 # extract top 20 most significant gene sets
) 

write.csv(mygs_gsea, "GSEA_334_T_segmentation.csv")


p <- plotGseaDotPlot(
  object = spata_obj,
  across = "segmentation",
  # across_subset = as.character(seq(0,5,1)),
  n_gsets = 20,
  pt_alpha = 0.8,
  transform_with = list("fdr" = c("log10")),
  by_group = FALSE # merge in one plot
) 

ggsave(plot = p, file=file.path("H:\\lab data\\Glioblastoma\\SPATA\\CNV\\334_T_segmentation_GSEA.pdf"),
       width = 11, height = 7)

#########################################################
# TRAJECTORIES
# open interactive application
spata_obj <- createTrajectories(object = spata_obj)

# get trajectory names 
getTrajectoryNames(object = spata_obj)

plotTrajectory(object = spata_obj, 
               trajectory_name = "5ALA_evolution",
               color_by = "seurat_clusters",
               pt_clrp = "npg",
               pt_alpha = 0.25, # reduce alpha to highlight the trajectory's course
               pt_alpha2 = 1,
               display_image = FALSE) +
  legendTop()

plotTrajectory(object = spata_obj, 
               trajectory_name = "5ALA_evolution",
               color_by = "nCount_Spatial",
               smooth = TRUE, 
               pt_alpha = 0.25, 
               display_image = FALSE) +
  legendTop()

plotTrajectoryFeatures(object = spata_obj,
                       trajectory_name = "5ALA_evolution",
                       features = "nCount_Spatial", 
                       smooth_method = "loess", 
                       smooth_span = 0.2, 
                       smooth_se = TRUE) 


plotTrajectoryFeaturesDiscrete(object = spata_obj,
                               trajectory_name = "5ALA_evolution",
                               discrete_feature = "seurat_clusters", 
                               clrp = "npg",
                               display_trajectory_parts = FALSE) 


# gene-set names
genes_of_interest <- c("CD44", "NF1", "PECAM1", "RBFOX3", "VEGFA")

# plot lineplot
plotTrajectoryGenes(object = spata_obj,
                    trajectory_name = "5ALA_evolution", 
                    genes = genes_of_interest,
                    smooth_span = 0.2,
                    smooth_se = TRUE, 
                    display_facets = TRUE, # use facet_wrap() to split the plot in four parts
                    nrow = 2 # align the sub plots in two rows 
)


plotTrajectoryGeneSets(
  object = spata_obj,
  trajectory_name = "5ALA_evolution",
  gene_sets = "mygs_GPM",
  display_trajectory_parts = FALSE) + # results in missing vertical lines 
  legendTop()

plotTrajectoryGeneSets(
  object = spata_obj,
  trajectory_name = "5ALA_evolution",
  gene_sets = "mygs_IWR",
  display_trajectory_parts = FALSE) + # results in missing vertical lines 
  legendTop()

plotTrajectoryGeneSets(
  object = spata_obj,
  trajectory_name = "5ALA_evolution",
  gene_sets = "mygs_MES",
  display_trajectory_parts = FALSE) + # results in missing vertical lines 
  legendTop()

plotTrajectoryGeneSets(
  object = spata_obj,
  trajectory_name = "5ALA_evolution",
  gene_sets = "mygs_MTC",
  display_trajectory_parts = FALSE) + # results in missing vertical lines 
  legendTop()


all_genes <- getGenes(spata_obj)
all_gene_sets <- getGeneSets(spata_obj)

# obtain an assessed trajectory data.frame for all genes
atdf_genes <- assessTrajectoryTrends(object = spata_obj,
                                     trajectory_name = "5ALA_evolution",
                                     variables = all_genes)

# obtain an assessed trajectory data.frame for all gene-sets
atdf_gene_sets <- assessTrajectoryTrends(object = spata_obj,
                                         trajectory_name = "5ALA_evolution",
                                         variables = all_gene_sets)

# output
atdf_genes

# compare the trend of a variable to different models
plotTrajectoryFit(object = spata_obj,
                  trajectory_name = "5ALA_evolution",
                  variable = "mygs_MES",
                  display_residuals = TRUE) +
  legendTop()

# example 1: extract all variables that follow the linear descending trend while moving towards the hypoxic area (See Figure 2.1)
# with an auc-evaluation equal to or lower than 2
descending_genes <-
  filterTrajectoryTrends(atdf = atdf_genes,
                         limit = 3,
                         trends = "One peak", 
                         variables_only = FALSE) # return a data.frame

descending_genes


### customized trajectory

trajectory_length <- getTrajectoryLength(spata_obj, trajectory_name = "5ALA_evolution", binwidth = 5)

trajectory_length


trajectory_direction <- 1:trajectory_length

linear_ascending <- scales::rescale(1:trajectory_length, to = c(0,1))

plotTrajectoryGenes(object = spata_obj, trajectory_name = "5ALA_evolution", genes = "CD44") + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = trajectory_direction, y = linear_ascending),
                     color = "blue",
                     size = 1)

traj_df <- getTrajectoryDf(object = spata_obj,
                           trajectory_name = "5ALA_evolution",
                           variables = c("mygs_ALA",
                                         "mygs_GPM",
                                         "mygs_MES",
                                         "mygs_IWR",
                                         "mygs_MTC",
                                         "mygs_Hypoxia",
                                         "mygs_Endo",
                                         "CD44"), # CD44 gene 
                           binwidth = 5, 
                           shift_wider = TRUE)

dplyr::select(traj_df, trajectory_part, trajectory_order, mygs_ALA) 

atdf_cust <- assessTrajectoryTrendsCustomized(object = spata_obj,
                                              trajectory_name = "5ALA_evolution",
                                              customized_trends = dplyr::select(traj_df, mygs_ALA, mygs_GPM,
                                                                                mygs_MES, mygs_IWR,
                                                                                mygs_MTC, mygs_Hypoxia,
                                                                                mygs_Endo, CD44), 
                                              variables = all_genes
)

atdf_cust


similar_to_mygs_ALA <- filterTrajectoryTrends(atdf_cust, limit = 0.8, trends = "mygs_GPM")

# print names of similar genes
similar_to_mygs_ALA

# plot both in comparison 
plotTrajectoryGeneSets(object = spata_obj, 
                       trajectory_name = "5ALA_evolution", 
                       gene_sets = c("mygs_GPM"), 
                       display_facets = TRUE, 
                       clrp = "default")

plotTrajectoryGenes(object = spata_obj, 
                    trajectory_name = "5ALA_evolution", 
                    genes = similar_to_mygs_ALA[1:4], 
                    display_facets = TRUE)

similar_to_CD44 <- filterTrajectoryTrends(atdf_cust, limit = 2, trends = "CD44")

# print names of similar gene sets
similar_to_CD44

# plot both in comparison 
plotTrajectoryGenes(object = spata_obj, 
                    trajectory_name = "5ALA_evolution", 
                    genes = c("CD44"), 
                    display_facets = TRUE, 
                    clrp = "default")

plotTrajectoryGenes(object = spata_obj, 
                    trajectory_name = "5ALA_evolution", 
                    genes = similar_to_CD44, 
                    display_facets = TRUE)


descending_genes_vec <- descending_genes$variables

hm_colors <- viridis::inferno(n = 100)

plotTrajectoryHeatmap(object = spata_obj, 
                      trajectory_name = "5ALA_evolution",
                      variables = descending_genes_vec,
                      arrange_rows = "maxima",
                      colors = hm_colors,
                      show_rownames = TRUE,
                      split_columns = FALSE, 
                      smooth_span = 0.5)


######## END of TRAJECTORY analysis ########

######################################################

# 1. get a data.frame that contains barcodes variables
coords_df <- getCoordsDf(spata_obj)

# 2. join this data.frame with additional information
joined_df <- 
  joinWith(object = spata_obj, 
           spata_df = coords_df,
           # genes = geneList$ALA_pos_signature_251, # expression values of the gene set Hallmark-Hypoxia
           features = "seurat_clusters", # cluster belonging
           verbose = FALSE)

# output 
joined_df

write.csv(joined_df, "334_Ti_coordinates.csv")

#################
# CHECK SPATA OBJ
plotSurface(object = spata_obj, color_by = "seurat_clusters", pt_size = 2)
plotSurface(object = spata_obj, color_by = "nCount_Spatial", smooth_span = 0.1, pt_size = 1)

# segment the "good quality area" as well as the "bad quality area"
spata_obj2 <- createSegmentation(object = spata_obj)


# display the current segmentation
#plotSegmentation(object = spata_obj, pt_size = 1.9)


monocle_clusters <- findMonocleClusters(object = spata_obj, 
                                        preprocess_method = "PCA", 
                                        reduction_method = c("UMAP", "PCA", "tSNE"), 
                                        cluster_method = c("leiden", "louvain"), 
                                        k = 5, 
                                        num_iter = 5)

# output
monocle_clusters
# output
examineClusterResults(data = monocle_clusters)
# add the cluster results
spata_obj <- 
  addFeatures(object = spata_obj, 
              feature_names = c("cluster_leiden_UMAP_k5", "cluster_leiden_tSNE_k5","cluster_louvain_PCA_k5"), 
              feature_df = monocle_clusters,
              overwrite = TRUE,
              key = "barcodes")

# feature names afterwards
getFeatureNames(spata_obj)

plotSurface(object = spata_obj,
            color_by = "cluster_leiden_UMAP_k5",
            pt_size = 1.9,
            pt_clrp = "jama") +
  ggplot2::labs(color = "Leiden UMAP")


genes_b <- c("CARTPT", "OLIG1", "GFAP", "SYNPR", "HOPX", "CCK")
genes_b <- as.character(genes_b)
library(dplyr)
names(plots)
plots <- plotSurfaceInteractive(object = spata_obj)
my_df <- 
  joinWith(object = spata_obj, gene_sets = geneList$Hypoxia, method_gs = "mean", smooth = FALSE, normalize = TRUE) %>% 
  dplyr::mutate(
    new_grouping = dplyr::case_when(
      HM_HYPOXIA< 0.25 ~ "group_a", 
      HM_HYPOXIA< 0.5 ~ "group_b", 
      HM_HYPOXIA< 0.75 ~ "group_c", 
      HM_HYPOXIA> 0.75 ~ "group_d"
    ) %>% as.factor()
  )

object <- addFeatures(object= spata_obj, feature_df = my_df, feature_names = "new_grouping")



####################
### DEA analysis ###
####################

plotSurface(object = spata_obj, 
            color_by = "seurat_clusters", 
            # pt_clrp = "jama", 
            pt_size = 2.0
)

spata_obj <- 
  runDeAnalysis(object = spata_obj,
                across = "seurat_clusters", # across which identity groups
                method_de = c("wilcox", "bimod") # with which methods
  )

# saveRDS(spata_obj, "334_Ti_AD_DE.RDS")

gsea_resut <- getDeaResultsDf(object = spata_obj, 
                              across = "seurat_clusters", 
                              method_de = "wilcox")
write.csv(gsea_resut, "243_T_Top_genes.csv")
cluster_of_interest <- as.character(seq(0,8,1))

plotDeaHeatmap(object = spata_obj, 
               across = "seurat_clusters", # the grouping variable
               across_subset = cluster_of_interest, # the identity groups of interest
               method_de = "wilcox", # the method with which the results were computed
               max_adj_pval = 0.05, # the adjusted p-value threshold
               n_lowest_pval = 20, 
               n_highest_lfc = 20, 
               clrp = "jama", 
               fontsize = 6) 

plotDeaDotPlot(
  object = spata_obj, 
  across = "seurat_clusters", 
  across_subset = cluster_of_interest,
  color_by = "avg_logFC",
  color_trans = "log10", # transform color scale with log10
  n_highest_lfc = 20,
  by_group = TRUE, # create separate plots for each specified group
  nrow = 2 # make sure that subplot are displayed in two rows
)

plotDeaDotPlot(
  object = spata_obj, 
  across = "seurat_clusters", 
  across_subset = cluster_of_interest,
  n_highest_lfc = 7,
  color_by = "avg_logFC",
  color_trans = "log10",
  by_group = FALSE # merge in one plot
)

genes_of_interest <-
  getDeaGenes(object = spata_obj,
              across = "seurat_clusters", # the grouping variable
              across_subset = cluster_of_interest, # the identity groups of interest
              method_de = "wilcox", # the method with which the results were computed
              max_adj_pval = 0.025, # the adjusted p-value threshold
              n_lowest_pval = 10, 
              n_highest_lfc = 10
  )


head(genes_of_interest) # first six


high_threshold_markers <- 
  getDeaGenes(object = spata_obj,
              across = "seurat_clusters", 
              # across_subset = cluster_of_interest, 
              method = "wilcox", 
              n_lowest_pval = 3
  ) # return only the two genes with the lowest p-value for every identity group

# print genes
high_threshold_markers

# visualize statistics
plotBoxplot(object = spata_obj, 
            variables = high_threshold_markers, 
            across = "seurat_clusters", 
            # across_subset = cluster_of_interest, 
            clrp = "jama", 
            nrow = 3
) + 
  legendTop()



plotSurfaceComparison(object = spata_obj, 
                      color_by = high_threshold_markers, 
                      smooth = TRUE, 
                      pt_size = 1
) + legendNone()

# compare to cluster localisation
plotSurface(object = spata_obj, 
            color_by = "seurat_clusters", 
            pt_clrp = "jama"
)

plotDeaGeneCount(object = spata_obj, across = "seurat_clusters", method_de = "wilcox", clrp = "jama")

plotDeaLogFC(object = spata_obj, across = "seurat_clusters", method_de = "wilcox", clrp = "jama")

plotDeaPvalues(object = spata_obj,
               across = "seurat_clusters",
               method_de = "wilcox",
               clrp = "jama",
               plot_type = "density", 
               scales = "free")

spata_obj <- 
  runGSEA(object = spata_obj, 
          across = "seurat_clusters",
          methods_de = "wilcox")


mygs_gsea <- getGseaDf(
  object = spata_obj, 
  across = "seurat_clusters",
  method_de = "wilcox", 
  n_gsets = 20 # extract top 20 most significant gene sets
) 

write.csv(mygs_gsea, "GSEA_334_T_all.csv")


plotGseaDotPlot(
  object = spata_obj,
  across = "seurat_clusters",
  # across_subset = as.character(seq(0,5,1)),
  n_gsets = 5,
  pt_alpha = 0.8,
  transform_with = list("fdr" = c("log10")),
  by_group = FALSE # merge in one plot
) 


####################

plotSurface(object = spata_obj, color_by = "seurat_clusters")

plotCnvResults(object = spata_obj, across = "seurat_clusters")

plotSurface(object = spata_obj, color_by = "Chr7", pt_clrsp = "Reds 3", c1 = 1)

plotSurface(object = spata_obj, color_by = "Chr10", pt_clrsp = "Blues 3", rev = FALSE, c1 = 1)

q(save = "no")



####################
### spatawrapper ###
####################

# remotes::install_github("heilandd/SPATAwrappers")
library(SPATAwrappers)
library(Seurat)

seuratObj <- SPATA2::transformSpataToSeurat(spata_obj)

plot.folder <- "H:\\lab data\\Glioblastoma\\SPATA\\CNV\\plot folder"

seuratObj <- 
  seuratObj %>% 
  SCTransform(assay="Spatial") %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:10)


seuratObj <- seuratObj %>% run.SNN.stability(assay="SCT", 
                                             reduction="pca",
                                             dims = 1:30, 
                                             resolution = seq(from=0.1, to=1.5, by=0.05), 
                                             cluster_id="Louvain", algorithm=1)

seuratObj <- seuratObj %>% run.SNN.stability(assay="SCT", reduction="pca",
                                             resolution = seq(from=0.1, to=1.5, by=0.05), 
                                             dims = 1:30, cluster_id="SLM",algorithm=3)



ggsave(Seurat::DimPlot(seuratObj, group.by = "Louvain"),
       filename=paste0("H:\\lab data\\Glioblastoma\\SPATA\\CNV/DimPlot_cluster_State1_Louvain.pdf"))

ggsave(Seurat::DimPlot(seuratObj, group.by = "SLM"),
       filename=paste0("H:\\lab data\\Glioblastoma\\SPATA\\CNV/DimPlot_cluster_State1_SLM.pdf"))


##Add data to matrix
#Create DF for all clusters

# PAM Cluster

k=length(unique(seuratObj@meta.data$Louvain))

pca <- seuratObj %>% Seurat::Reductions("pca")
pca <- pca@cell.embeddings %>% as.data.frame()
pca <- pca[,1:30]
pam_PCA <- cluster::pam(pca, k=k)

pca <- seuratObj %>% Seurat::Reductions("pca")
pca <- pca@cell.embeddings %>% as.data.frame()
pca <- pca[,1:30]
HC_PCA <- hclust(dist(pca))
HC_PCA <- cutree(HC_PCA, k=k)


cluster_summary=data.frame(barcodes=seuratObj@meta.data %>% rownames(), 
                           PCA_1=c(as.numeric(seuratObj@meta.data$Louvain)+1),
                           PCA_2=c(as.numeric(seuratObj@meta.data$SLM)+1),
                           PCA_4=pam_PCA$clustering %>% as.numeric(),
                           PCA_5=HC_PCA %>% as.numeric())

write.csv(cluster_summary, file=paste0(plot.folder, "/cluster_summary.csv"))

# ggsave(run.fscore(cluster_summary[-c(1)]),
#       filename=paste0(plot.folder, "/Fscore.png"))

# write.csv(run.fscore.mat(cluster_summary[-c(1)]), file=paste0(plot.folder, "/Fscore.csv"))


prefix = "PCA_"
out <- clustree::clustree(cluster_summary, 
                          prefix = prefix)

out.plot <- out$data %>% 
  select(!!sym(prefix), sc3_stability) %>% 
  group_by(!!sym(prefix)) %>% 
  summarise(mean(sc3_stability), sd(sc3_stability)) %>% 
  as.data.frame()


## Validate stabolity
best.cluster <- which.max(out.plot$`mean(sc3_stability)`)

cluster_summary_UMAP=data.frame(barcodes=seuratObj@meta.data %>% rownames(), 
                                PCA_1=c(as.numeric(seuratObj@meta.data$Louvain)+1),
                                PCA_2=c(as.numeric(seuratObj@meta.data$SLM)+1),
                                PCA_3=pam_PCA$clustering %>% as.numeric(),
                                UMAP1=seuratObj@reductions$umap@cell.embeddings[,1],
                                UMAP2=seuratObj@reductions$umap@cell.embeddings[,2])

ggsave(clustree::clustree_overlay(cluster_summary_UMAP, 
                                  prefix = prefix, 
                                  x_value = "UMAP1", 
                                  y_value = "UMAP2"),
       filename=paste0(plot.folder, "/ClusterTree_UMAP.pdf"))




# Merge Clusters together
stable <- 
  cluster_summary_UMAP %>% 
  mutate(stable = PCA_1+PCA_2+PCA_3 ) %>% 
  count(stable) %>% 
  filter(n>200) %>% 
  pull(stable)

# Select stable clusters
cluster_summary_UMAP <- 
  cluster_summary_UMAP %>% 
  mutate(stable = PCA_1+PCA_2+PCA_3 ) %>% 
  filter(stable %in% {{stable}})

plot <- clustree::clustree_overlay(cluster_summary_UMAP, 
                                   prefix = prefix, 
                                   x_value = "UMAP1", 
                                   y_value = "UMAP2")

ggsave(plot,
       filename=paste0(plot.folder, "/ClusterTree_UMAP_stable.pdf"))



for(i in 1:length(unique(cluster_summary_UMAP$stable))){cluster_summary_UMAP$stable[cluster_summary_UMAP$stable==unique(cluster_summary_UMAP$stable)[i]] <- i }

#ggplot(data=cluster_summary_UMAP, aes(x=UMAP1, y=UMAP2, color=as.factor(stable)))+geom_point()+theme_classic()

## Add consensus clusers to seurat

seuratObj <- subset(seuratObj, cells=cluster_summary_UMAP$barcodes)
seuratObj <- 
  seuratObj %>% 
  SCTransform(assay="Spatial") %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:10)

seuratObj@meta.data <- 
  seuratObj@meta.data[cluster_summary_UMAP$barcodes, ] %>% 
  mutate(consensus=cluster_summary_UMAP$stable)

ggsave(Seurat::DimPlot(seuratObj, reduction = "pca", group.by = "consensus" ),
       filename=paste0(plot.folder, "/consensus_PCA.pdf"))
ggsave(Seurat::DimPlot(seuratObj, reduction = "umap", group.by = "consensus" ),
       filename=paste0(plot.folder, "/consensus_UMAP.pdf"))


seuratObj <- SetIdent(seuratObj, value=seuratObj@meta.data$consensus)

markers <- FindAllMarkers(seuratObj)

final.clusters <- markers %>% count(cluster) %>% filter(n>50) %>% pull(cluster)
markers.final <- markers %>% filter(cluster %in% final.clusters)

# Remove duplicated genes 
markers.final$rows <- 1:nrow(markers.final)
dup <- markers.final[duplicated(markers.final$gene), ]$gene
remove <- purrr::map(.x=dup, .f=function(i){
  keep <- markers.final %>% 
    filter(gene=={{i}}) %>% 
    arrange(desc(avg_log2FC)) %>% 
    head(1) %>% 
    pull(rows)
  
  remove <- markers.final %>% 
    filter(gene=={{i}}) %>% 
    filter(rows!={{keep}}) %>% 
    pull(rows)
  
  return(remove)
  
}) %>% unlist()


markers.final <- markers.final[-c(remove), ]
final.clusters <- markers.final %>% count(cluster) %>% filter(n>20) %>% pull(cluster)
markers.final <- markers.final %>% filter(cluster %in% final.clusters)
rownames(markers.final) <- markers.final$gene


## Add new cluster to seurat
cluster <- as.numeric(unique(markers.final$cluster))
bc <- 
  seuratObj@meta.data[seuratObj@meta.data$consensus %in% cluster, ] %>% rownames()

seuratObj <- subset(seuratObj, cells=bc)
seuratObj <- 
  seuratObj %>% 
  SCTransform(assay="SCT") %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:10)

ggsave(Seurat::DimPlot(seuratObj, reduction = "pca", group.by = "consensus" ),
       filename=paste0(plot.folder, "/consensus_PCA_final.pdf"))
ggsave(Seurat::DimPlot(seuratObj, reduction = "umap", group.by = "consensus" ),
       filename=paste0(plot.folder, "/consensus_UMAP_final.pdf"))




p1 <- DimPlot(seuratObj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(seuratObj, label = TRUE, label.size = 3)
p1 + p2

spata_obj <- transformSeuratToSpata(seuratObj,sample_name ="334_T_seu",method = "spatial",
                                                        assay_name="Spatial")


# plots <- plotSurfaceInteractive(object = x)

plots <- plotSurfaceInteractive(object = spata_obj)

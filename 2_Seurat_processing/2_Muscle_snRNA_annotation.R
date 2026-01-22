setwd('/Volumes/BB_HQ_4/Collaboration/Cortes_lab/Muscle_snRNAseq/Seurat_Processing')
options(stringsAsFactors = F)
options (future.globals.maxSize = 32000 * 1024^2)

# Load packages
library('Seurat')    # 
library(sctransform) # 
library(clustree)    # 
library(scales)      # 
library(dplyr)      # 
library(readxl)
library(Polychrome)

##########  Cell identity annotation packages ##########  
# https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html
library(scSorter)

# https://github.com/TianLab-Bioinfo/scMAGIC
library(scMAGIC)
########################################################

#####################################################################################################################
# 2024-11-29
# annotate cells with:
#    - scSorter/scbase and scMayoMap
#    - scMAGIC
#####################################################################################################################


#####################################################################################################################
#### 1. Load Cleaned up Seurat Objects and merge data

############ load integrated object
load('2024-11-29_Muscle_Exercise_TFEB_Seurat_object_SCT_with_UMAP.RData')
muscle.clean
# An object of class Seurat 
# 21160 features across 209939 samples within 1 assay 
# Active assay: RNA (21160 features, 5000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony


############ Add clustering granularity to get all important cell types resolved
muscle.clean <- FindNeighbors(muscle.clean, dims = 1:18)  # number of PCs found to make sense
muscle.clean <- FindClusters(muscle.clean, resolution = c(0.6, 1.2, 1.5))
# Number of communities: 25
# Number of communities: 33
# Number of communities: 40


pdf(paste(Sys.Date(),"Muscle_Exercise_TFEB_Singlets_UMAP_snn_1.2.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(muscle.clean, reduction = "umap", group.by = "SCT_snn_res.1.2", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Muscle_Exercise_TFEB_Singlets_UMAP_snn_1.5.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(muscle.clean, reduction = "umap", group.by = "SCT_snn_res.1.5", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Muscle_Exercise_TFEB_Singlets_UMAP_snn_1.2_LABEL.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(muscle.clean, reduction = "umap", group.by = "SCT_snn_res.1.2", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024), label = T)
dev.off()

pdf(paste(Sys.Date(),"Muscle_Exercise_TFEB_Singlets_UMAP_snn_1.5_LABEL.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(muscle.clean, reduction = "umap", group.by = "SCT_snn_res.1.5", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024), label = T)
dev.off()

pdf(paste(Sys.Date(),"Muscle_Exercise_TFEB_Singlets_UMAP_snn_0.6_LABEL.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(muscle.clean, reduction = "umap", group.by = "SCT_snn_res.0.6", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024), label = T)
dev.off()


pdf(paste(Sys.Date(),"Muscle_RidgePlot_XY_linked_genes.pdf", sep = "_"), height = 6, width = 10)
RidgePlot(muscle.clean, features = c("Ddx3y","Xist"), group.by = "Group", cols = c(alpha("lightpink2", alpha = 0.5    ), alpha("deeppink", alpha = 0.5    ), alpha("mediumpurple4", alpha = 0.5    ) ,
                                                                                   alpha("mediumturquoise", alpha = 0.5 ), alpha("deepskyblue", alpha = 0.5 ), alpha("royalblue3", alpha = 0.5 ) ))
dev.off()


################################################################################################################################################################



################################################################################################################################################################
##### 2. Try using scSorter "gates" to annotate cell types
# https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html

# Use mouse markers (singleCellBase)


####################################################################################################################################################
#####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Mouse markers from singleCellBase      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#####
####################################################################################################################################################


# ###%%%%%%%%%%%%%%%%%%%### A. Filter and prep singleCellBase marker data  ###%%%%%%%%%%%%%%%%%###

############ Section 1 - Generate annotation object for scSorter ############
# Read in marker file
scbase.markers <- read.csv("../Reference/Markers/2024-11-29_singleCellBase_20230904_mouse.txt", sep = "\t", header = T)
scbase.markers <- scbase.markers[scbase.markers$tissue_type %in% c("Gastrocnemius muscle","Limb muscle", "Muscle"),]

unique(scbase.markers$cell_type)
# [1] "Tumor cells"                    "Dendritic cells"                "Endothelial cells"              "Fibroblasts"                   
# [5] "Muscle cells"                   "Myeloid cells"                  "Myofibroblasts"                 "Neutrophils"                   
# [9] "Plasmacytoid dentritic cells"   "T cells"                        "Stem cells"                     "Satellite cells"               
# [13] "B cells"                        "Ec"                             "Faps"                           "Mph"                           
# [17] "Myonuclei"                      "Myonucler"                      "Pro-inflammatory maccrophages"  "Prolif ic"                     
# [21] "Regmyon"                        "Tc"                             "Tenocytes"                      "Progenitor cells"              
# [25] "Glial cells"                    "Macrophages"                    "Myoblasts"                      "Antigen-presenting cells (APC)"
# [29] "Monocytes"                      "Pericytes"                      "Schwann cells"       


# reduce cell type redundancy
scbase.markers$cell_type[scbase.markers$cell_type %in% "Ec"]                             <- "Endothelial cells"
scbase.markers$cell_type[scbase.markers$cell_type %in% "Tc"]                             <- "T cells"
scbase.markers$cell_type[scbase.markers$cell_type %in% "Faps"]                           <- "FAPs"
scbase.markers$cell_type[scbase.markers$cell_type %in% "Progenitor cells"]               <- "FAPs"
scbase.markers$cell_type[scbase.markers$cell_type %in% "Stem cells"]                     <- "Satellite cells"
scbase.markers$cell_type[scbase.markers$cell_type %in% "Myonucler"]                      <- "Myonuclei"
scbase.markers$cell_type[scbase.markers$cell_type %in% "Muscle cells"]                   <- "Myonuclei"
scbase.markers$cell_type[scbase.markers$cell_type %in% "Mph"]                            <- "Macrophages"
scbase.markers$cell_type[scbase.markers$cell_type %in% "Pro-inflammatory maccrophages"]  <- "Macrophages"
scbase.markers$cell_type[scbase.markers$cell_type %in% "Plasmacytoid dentritic cells"]   <- "Dendritic cells"

scbase.markers$cell_type[scbase.markers$cell_subtype %in% "Faps"]  <- "FAPs"
scbase.markers$cell_type[scbase.markers$cell_subtype %in% "Muscle stem cells"]  <- "Satellite cells"
scbase.markers$cell_type[scbase.markers$cell_subtype %in% "Monocytes/macrophages"]  <- "Macrophages"

scbase.cell.types <- unique(scbase.markers$cell_type)
scbase.cell.types <- setdiff(scbase.cell.types, c("Tumor cells","Prolif ic","Myeloid cells"))
scbase.cell.types
# [1] "Tumor cells"                    "Dendritic cells"                "Endothelial cells"              "Fibroblasts"                   
# [5] "Muscle cells"                   "Myeloid cells"                  "Myofibroblasts"                 "Neutrophils"                   
# [9] "Plasmacytoid dentritic cells"   "T cells"                        "Stem cells"                     "Satellite cells"               
# [13] "B cells"                        "Ec"                             "Faps"                           "Mph"                           
# [17] "Myonuclei"                      "Myonucler"                      "Pro-inflammatory maccrophages"  "Prolif ic"                     
# [21] "Regmyon"                        "Tc"                             "Tenocytes"                      "Progenitor cells"              
# [25] "Glial cells"                    "Macrophages"                    "Myoblasts"                      "Antigen-presenting cells (APC)"
# [29] "Monocytes"                      "Pericytes"                      "Schwann cells"       


# create marker list object to summarize all marker genes
my.marker.list         <- vector(mode = "list",length = length(scbase.cell.types))
names(my.marker.list)  <- scbase.cell.types

for (i in 1:length(scbase.cell.types)) {
  
  # get unique markers from singleCellBase
  
  my.marker.list[[i]] <- unique( unlist(strsplit(scbase.markers$gene_symbol[scbase.markers$cell_type %in% scbase.cell.types[i]], ", ")))
  
}

names(my.marker.list)

# overlap with genes detected in our data
my.detected <- rownames(muscle.clean)
my.marker.list.v2 <- vector(mode = "list",length = length(scbase.cell.types))

for (i in 1:length(my.marker.list)) {
  my.marker.list.v2[[i]] <- intersect(my.marker.list[[i]],my.detected)
}
names(my.marker.list.v2) <- names(my.marker.list)

# check marker lists
lapply(my.marker.list.v2, length)


### Generate anno table
# initialize
anno           <- data.frame(cbind(rep(names(my.marker.list.v2)[1],length(my.marker.list.v2[[1]])), my.marker.list.v2[[1]], rep(2,length(my.marker.list.v2[[1]]))))
colnames(anno) <- c("Type", "Marker", "Weight")

for (i in 2:length(my.marker.list.v2)) {
  anno <- rbind(anno,
                cbind("Type" = rep(names(my.marker.list.v2)[i],length(my.marker.list.v2[[i]])),
                      "Marker" = my.marker.list.v2[[i]],
                      "Weight" = rep(2,length(my.marker.list.v2[[i]]))))
  
}

############ Section 2 - Pre-processing the data for scSorter ############

# get highly variable genes and filter out genes with non-zero expression in less than 5% of total cells.
topgenes        <- head(VariableFeatures(muscle.clean), 2500)
expr            <- GetAssayData(muscle.clean)
topgene_filter  <- rowSums(expr[topgenes, ]!=0) > ncol(expr)*.05
topgenes        <- topgenes[topgene_filter]

# At last, we subset the preprocessed expression data.
# Keep all genes that can be markers, even if not for cell types in that tissue
picked_genes    <- unique(c(anno$Marker, topgenes))
print(length(picked_genes)) # anno  ### 697
expr.sctype     <- expr[rownames(expr) %in% picked_genes, ]

# Now, we are ready to run scSorter.


############ Section 3 - Running scSorter ############
# run scSorter
rts.muscle  <- scSorter(expr.sctype , anno)
save(rts.muscle, file = paste0(Sys.Date(),"_scSorter_out_Muscle_singleCellBase.RData"))

############ Section 4 - gate ScSorter calls back to Seurat  ############
# get calls from scSorter calls
muscle.clean@meta.data$ScSorter_scBase <- "other" # initialize
muscle.clean@meta.data[colnames(expr.sctype ), ]$ScSorter_scBase <- rts.muscle$Pred_Type


pdf(paste(Sys.Date(),"Muscle_UMAP_color_by_ScSorter_singleCellBase_call.pdf", sep = "_"), height = 5, width = 10)
DimPlot(muscle.clean, reduction = "umap", group.by = "ScSorter_scBase", raster = F)
dev.off()

pdf(paste(Sys.Date(),"Muscle_UMAP_color_by_ScSorter_singleCellBase_call_RASTER.pdf", sep = "_"), height = 5, width = 10)
DimPlot(muscle.clean, reduction = "umap", group.by = "ScSorter_scBase", raster = T, raster.dpi = c(600,600))
dev.off()

# save annotated object
save(muscle.clean, file = paste0(Sys.Date(),"_Seurat_object_with_ScSorter_singleCellBase_Raw_output.RData"))
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###


# ###%%%%%%%%%%%%%%%%%%%### B. Filter and prep scMayoMap marker data  ###%%%%%%%%%%%%%%%%%###

############ Section 1 - Generate annotation object for scSorter ############
# Read in marker file
scMayoMap.markers <- readxl::read_xlsx("../Reference/Markers/scMayoMap_12915_2023_1728_MOESM1_ESM.xlsx", skip = 1)
scMayoMap.muscle  <- data.frame(scMayoMap.markers[scMayoMap.markers$tissue %in% "muscle",])

unique(scMayoMap.muscle$celltype)
# [1] "Adventitial cell"                     "B cell"                               "Endothelial cell"                    
# [4] "Erythroblast"                         "Fibroblast"                           "Granulocyte monocyte progenitor cell"
# [7] "Inflammatory cell"                    "Macrophage"                           "Mesenchymal progenitor cell"         
# [10] "Mesenchymal stem cell"                "Muscle-derived cell"                  "Myoblast"                            
# [13] "Myocyte"                              "Myocyte progenitor cell"              "Myoepithelial cell"                  
# [16] "Myofibroblast"                        "Myogenic endothelial cell"            "Neutrophil"                          
# [19] "Pericyte"                             "Progenitor cell"                      "Satellite cell"                      
# [22] "Schwann cell"                         "Smooth muscle cell"                   "T cell"                              
# [25] "Tenocyte"    

####################################
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# create marker list object to summarize all marker genes
scMayoMap.cell.types    <- unique(scMayoMap.muscle$celltype)
my.marker.list         <- vector(mode = "list",length = length(scMayoMap.cell.types))
names(my.marker.list)  <- scMayoMap.cell.types

for (i in 1:length(scMayoMap.cell.types)) {
  
  ix.c <- which (scMayoMap.muscle$celltype %in% scMayoMap.cell.types[i])
  
  # get unique markers from singleCellBase
  my.marker.list[[i]] <- unique(firstup(tolower(as.character(na.omit(as.character(scMayoMap.muscle[ix.c,-c(1:2)]))))))
  
}


names(my.marker.list)

# overlap with genes detected in our data
my.detected <- rownames(muscle.clean)
my.marker.list.v2 <- vector(mode = "list",length = length(scbase.cell.types))

for (i in 1:length(my.marker.list)) {
  my.marker.list.v2[[i]] <- intersect(my.marker.list[[i]],my.detected)
}
names(my.marker.list.v2) <- names(my.marker.list)

# check marker lists
lapply(my.marker.list.v2, length)


### Generate anno table
# initialize
anno           <- data.frame(cbind(rep(names(my.marker.list.v2)[1],length(my.marker.list.v2[[1]])), my.marker.list.v2[[1]], rep(2,length(my.marker.list.v2[[1]]))))
colnames(anno) <- c("Type", "Marker", "Weight")

for (i in 2:length(my.marker.list.v2)) {
  anno <- rbind(anno,
                cbind("Type" = rep(names(my.marker.list.v2)[i],length(my.marker.list.v2[[i]])),
                      "Marker" = my.marker.list.v2[[i]],
                      "Weight" = rep(2,length(my.marker.list.v2[[i]]))))
  
}

############ Section 2 - Pre-processing the data for scSorter ############

# get highly variable genes and filter out genes with non-zero expression in less than 5% of total cells.
topgenes        <- head(VariableFeatures(muscle.clean), 2500)
expr            <- GetAssayData(muscle.clean)
topgene_filter  <- rowSums(expr[topgenes, ]!=0) > ncol(expr)*.05
topgenes        <- topgenes[topgene_filter]

# At last, we subset the preprocessed expression data.
# Keep all genes that can be markers, even if not for cell types in that tissue
picked_genes    <- unique(c(anno$Marker, topgenes))
print(length(picked_genes)) # anno  ### 1004
expr.sctype     <- expr[rownames(expr) %in% picked_genes, ]

# Now, we are ready to run scSorter.


############ Section 3 - Running scSorter ############
# run scSorter
rts.muscle  <- scSorter(expr.sctype , anno)
save(rts.muscle, file = paste0(Sys.Date(),"_scSorter_out_Muscle_scMayoMap.RData"))

############ Section 4 - gate ScSorter calls back to Seurat  ############
# get calls from scSorter calls
muscle.clean@meta.data$ScSorter_scMayoMap <- "other" # initialize
muscle.clean@meta.data[colnames(expr.sctype ), ]$ScSorter_scMayoMap <- rts.muscle$Pred_Type

pdf(paste(Sys.Date(),"Muscle_UMAP_color_by_ScSorter_scMayoMap_call.pdf", sep = "_"), height = 5, width = 15)
DimPlot(muscle.clean, reduction = "umap", group.by = "ScSorter_scMayoMap", raster = F)
dev.off()

pdf(paste(Sys.Date(),"Muscle_UMAP_color_by_ScSorter_scMayoMap_call_RASTER.pdf", sep = "_"), height = 5, width = 15)
DimPlot(muscle.clean, reduction = "umap", group.by = "ScSorter_scMayoMap", raster = T, raster.dpi = c(600,600))
dev.off()

# save annotated object
save(muscle.clean, file = paste0(Sys.Date(),"_Seurat_object_with_ScSorter_singleCellBase_scMayoMap_Raw_output.RData"))
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###

###############################################################################################################################################################


###############################################################################################################################################################
#### 3. Try scMagic
# https://github.com/TianLab-Bioinfo/scMAGIC

# ###%%%%%%%%%%%%%%###  0. Before running scMAGIC, please firstly run the following codes  ###%%%%%%%%%%%%%%###
library(reticulate)
py_config()  # your python environment
print(py_module_available('numpy')) # whether the "numpy" has been installed
np      <- import("numpy")
np.exp2 <- np$exp2
np.max  <- np$max
np.sum  <- np$sum

# For Windows users, please set "method_HVGene = 'SciBet_R'"
#### Also for MacOSX !!!!
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###

# ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%### A. Filter and prep Cosgrove training data  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
# Downloaded data from dryad webiste
load("/Volumes/BB_HQ_4/Collaboration/Cortes_lab/Muscle_snRNAseq/Reference/Atlas/myo_slim_seurat_v1-1.RData")
myo.slim.seurat

table(myo.slim.seurat$prefilter_IDs)
# B Cells     Dendritic   Endothelial          FAPs   Macrophages     Monocytes         MuSCs     Myoblasts     Myonuclei        Neural 
# 29            17            32            82            11           881         39368          1923         41873            12 
# Neutrophils      NK Cells Smooth Muscle       T Cells     Tenocytes 
# 10            12            20           103            10 

# format reference dataset
ref.mtx     <- myo.slim.seurat@assays$RNA@counts
ref.labels  <- myo.slim.seurat$prefilter_IDs

# run scMAGIC
muscle.scMAGIC <- scMAGIC(exp_sc_mat         =   muscle.clean@assays$RNA@counts    ,
                          exp_ref_mat        =   ref.mtx         ,
                          exp_ref_label      =   ref.labels      ,
                          atlas              =   "MCA"            ,
                          method_HVGene      =   'SciBet_R'      ,
                          num_threads        =   1               ,
                          method_findmarker  =   'Seurat'        ,
                          cluster_num_pc     =   18              , # from PCA analysis
                          min_cell           =   10              ,
                          method1            =   'spearman'      ,
                          percent_high_exp   =    0.8            )

save(muscle.scMAGIC, file = paste0(Sys.Date(),"_scMAGIC_Cosgrove_MCA_RawOutput.RData"))

#######
# gate back to killi brain seurat object
muscle.clean@meta.data$scMAGIC_Cosgrove <- NA

# grab cell ids
my.cell.ids <- rownames(muscle.scMAGIC)
my.cell.ids <- my.cell.ids[-grep('NA',my.cell.ids)] # remove NAs

# populate meta-data
muscle.clean@meta.data[my.cell.ids,]$scMAGIC_Cosgrove <- as.vector(muscle.scMAGIC[my.cell.ids,])

# Plot Raw output
pdf(paste(Sys.Date(),"Muscle_UMAP_color_by_scMAGIC_Cosgrove_call_RASTER.pdf", sep = "_"), height = 5, width = 6)
DimPlot(muscle.clean, reduction = "umap", group.by = "scMAGIC_Cosgrove", raster = T, raster.dpi = c(600,600))
dev.off()

# Save object with raw annotations
save(muscle.clean, file = paste0(Sys.Date(),"_Seurat_object_with_scSorter_scBase_scMAGIC_Cosgrove_MCA_Raw_output.RData"))
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###

################################################################################################################################################################


################################################################################################################################################################
##### 4. Calculate cluster markers and plot potential marker gene expression


############### SCT SNN 0.6
# find markers for every cluster compared to all remaining cells, report only the positive ones
muscle.clean <- SetIdent(object = muscle.clean, value = 'SCT_snn_res.0.6')

muscle.markers_0.6 <- FindAllMarkers(muscle.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

# get top 5 and top 10
muscle.markers_0.6 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)  -> top5_0.6
muscle.markers_0.6 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10_0.6

# plot heatmap
my.heat.2_0.6 <- DoHeatmap(muscle.clean, features = top5_0.6$gene, group.by = 'SCT_snn_res.0.6',   size = 3) + scale_fill_gradientn(colors = c("blue", "white", "red"))  + theme(axis.text.y = element_text(size = 1))

png(paste0(Sys.Date(),"Top5_marker_heatmap_SCT_snn_res_0.6.png"), width = 40, height = 20, units = "cm", res = 300)
my.heat.2_0.6
dev.off()

# Dotplot of top5 findmarkers markers by cluster
pdf(paste(Sys.Date(),"Muscle_Dotplot_TOP5_markers_by_SCT_snn_res.0.6_CLEAN.pdf", sep = "_"), height = 12, width = 18)
DotPlot(muscle.clean,
        features = unique(top5_0.6$gene),
        group.by = "SCT_snn_res.0.6") + RotatedAxis()
dev.off()


# write markers to file
write.table(data.frame(top5_0.6), file = paste(Sys.Date(),"Seurat_top5_markers_SCT_snn_res.0.6_Clustering.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(data.frame(top10_0.6), file = paste(Sys.Date(),"Seurat_top10_markers_SCT_snn_res.0.6_Clustering.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(muscle.markers_0.6[muscle.markers_0.6$p_val_adj < 0.05,], file = paste(Sys.Date(),"Seurat_ALL_FDR5_markers_SCT_snn_res.0.6_Clustering.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
save(muscle.markers_0.6, file = paste(Sys.Date(),"Seurat_markers_SCT_snn_res.0.6.RData", sep = "_"))

###############################################################################################################################################################


################################################################################################################################################################
##### 5. Parse cell type predictions to assign an identity to each cell cluster

# load('2024-12-10_Seurat_object_with_ScSorter_singleCellBase_scMayoMap_Raw_output.RData')

################################################################################
# Parse results from each prediction by SNN cluster group (to apply majority and decrease noise)
# we are using the 0.6 clustering resolution to improve granularity of assignment
my.cluster.annot                 <- data.frame(matrix(0, length(unique(muscle.clean@meta.data$SCT_snn_res.0.6)), 16) )
rownames(my.cluster.annot)       <- paste0("Cluster_",0:(length(unique(muscle.clean@meta.data$SCT_snn_res.0.6))-1))
colnames(my.cluster.annot)       <- c("SCT_snn_res.0.6",
                                      "Top_ScSorter_scBase"         , "Top_ScSorter_scBase_Perc"      ,
                                      "Second_ScSorter_scBase"      , "Second_ScSorter_scBase_Perc"   ,
                                      "Top_ScSorter_scMayoMap"      , "Top_ScSorter_scMayoMap_Perc"   ,
                                      "Second_ScSorter_scMayoMap"   , "Second_ScSorter_scMayoMap_Perc",
                                      "Top_scMAGIC_Cosgrove"        , "Top_scMAGIC_Cosgrove_Perc"     ,
                                      "Second_scMAGIC_Cosgrove"     , "Second_scMAGIC_Cosgrove_Perc"  ,
                                      "Cell_Cycle_Phase"            , "Cell_Cycle_Phase_Perc"         ,
                                      "Top_5_marker_genes")

my.cluster.annot$SCT_snn_res.0.6 <- 0:(length(unique(muscle.clean@meta.data$SCT_snn_res.0.6))-1) # initialize

top5.genes <- data.frame(top5_0.6)

# function to parse predictions
get_top2_preds <- function (clust.preds, my.colname){
  
  # tabulate by cell type called by algorithm
  tab.res <- table(meta.snn.sub[,my.colname])
  
  # Sort to get top 2
  tab.res.sort <- sort(tab.res, decreasing = T)
  
  # parse information of top 2 most frequent predictions in the cluster
  a  <- names(tab.res.sort)[1]
  b  <- 100*round(tab.res.sort[1]/sum(tab.res.sort), digits = 4)
  
  c  <- names(tab.res.sort)[2]
  d  <- 100*round(tab.res.sort[2]/sum(tab.res.sort), digits = 4)
  
  return(c(a,b,c,d))
}

####
for (i in 1:nrow(my.cluster.annot)) {
  
  my.snn.clust <- my.cluster.annot$SCT_snn_res.0.6[i]
  
  # subset cluster from metadata dataframe to extract info
  meta.snn.sub <- muscle.clean@meta.data[muscle.clean@meta.data$SCT_snn_res.0.6 == my.snn.clust,]
  
  # parse information for ScSorter singleCellBase
  scSorter.parse <- get_top2_preds(meta.snn.sub, "ScSorter_scBase")
  my.cluster.annot[i,]$Top_ScSorter_scBase             <- scSorter.parse[1]
  my.cluster.annot[i,]$Top_ScSorter_scBase_Perc        <- scSorter.parse[2]
  my.cluster.annot[i,]$Second_ScSorter_scBase          <- scSorter.parse[3]
  my.cluster.annot[i,]$Second_ScSorter_scBase_Perc     <- scSorter.parse[4]
  
  # parse information for ScSorter scMayoMap
  scSorter.scMayoMap.parse <- get_top2_preds(meta.snn.sub, "ScSorter_scMayoMap")
  my.cluster.annot[i,]$Top_ScSorter_scMayoMap             <- scSorter.scMayoMap.parse[1]
  my.cluster.annot[i,]$Top_ScSorter_scMayoMap_Perc        <- scSorter.scMayoMap.parse[2]
  my.cluster.annot[i,]$Second_ScSorter_scMayoMap          <- scSorter.scMayoMap.parse[3]
  my.cluster.annot[i,]$Second_ScSorter_scMayoMap_Perc     <- scSorter.scMayoMap.parse[4]
  
  # parse information for scMAGIC_Cosgrove
  scMAGIC_Cosgrove.parse <- get_top2_preds(meta.snn.sub, "scMAGIC_Cosgrove")
  my.cluster.annot[i,]$Top_scMAGIC_Cosgrove                <- scMAGIC_Cosgrove.parse[1]
  my.cluster.annot[i,]$Top_scMAGIC_Cosgrove_Perc           <- scMAGIC_Cosgrove.parse[2]
  my.cluster.annot[i,]$Second_scMAGIC_Cosgrove             <- scMAGIC_Cosgrove.parse[3]
  my.cluster.annot[i,]$Second_scMAGIC_Cosgrove_Perc        <- scMAGIC_Cosgrove.parse[4]
  
  # cell cycle information
  cell.cyc      <- table(meta.snn.sub[,"Phase"])
  cell.cyc.perc <- 100*round(cell.cyc/sum(cell.cyc), digits = 4)
  
  my.cluster.annot[i,]$Cell_Cycle_Phase        <- names(cell.cyc.perc)[which.max(cell.cyc.perc)]
  my.cluster.annot[i,]$Cell_Cycle_Phase_Perc   <- cell.cyc.perc[which.max(cell.cyc.perc)]
  
  
  # cluster markers
  my.cluster.annot[i,]$Top_5_marker_genes <- paste(top5.genes[top5.genes$cluster == my.snn.clust,]$gene, collapse = ",")
  
  
}

write.table(my.cluster.annot, file = paste(Sys.Date(),"Parsed_Cell_Annotation_Results_by_SCT_snn_res.0.6_Clusters_with_markers_ALL.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = T)


FeaturePlot(muscle.clean, reduction = "umap", features = "Ptprc", pt.size	= 2, raster = T, raster.dpi = c(1024, 1024), label = T)
FeaturePlot(muscle.clean, reduction = "umap", features = "Pax7", pt.size	= 2, raster = T, raster.dpi = c(1024, 1024), label = T)
FeaturePlot(muscle.clean, reduction = "umap", features = "Cd3d", pt.size	= 2, raster = T, raster.dpi = c(1024, 1024), label = T)


##### add additional marker plotting to help annotate
pdf(paste(Sys.Date(),"Muscle_Dotplot_PROLIFERATION_markers_by_SCT_snn_res.0.6.pdf", sep = "_"), height = 8, width = 3.5)
DotPlot(muscle.clean,
        features = c("Pcna","Mki67"),                   # Proliferation
        group.by = "SCT_snn_res.0.6") + RotatedAxis()
dev.off()


pdf(paste(Sys.Date(),"Muscle_Dotplot_Misc_markers_by_SCT_snn_res.0.6.pdf", sep = "_"), height = 8, width = 4.5)
DotPlot(muscle.clean, features = c("Ptprc", "Mbp", "Pax7","Pdgfra","Dcn","Cd34","Col5a3", "Igfbp7"), group.by = "SCT_snn_res.0.6") + RotatedAxis()
dev.off()

################################################################################################################################################################

################################################################################################################################################################
##### 6. Assign cell type predictions to each cell cluster

# Used the predictions combined with marker gene expression to annotate in Excel

# Import predictions
my.annot <- read_xlsx("2024-12-10_Parsed_Cell_Annotation_Results_by_SCT_snn_res.0.6_Clusters_with_markers_ALL_V2.xlsx")

# get calls from scSorter calls
muscle.clean@meta.data$Cell_Identity <- NA # initialize

for (i in 0:24) {
  muscle.clean@meta.data[muscle.clean@meta.data$SCT_snn_res.0.6 == i,]$Cell_Identity <- my.annot$`Consistent Annotation`[my.annot$SCT_snn_res.0.6 == i]
}


table(muscle.clean@meta.data$Cell_Identity)
# Endothelial               FAPs             Immune          Myonuclei           Pericyte      Schwann cells Smooth muscle cell 
# 2345               3795                383              51720                532                259                535 


length(unique(muscle.clean@meta.data$Cell_Identity)) # 7

# create your own color palette based on `seedcolors`
set.seed(123123) # stabilize
P16 = createPalette(8+3,  c("#ff0000", "#00ff00", "#0000ff"))
swatch(P16)

write.table(cbind("col" = P16[-c(1:4)], "cell_type" = sort(unique(muscle.clean@meta.data$Cell_Identity))), file = paste0(Sys.Date(),"_color_palette_annotation.txt"), sep = "\t", row.names = F)

# Plot Raw output
pdf(paste(Sys.Date(),"Muscle_UMAP_color_by_deNovo_annotation_RASTER.pdf", sep = "_"), height = 3.5, width = 5)
DimPlot(muscle.clean, reduction = "umap", group.by = "Cell_Identity", raster = T,
        raster.dpi = c(850,850), cols = as.vector(P16[-c(1:4)]))
dev.off()

# Save object with raw annotations
save(muscle.clean, file = paste0(Sys.Date(),"_Seurat_object_with_Manual_annotation.RData"))
###############################################################################################################################################################


################################################################################################################################################################
##### 7. Marker plotting and QC

# Load object
load('2024-12-13_Seurat_object_with_Manual_annotation.RData')

# recalculate marker genes on the manual annotations
muscle.clean <- SetIdent(object = muscle.clean, value = 'Cell_Identity')

muscle.markers.celltype <- FindAllMarkers(muscle.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

# get top 5 and top 10
muscle.markers.celltype %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)  -> top5_annot
muscle.markers.celltype %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10_annot

# reorder
annot.muscle.markers_celltype.v2 <- muscle.markers.celltype[,c("cluster","avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj", "gene")]
annot.muscle.markers_celltype.v2 <- annot.muscle.markers_celltype.v2[order(annot.muscle.markers_celltype.v2$cluster, decreasing = F),]

# get top 5 and top 10
annot.muscle.markers_celltype.v2 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)  -> ann.top5_celltype
annot.muscle.markers_celltype.v2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> ann.top10_celltype

# write annotated markers to file
write.table( data.frame(ann.top5_celltype ), file = paste(Sys.Date(),"Seurat_top5_markers_SCT_snn_res.celltype_Clustering_ANNOT.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table( data.frame(ann.top10_celltype), file = paste(Sys.Date(),"Seurat_top10_markers_SCT_snn_res.celltype_Clustering_ANNOT.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(annot.muscle.markers_celltype.v2[annot.muscle.markers_celltype.v2$p_val_adj < 0.05,], file = paste(Sys.Date(),"Seurat_ALL_FDR5_markers_SCT_snn_res.celltype_Clustering_ANNOT.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)


############
# reorder  cell types for ease of plotting
muscle.clean@meta.data$Cell_Identity <- factor(x = muscle.clean@meta.data$Cell_Identity, 
                                               levels = rev(c("Myonuclei",
                                                              "FAPs",
                                                              "Endothelial",
                                                              "Pericyte" ,
                                                              "Immune",
                                                              "Schwann cells",
                                                              "Smooth muscle cell" )))

muscle.clean <- SetIdent(object = muscle.clean, value = 'Cell_Identity')



# ###### Plot marker genes - mix of previously known and from the marker calculation above
pdf(paste(Sys.Date(),"Muscle_Dotplot_Known_Cell_type_markers_by_ManualAnnotation.pdf", sep = "_"), height = 6, width = 10)
DotPlot(muscle.clean,
        features = c( "Acta1", "Ttn","Neb",                        # Myonuclei
                      "Pdgfra","Abca8a","Dcn",                     # FAPs
                      "Pecam1", "Flt1","Mecom",                    # Endothelial
                      "Dlc1","Rgs5", "Pdgfrb",                     # Pericyte 
                      "Ptprc", "Mrc1",                             # Immune
                      "Mbp","S100b","Plp1",                        # Schwann cells
                      "Myh11", "Acta2","Myl9"                      # Smooth muscle cell
        ), 
        cols = c("gainsboro", "firebrick2"),
        col.min = -0.5, col.max = 2.5,
        group.by = "Cell_Identity") + RotatedAxis() + scale_size_area(limits = c(0,100))
dev.off()



# #### for heatmap, needs more balance: subset 200 cells per cell type to plot (otherwise, readability is terrible)
# https://satijalab.org/seurat/v3.0/multimodal_vignette.html
small.plot <- subset(muscle.clean, downsample = 250)

# Get normalized RNA values
small.plot <- NormalizeData(small.plot)
small.plot <- FindVariableFeatures(small.plot, selection.method = "vst", nfeatures = 5000)
small.plot <- ScaleData(small.plot, features = rownames(small.plot))

png(paste(Sys.Date(),"Muscle_Heatmap_Known_Cell_type_markers_by_ManualAnnotation_250CELLS_DS.png", sep = "_"), width = 14, height = 10, res = 600, units = 'cm')
DoHeatmap(small.plot, features = c( "Acta1", "Ttn","Neb",                       # Myonuclei
                                    "Pdgfra","Abca8a","Dcn",                     # FAPs
                                    "Pecam1", "Flt1","Mecom",                    # Endothelial
                                    "Dlc1","Rgs5", "Pdgfrb",                            # Pericyte 
                                    "Ptprc", "Mrc1",                             # Immune
                                    "Mbp","S100b","Plp1",                        # Schwann cells
                                    "Myh11", "Acta2","Myl9"                      # Smooth muscle cell
), size = 3, group.by = "Cell_Identity") + scale_fill_gradientn(colors = c("darkblue", "white", "red"))
dev.off()
################################################################################################################################################################


################################################################################################################################################################
##### 8. Annotate types of Myofibers
myoN.seurat <- subset(muscle.clean, subset = Cell_Identity %in% "Myonuclei")    # Myonuclei, 51720 cells

# recalculate
myoN.seurat <- FindNeighbors(myoN.seurat, dims = 1:18)  # same number of PCs as full object
myoN.seurat <- FindClusters(myoN.seurat, resolution = c(0.01,0.1,0.25,0.6))
# Number of communities: 5
# Number of communities: 9
# Number of communities: 12
# Number of communities: 17


# https://onlinelibrary.wiley.com/doi/10.1002/jcsm.13023
# https://www.nature.com/articles/s41467-024-53510-z
# Myh4 (Type IIb myonuclei), Myh1 (Type IIx myonuclei), Myh2 (Type IIa myonuclei, and Myh7 (Type I myonuclei) 
# Type I      : Myh7 (slow)
# Type IIa    : Myh2
# Type IIb    : Myh4
# Type IIx    : Myh1


pdf(paste(Sys.Date(),"MYONUCLEI_ONLY_Muscle_Exercise_TFEB_Singlets_UMAP_snn_0.1_LABEL.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(myoN.seurat, reduction = "umap", group.by = "SCT_snn_res.0.1", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024), label = T)
dev.off()

myoN.seurat <- SetIdent(object = myoN.seurat, value = 'SCT_snn_res.0.1')

pdf(paste(Sys.Date(),"MYONUCLEI_ONLY_UMAP_TypeMarker_Muscle_Exercise_TFEB_Singlets_UMAP_snn_0.1_LABEL.pdf", sep = "_"), height = 5, width = 6.5)
FeaturePlot(myoN.seurat, reduction = "umap", features = c("Myh7",
                                                          "Myh2", 
                                                          "Myh4", 
                                                          "Myh1"), pt.size	= 2, raster = T, raster.dpi = c(1024, 1024), label = T)
dev.off()

pdf(paste(Sys.Date(),"MYONUCLEI_ONLY_RidgePlot_TypeMarker_Muscle_Exercise_TFEB_Singlets_UMAP_snn_0.1_LABEL.pdf", sep = "_"), height = 5, width = 6.5)
RidgePlot(myoN.seurat,features = c("Myh7",
                                   "Myh2", 
                                   "Myh4", 
                                   "Myh1"))
dev.off()


### Port data
myoN.seurat@meta.data$Myonuclei_Type <- NA
myoN.seurat@meta.data$Myonuclei_Type[myoN.seurat@meta.data$SCT_snn_res.0.1 == 0] <- "Myonuclei_IIb" # Myh4 high
myoN.seurat@meta.data$Myonuclei_Type[myoN.seurat@meta.data$SCT_snn_res.0.1 == 1] <- "Myonuclei_IIx" # Myh1 high
myoN.seurat@meta.data$Myonuclei_Type[myoN.seurat@meta.data$SCT_snn_res.0.1 == 2] <- "Myonuclei_IIa" # Myh2 high
myoN.seurat@meta.data$Myonuclei_Type[myoN.seurat@meta.data$SCT_snn_res.0.1 == 3] <- "Myonuclei_IIb" # Myh4 high
myoN.seurat@meta.data$Myonuclei_Type[myoN.seurat@meta.data$SCT_snn_res.0.1 == 4] <- "Myonuclei_IIb" # Myh4 high
myoN.seurat@meta.data$Myonuclei_Type[myoN.seurat@meta.data$SCT_snn_res.0.1 == 5] <- "Myonuclei_IIb" # Myh4 high
myoN.seurat@meta.data$Myonuclei_Type[myoN.seurat@meta.data$SCT_snn_res.0.1 == 6] <- "Myonuclei_IIb" # Myh4 high
myoN.seurat@meta.data$Myonuclei_Type[myoN.seurat@meta.data$SCT_snn_res.0.1 == 7] <- "Myonuclei_IIb" # Myh4 high
myoN.seurat@meta.data$Myonuclei_Type[myoN.seurat@meta.data$SCT_snn_res.0.1 == 8] <- "Myonuclei_IIb" # Myh4 high

pdf(paste(Sys.Date(),"MYONUCLEI_ONLY_UMAP_MyoN_Type_Muscle_Exercise_TFEB_Singlets_UMAP.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(myoN.seurat, reduction = "umap", group.by = "Myonuclei_Type", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024), label = T)
dev.off()

pdf(paste(Sys.Date(),"MYONUCLEI_ONLY_Dotplot_Known_Cell_type_markers_by_ManualAnnotation.pdf", sep = "_"), height = 3, width = 6)
DotPlot(myoN.seurat,
        features = c("Myh7",  # Type I      : Myh7 (slow)
                     "Myh2",  # Type IIa    : Myh2
                     "Myh4",  # Type IIb    : Myh4
                     "Myh1"   # Type IIx    : Myh1
        ), 
        cols = c("gainsboro", "firebrick2"),
        col.min = -0.5, col.max = 2.5,
        group.by = "Myonuclei_Type") + RotatedAxis() + scale_size_area(limits = c(0,100))
dev.off()

# Save MYONUCLEI object with FINAL annotations
save(myoN.seurat, file = paste0(Sys.Date(),"_Seurat_MYONUCLEI_object_with_Manual_annotation.RData"))
################################################################################################################################################################

################################################################################################################################################################
##### 9. Add myofiber type to parent object

load('2024-12-13_Seurat_object_with_Manual_annotation.RData')
load('2024-12-13_Seurat_MYONUCLEI_object_with_Manual_annotation.RData')

# populate meta-data
muscle.clean@meta.data$Cell_Identity_FINAL               <- as.character(muscle.clean@meta.data$Cell_Identity)
muscle.clean@meta.data[rownames(myoN.seurat@meta.data),]$Cell_Identity_FINAL <- as.character(myoN.seurat@meta.data$Myonuclei_Type)

muscle.clean@meta.data$Cell_Identity_FINAL               <- factor( x = muscle.clean@meta.data$Cell_Identity_FINAL,
                                                                    levels = c("Myonuclei_IIa",
                                                                               "Myonuclei_IIb",
                                                                               "Myonuclei_IIx",
                                                                               "FAPs",
                                                                               "Endothelial",
                                                                               "Pericyte" ,
                                                                               "Immune",
                                                                               "Schwann cells",
                                                                               "Smooth muscle cell" )) 

# calculate membership
table(muscle.clean@meta.data$Cell_Identity_FINAL)
# Endothelial               FAPs             Immune      Myonuclei_IIa      Myonuclei_IIb      Myonuclei_IIx           Pericyte 
# 2345               3795                383               3075              37094              11551                532 
# Schwann cells Smooth muscle cell 
# 259                535 


# create your own color palette based on `seedcolors`
set.seed(987654321) # stabilize
P16 = createPalette(9+3,  c("#ff0000", "#00ff00", "#0000ff","#800080","#ff9900"))
swatch(P16)

write.table(cbind("col" = P16[-c(1:3)], "cell_type" = sort(unique(muscle.clean@meta.data$Cell_Identity_FINAL))), file = paste0(Sys.Date(),"_color_palette_annotation_FINAL.txt"), sep = "\t", row.names = F)

# Plot Raw output
pdf(paste(Sys.Date(),"Muscle_UMAP_color_by_deNovo_annotation_IncludingFiberType_RASTER.pdf", sep = "_"), height = 8, width = 9.5)
DimPlot(muscle.clean, reduction = "umap", group.by = "Cell_Identity_FINAL", raster = T,
        raster.dpi = c(1024,1024), 
        cols = as.vector(P16[-c(1:3)])
)
dev.off()

# Save object with FINAL annotations
save(muscle.clean, file = paste0(Sys.Date(),"_Seurat_object_with_Manual_annotation_FINAL.RData"))


muscle.clean@meta.data$Cell_Identity_FINAL               <- factor( x = muscle.clean@meta.data$Cell_Identity_FINAL,
                                                                    levels = rev(c("Myonuclei_IIa",
                                                                               "Myonuclei_IIb",
                                                                               "Myonuclei_IIx",
                                                                               "FAPs",
                                                                               "Endothelial",
                                                                               "Pericyte" ,
                                                                               "Immune",
                                                                               "Schwann cells",
                                                                               "Smooth muscle cell" )))

# ###### Plot marker genes - mix of previously known and from the marker calculation above
pdf(paste(Sys.Date(),"Muscle_Dotplot_Known_Cell_type_markers_by_ManualAnnotation_with_MyofiberType.pdf", sep = "_"), height = 7, width = 10)
DotPlot(muscle.clean,
        features = c( "Acta1", "Ttn","Neb",                       # Myonuclei
                      "Myh2","Myh4","Myh1",
                      "Pdgfra","Abca8a","Dcn",                     # FAPs
                      "Pecam1", "Flt1","Mecom",                    # Endothelial
                      "Dlc1","Rgs5", "Pdgfrb",                            # Pericyte 
                      "Ptprc", "Mrc1",                             # Immune
                      "Mbp","S100b","Plp1",                        # Schwann cells
                      "Myh11", "Acta2","Myl9"                      # Smooth muscle cell
        ), 
        cols = c("gainsboro", "firebrick2"),
        col.min = -0.5, col.max = 2.5,
        group.by = "Cell_Identity_FINAL") + RotatedAxis() + scale_size_area(limits = c(0,100))
dev.off()
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###


###############################################################################################################################################################
sink(file = paste(Sys.Date(),"_Muscle_ANNOT_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()



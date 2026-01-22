setwd('/Volumes/BB_HQ_4/Collaboration/Cortes_lab/Muscle_snRNAseq/Seurat_Processing')
options(stringsAsFactors = F)
options (future.globals.maxSize = 32000 * 1024^2)

# Load packages
library('Seurat')    # 
library(sctransform) # 
library(clustree)    # 
library(scales)      # 
library(harmony)     #


##################################################################################
# 2024-11-29
# Clean dataset before downstream processing
##################################################################################


#####################################################################################################################
#### 1. Load Cleaned up Seurat Object

# load up cleaned up singlets
load("../Preprocessing/2024-11-29_Muscle_data_Cortes_lab_Seurat_object_SINGLETS_ONLY.RData")

muscle.singlets
# An object of class Seurat 
# 57541 features across 59569 samples within 2 assays 
# Active assay: SCT (23845 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

# bring RNA as main assay again
DefaultAssay(muscle.singlets) <- "RNA"

table(muscle.singlets@meta.data$Group)
# F_runner F_sedentary      F_TFEB    M_runner M_sedentary      M_TFEB 
#    10189        9678        8733       12399       11367        7203 

# will need to clean up and rerun SCT 
muscle.singlets[['SCT']] <- NULL
muscle.singlets


#### Filter genes with expression that is too sparse
min.value = 0
min.cells = 250
genes.use <- rownames(muscle.singlets@assays$RNA)
num.cells <- Matrix::rowSums(muscle.singlets@assays$RNA@counts > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])

# remove low/null genes
muscle.singlets <- subset(muscle.singlets, features = genes.use)
muscle.singlets
# An object of class Seurat 
# 13185 features across 59569 samples within 1 assay 
# Active assay: RNA (13185 features, 0 variable features)
# 2 layers present: counts, data
# 2 dimensional reductions calculated: pca, umap

# remove irrelevant columns (previous SCT on unmerged objects)
sct.cols <- grep("SCT",colnames(muscle.singlets@meta.data))

muscle.singlets@meta.data <- muscle.singlets@meta.data[,-c(sct.cols)]

save(muscle.singlets, file = paste0(Sys.Date(),"_Muscle_Cortes_Lab_Seurat_object_preNorm.RData"))
################################################################################################################################################################


################################################################################################################################################################
#### 2. Scaling the data and removing unwanted sources of variation

# Run SCT 
muscle.clean <- SCTransform(object = muscle.singlets, vars.to.regress =  c("nFeature_RNA", "nCount_RNA", "percent.mito", "Phase"), variable.features.n = 5000)

# Identify the 20 most highly variable genes
muscle.var.top20 <- head(VariableFeatures(muscle.clean), 20)
# [1] "Meg3"    "Tnnc2"   "mt-Co3"  "Rian"    "mt-Atp6" "mt-Co1"  "Kcnh7"   "Gpc6"    "Ebf1"    "mt-Co2"  "Csmd1"   "mt-Nd4"  "Dlc1"    "F13a1"   "Myh2"    "Ano4"   
# [17] "mt-Cytb" "C7"      "Erbb4"   "Tenm2"

muscle.clean
# An object of class Seurat 
# 26370 features across 59569 samples within 2 assays 
# Active assay: SCT (13185 features, 5000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

table(muscle.clean@meta.data$Treat.gp, muscle.clean@meta.data$Sex)
#          F     M
# RUN  10189 12399
# SED   9678 11367
# TFEB  8733  7203

table(muscle.clean@meta.data$Group)
# F_runner F_sedentary      F_TFEB    M_runner M_sedentary      M_TFEB 
# 10189        9678        8733       12399       11367        7203 

save(muscle.clean, file = paste0(Sys.Date(),"_Muscle_Exercise_TFEB_Seurat_object_SCT.RData"))

# clean up unnecessary objects from memory
rm(muscle.singlets)
###############################################################################################################################################################


###############################################################################################################################################################
##### 3. Run dimensionality reduction

# Run dimensionality reduction PCA
muscle.clean <- RunPCA(muscle.clean, npcs = 40)

# Determine the ‘dimensionality’ of the dataset
# Approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
pdf(paste0(Sys.Date(), "_Muscle_Exercise_TFEB_elbowplot.pdf"), height = 5, width= 6)
ElbowPlot(muscle.clean, ndims = 40)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- muscle.clean[["pca"]]@stdev / sum(muscle.clean[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 33

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 18

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 18

# Based on these metrics, first 14 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_Muscle_Exercise_TFEB_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()
###############################################################################


# Calculate UMAP
muscle.clean <- RunUMAP(muscle.clean, dims = 1:pcs)
muscle.clean
# An object of class Seurat 
# 26370 features across 59569 samples within 2 assays 
# Active assay: SCT (13185 features, 5000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

save(muscle.clean, file = paste0(Sys.Date(),"_Muscle_Exercise_TFEB_Seurat_object_SCT_with_UMAP.RData"))

################  Plot summary UMAPs  ################

# QC UMAPs
pdf(paste(Sys.Date(),"Muscle_Exercise_TFEB_Singlets_UMAP_color_by_Sex.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(muscle.clean, reduction = "umap", group.by = "Sex", pt.size	= 2, shuffle = T,
        cols = c(alpha("deeppink", alpha = 0.5 ), alpha("deepskyblue", alpha = 0.5 ) ) , raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Muscle_Exercise_TFEB_Singlets_UMAP_color_by_Treatment.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(muscle.clean, reduction = "umap", group.by = "Treat.gp", pt.size	= 2, shuffle = T,
        cols = c(alpha("goldenrod1", alpha = 0.5 ), alpha("firebrick3", alpha = 0.5 ), alpha("steelblue2", alpha = 0.5 )  ) , raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Muscle_Exercise_TFEB_Singlets_UMAP_color_by_Group.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(muscle.clean, reduction = "umap", group.by = "Group", pt.size	= 2, shuffle = T,
        cols = c(alpha("lightpink2", alpha = 0.5    ), alpha("deeppink", alpha = 0.5    ), alpha("mediumpurple4", alpha = 0.5    ) ,
                 alpha("mediumturquoise", alpha = 0.5 ), alpha("deepskyblue", alpha = 0.5 ), alpha("royalblue3", alpha = 0.5 ) ) ,
        raster = T, raster.dpi = c(1024, 1024))
dev.off()

###############################################################################################################################################################


###############################################################################################################################################################
sink(file = paste(Sys.Date(),"_Muscle_Exercise_TFEB_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()


setwd('/Volumes/BB_HQ_4/Collaboration/Cortes_lab/Muscle_snRNAseq/Seurat_Processing')
options(stringsAsFactors = F)
options (future.globals.maxSize = 32000 * 1024^2)

# Load packages
library('Seurat')         # 
library(sctransform)      # 
library("singleCellTK")   # 
library("scater")         # 


#####################################################################################################################
# 2024-12-13
# Plot variance explained
#
#####################################################################################################################

#####################################################################################################################
#### 1. Load Cleaned up annotated Seurat Object

# Import final annotation
load('../Seurat_Processing/2024-12-13_Seurat_object_with_Manual_annotation_FINAL.RData')

muscle.clean
# An object of class Seurat 
# 26370 features across 59569 samples within 2 assays 
# Active assay: SCT (13185 features, 5000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

View(muscle.clean@meta.data)

muscle.clean.sce <- as.SingleCellExperiment(muscle.clean)

# Computing variance explained on the log-counts
# https://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html
vars <- getVarianceExplained(muscle.clean.sce,
                             variables=c("Sex", "Treat.gp", "Group", "Phase", "Cell_Identity_FINAL"))
head(vars)
#                 Sex    Treat.gp      Group       Phase Cell_Identity_FINAL
# Xkr4   0.0169601546 0.000647253 0.02075328 0.018848566          7.87922361
# Rp1    0.0002441945 0.017935679 0.02630331 0.008087985          4.14806147
# Sox17  0.0105624179 0.006911576 0.08894705 0.021420517         10.99758648
# Mrpl15 0.0205797262 0.132918104 0.36586348 0.015708451          0.36855629
# Lypla1 0.0953325453 0.085810982 0.34893097 0.053778494          0.04728759
# Tcea1  0.0103906237 0.002585575 0.15691257 0.036530826          0.09232339

pdf(paste0(Sys.Date(),"_Variance_Explained_plot_Muscle_Dataset.pdf"), width = 5, height = 4)
plotExplanatoryVariables(vars[,c("Sex", "Treat.gp", "Cell_Identity_FINAL")])
dev.off()

pdf(paste0(Sys.Date(),"_Variance_Explained_plot_Muscle_Dataset_VIOLINS.pdf"), width = 5, height = 4)
vioplot::vioplot(vars[,c("Sex", "Treat.gp", "Cell_Identity_FINAL")], las = 2, col = c("orange","purple","darkblue"), ylim = c(0,100), ylab = "% variance explained")
dev.off()

###############################################################################################################################################################
sink(file = paste(Sys.Date(),"_Muscle_VarExp_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()

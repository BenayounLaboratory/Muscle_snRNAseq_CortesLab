setwd('/Volumes/BB_HQ_4/Collaboration/Cortes_lab/Muscle_snRNAseq/AUGUR/')
options(stringsAsFactors=FALSE)

# scRNA tools
library(Seurat)
library(Augur)
library(viridis)

# 2024-12-13
# run AUGUR

######################################################################
# load annotated Seurat Objects
load('../Seurat_Processing/2024-12-13_Seurat_object_with_Manual_annotation_FINAL.RData')

# Split Seurat Objects
muscle.clean.f <- subset(muscle.clean, subset = Sex %in% "F")    # 28600 cells
muscle.clean.m <- subset(muscle.clean, subset = Sex %in% "M")    # 30969 cells

muscle.clean.f
# An object of class Seurat 
# 26370 features across 28600 samples within 2 assays 
# Active assay: SCT (13185 features, 5000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

muscle.clean.m
# An object of class Seurat 
# 26370 features across 30969 samples within 2 assays 
# Active assay: SCT (13185 features, 5000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

#######################

augur.muscle.all <-  calculate_auc(muscle.clean,
                                   cell_type_col =   "Cell_Identity_FINAL",
                                   label_col     =   "Group",
                                   n_threads     =   1,
                                   min_cells     =   100)

augur.muscle.f   <-  calculate_auc(muscle.clean.f,
                                   cell_type_col =   "Cell_Identity_FINAL",
                                   label_col     =   "Group",
                                   n_threads     =   1,
                                   min_cells     =   100)

augur.muscle.m   <-  calculate_auc(muscle.clean.m,
                                   cell_type_col =   "Cell_Identity_FINAL",
                                   label_col     =   "Group",
                                   n_threads     =   1,
                                   min_cells     =   100)


pdf(paste0(Sys.Date(),"_Augur_Muscle_scRNAseq_ALL.pdf"), width = 3, height = 3)
plot_umap(augur.muscle.all, muscle.clean, cell_type_col = "Cell_Identity_FINAL", palette = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","lightgrey","#CCCCFF","#9999FF","#333399")))(50))
dev.off()

pdf(paste0(Sys.Date(),"_Augur_Muscle_scRNAseq_Female.pdf"), width = 3, height = 3)
plot_umap(augur.muscle.f, muscle.clean.f, cell_type_col = "Cell_Identity_FINAL", palette = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","lightgrey","#CCCCFF","#9999FF","#333399")))(50))
dev.off()

pdf(paste0(Sys.Date(),"_Augur_Muscle_scRNAseq_Male.pdf"), width = 3, height = 3)
plot_umap(augur.muscle.m, muscle.clean.m, cell_type_col = "Cell_Identity_FINAL", palette = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","lightgrey","#CCCCFF","#9999FF","#333399")))(50))
dev.off()

pdf(paste0(Sys.Date(),"_Augur_Muscle_scRNAseq_ALL_Lollipop.pdf"), width = 3, height = 3)
plot_lollipop(augur.muscle.all )
dev.off()

pdf(paste0(Sys.Date(),"_Augur_Muscle_scRNAseq_Female_Lollipop.pdf"), width = 3, height = 3)
plot_lollipop(augur.muscle.f )
dev.off()

pdf(paste0(Sys.Date(),"_Augur_Muscle_scRNAseq_Male_Lollipop.pdf"), width = 3, height = 3)
plot_lollipop(augur.muscle.m)
dev.off()
###################################################################################################################


#######################
sink(file = paste(Sys.Date(),"_Muscle_scRNAseq_AUGUR_session_Info.txt", sep =""))
sessionInfo()
sink()

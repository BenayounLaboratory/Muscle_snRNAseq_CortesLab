setwd('/Volumes/BB_HQ_4/Collaboration/Cortes_lab/Muscle_snRNAseq/scProportionTest')
options(stringsAsFactors = F)

library(Seurat)
library(scProportionTest)
library(ggplot2)

# 2024-12-13
# Run scProportionTest on each sex

# 2025-03-13
# replot with extended x axis

######################################################################
# load annotated Seurat Objects
load('../Seurat_Processing/2024-12-13_Seurat_object_with_Manual_annotation_FINAL.RData')

muscle.clean
# An object of class Seurat 
# 26370 features across 59569 samples within 2 assays 
# Active assay: SCT (13185 features, 5000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

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

# update to include sample ID
table(muscle.clean.f$Sample, muscle.clean.f$Cell_Identity_FINAL)
#         Mstnuclei_IIa Mstnuclei_IIb Mstnuclei_IIx FAPs Endothelial Pericyte Immune Schwann cells Smooth muscle cell
# F_RUN_1           400          3244          1052  285         223       17     34            24                 49
# F_RUN_2           246          3144           967  183         171       35     41            16                 58
# F_SED_1           187          3579           943  306         212       34     44            22                 34
# F_SED_2           328          2639           927  223         103       15     17            15                 50
# F_TFB_1           227          2873           901  306         160       60     29            15                 42
# F_TFB_2           229          2464          1018  146         127       44     16            18                 58

table(muscle.clean.m$Sample, muscle.clean.m$Cell_Identity_FINAL)
#         Mstnuclei_IIa Mstnuclei_IIb Mstnuclei_IIx FAPs Endothelial Pericyte Immune Schwann cells Smooth muscle cell
# M_RUN_1           199          2931           822  332         174       23     33            18                 32
# M_RUN_2           384          4965          1298  656         320       98     30            23                 61
# M_SED_1           238          4317          1012  480         257      105     50            22                 61
# M_SED_2           184          3358           713  252         222       29     26            20                 21
# M_TFB_1           214          1638           868  298         171       31     38            39                 38
# M_TFB_2           239          1942          1030  328         205       41     25            27                 31


###################################################################
# create prop test objects
muscle.clean.f.prop_test <- sc_utils(muscle.clean.f)
muscle.clean.m.prop_test <- sc_utils(muscle.clean.m)

# Once the object is created, the permutation testing and bootstrapping can be run.
set.seed(1234567890)

### Females
f.prop_test.ex <- permutation_test(muscle.clean.f.prop_test,
                                   cluster_identity = "Cell_Identity_FINAL",
                                   sample_1         = "F_sedentary",
                                   sample_2         = "F_runner",
                                   sample_identity  = "Group")

f.prop_test.tfeb <- permutation_test(muscle.clean.f.prop_test,
                                     cluster_identity = "Cell_Identity_FINAL",
                                     sample_1         = "F_sedentary",
                                     sample_2         = "F_TFEB",
                                     sample_identity  = "Group")

### Males
m.prop_test.ex <- permutation_test(muscle.clean.m.prop_test,
                                   cluster_identity  = "Cell_Identity_FINAL",
                                   sample_1          = "M_sedentary",
                                   sample_2          = "M_runner",
                                   sample_identity   = "Group")

m.prop_test.tfeb <- permutation_test(muscle.clean.m.prop_test,
                                     cluster_identity = "Cell_Identity_FINAL",
                                     sample_1         = "M_sedentary",
                                     sample_2         = "M_TFEB",
                                     sample_identity  = "Group")


### Parse and Plot

#############################################
## Modify function
permutation_plot_mod <- function (se.data, st.data, FDR_threshold = 0.05, cols_vals = c("dodgerblue"   , "firebrick","grey"), my.title) {
  
  # parse data
  plot_data.se <-se.data@results$permutation
  plot_data.st <-st.data@results$permutation
  
  # add category
  plot_data.se$comp <- "Exercise"
  plot_data.st$comp <- "TFEb"
  
  # merge data from comparisons
  colnames(plot_data.se) <- colnames(plot_data.st)
  plot_data <- rbind(plot_data.se,
                     plot_data.st)
  
  # fix order alphabetically
  plot_data$clusters <- factor(plot_data$clusters, levels = rev(plot_data.se$clusters))
  
  # fix plotting order
  plot_data$comp <- factor(plot_data$comp, levels = c("Exercise","TFEb"))
  
  # get significance and color scale
  plot_data$significance <- "n.s."
  plot_data$significance[plot_data$FDR < FDR_threshold & plot_data$boot_mean_log2FD > 0] <- paste("FDR <", FDR_threshold, "(Increased)")
  plot_data$significance[plot_data$FDR < FDR_threshold & plot_data$boot_mean_log2FD < 0] <- paste("FDR <", FDR_threshold, "(Decreased)")
  plot_data$significance <- factor(plot_data$significance, levels = c(paste("FDR <", FDR_threshold, "(Decreased)"), paste("FDR <", FDR_threshold, "(Increased)"),  "n.s."))
  
  # fix range for symmetry
  max_range <- 2
  
  # plot
  p <- ggplot(plot_data, aes(x = clusters, y = obs_log2FD)) +   theme_bw() +
    geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance, shape = comp),position = position_dodge(width = -1)) +
    geom_hline(yintercept = 0) + scale_color_manual(values = cols_vals) + ylim(c(-max_range,max_range)) +
    coord_flip()  + ggtitle(my.title)
  p
  return(p)
}

#############################################

# A point-range plot of the results can then be created.
perm.f     <- permutation_plot_mod(f.prop_test.ex, f.prop_test.tfeb, my.title = "Female Exercise/TFEb")
perm.m     <- permutation_plot_mod(m.prop_test.ex, m.prop_test.tfeb, my.title = "Male Exercise/TFEb")

pdf(paste(Sys.Date(),"scProportionTest_BrainAging.pdf",sep = "_"), height = 6.5, width = 12)
gridExtra::grid.arrange(perm.f, perm.m, nrow = 1)
dev.off()
# ############################################################################################################

#######################
sink(file = paste(Sys.Date(),"_scProportionTest_Muscle_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()

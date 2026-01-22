setwd('/Volumes/BB_HQ_4/Collaboration/Cortes_lab/Muscle_snRNAseq/DEG_analysis/')
options(stringsAsFactors = F)

#### Packages
library('Seurat')         # 
library(sctransform)      # 
library("singleCellTK")   # 

library('muscat')         # 
library('DESeq2')         # 
library('sva')            # 
library('limma')          # 

library(ggplot2)          # 
library(scales)           # 
library("bitops")         # 
library(Vennerable)       # 
library(data.table)       #

library(ComplexHeatmap)   #
library(circlize)         #


theme_set(theme_bw())   

# 2024-12-13
# Process scRNAseq muscle cohorts for differential gene analysis

# 2025-03-04
# get some summary plots

###############################################################################################
# 0. preprocess Seurat object for use with muscat

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


#### Make covariate table
Female.meta <- unique(muscle.clean.f@meta.data[,c("Sample", "Treat.gp")])
Male.meta   <- unique(muscle.clean.m@meta.data[,c("Sample", "Treat.gp")])

write.table(Female.meta, file = paste0(Sys.Date(),"_Female_sample_metadata_table.txt"), sep = "\t", row.names = F, quote = F)
write.table(Male.meta, file = paste0(Sys.Date(),"_Male_sample_metadata_table.txt"), sep = "\t", row.names = F, quote = F)

# bring RNA as main assay for processing
DefaultAssay(muscle.clean.f) <- "RNA"
DefaultAssay(muscle.clean.m) <- "RNA"

# convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
muscle.clean.Female.sce <- as.SingleCellExperiment(muscle.clean.f)
muscle.clean.Male.sce   <- as.SingleCellExperiment(muscle.clean.m)
save(muscle.clean.Female.sce, muscle.clean.Male.sce, 
     file = paste(Sys.Date(),"SingleCellExperimnents_objects_PER_SEX.RData",sep = "_"))

rm(muscle.clean) # free up some memory
rm(muscle.clean.f, muscle.clean.m) # free up some memory
###############################################################################################


###############################################################################################
# 1. Run muscat for pseudobulking and extraction of samples

# Clean workspace and reload necessary objects
load('2024-12-13_SingleCellExperimnents_objects_PER_SEX.RData')

##################################################
####### Data preparation   ++++   Female   #######
##################################################

muscle.Female.sce.cl <- prepSCE(muscle.clean.Female.sce, 
                                kid    = "Cell_Identity_FINAL",  # population assignments
                                gid    = "Treat.gp"           ,  # group IDs (ctrl/stim)
                                sid    = "Sample"             ,  # sample IDs (ctrl/stim.1234)
                                drop   = TRUE                 )  # drop all other colData columns

# store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk  <- length(kids <- levels(muscle.Female.sce.cl$cluster_id))
ns  <- length(sids <- levels(muscle.Female.sce.cl$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(muscle.Female.sce.cl$cluster_id, muscle.Female.sce.cl$sample_id))

# Aggregation of single-cell to pseudobulk data
pb.f <- aggregateData(muscle.Female.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# one list item per cell type
assayNames(pb.f)
# [1] "Myonuclei_IIa"      "Myonuclei_IIb"      "Myonuclei_IIx"      "FAPs"               "Endothelial"        "Pericyte"           "Immune"            
# [8] "Schwann cells"      "Smooth muscle cell"

# Number of cells in each sample and cell type
cell.per.samp.tab.f <- t(table(muscle.Female.sce.cl$cluster_id, muscle.Female.sce.cl$sample_id))

# extract pseudobulk information
counts.pb.tmp.f <- pb.f@assays@data

# get the genes with no reads in at least half the samples out, they mess up the algorithm
for (i in 1:length(counts.pb.tmp.f)) {
  my.good <- which(apply(counts.pb.tmp.f[[i]]>0, 1, sum) >= nrow(cell.per.samp.tab.f)/2) # see deseq2 vignette, need to remove too low genes
  counts.pb.tmp.f[[i]] <- counts.pb.tmp.f[[i]][my.good,]
}


###############################################
####### Data preparation   ++++   Male  #######
###############################################

muscle.Male.sce.cl <- prepSCE(muscle.clean.Male.sce, 
                              kid    = "Cell_Identity_FINAL",  # population assignments
                              gid    = "Treat.gp"           ,  # group IDs (ctrl/stim)
                              sid    = "Sample"             ,  # sample IDs (ctrl/stim.1234)
                              drop   = TRUE                 )  # drop all other colData columns

# store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk  <- length(kids <- levels(muscle.Male.sce.cl$cluster_id))
ns  <- length(sids <- levels(muscle.Male.sce.cl$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(muscle.Male.sce.cl$cluster_id, muscle.Male.sce.cl$sample_id))

# Aggregation of single-cell to pseudobulk data
pb.m <- aggregateData(muscle.Male.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# one list item per cell type
assayNames(pb.m)
# [1] "Myonuclei_IIa"      "Myonuclei_IIb"      "Myonuclei_IIx"      "FAPs"               "Endothelial"        "Pericyte"           "Immune"            
# [8] "Schwann cells"      "Smooth muscle cell"

# Number of cells in each sample and cell type
cell.per.samp.tab.m <- t(table(muscle.Male.sce.cl$cluster_id, muscle.Male.sce.cl$sample_id))

# extract pseudobulk information for samples that pass the cell number cutoff
counts.pb.tmp.m <- pb.m@assays@data

# get the genes with no reads in at least half the samples out, they mess up the algorithm
for (i in 1:length(counts.pb.tmp.m)) {
  my.good <- which(apply(counts.pb.tmp.m[[i]]>0, 1, sum) >= nrow(cell.per.samp.tab.m)/2) # see deseq2 vignette, need to remove too low genes
  counts.pb.tmp.m[[i]] <- counts.pb.tmp.m[[i]][my.good,]
}
########################################################


####################################################################
#######       Data preparation   ++++   QC Cell types        #######
####################################################################

# cell types with at least 20 cells in all samples
Female.celltype.qc <- colnames(cell.per.samp.tab.f)[colSums(cell.per.samp.tab.f  >= 20) == nrow(cell.per.samp.tab.f)]
Female.celltype.qc
# [1] "Myonuclei_IIa"      "Myonuclei_IIb"      "Myonuclei_IIx"      "FAPs"               "Endothelial"        "Smooth muscle cell"

# cell types with at least 20 cells in all samples
Male.celltype.qc <- colnames(cell.per.samp.tab.m)[colSums(cell.per.samp.tab.m  >= 20) == nrow(cell.per.samp.tab.m)]
Male.celltype.qc
#[1] "Myonuclei_IIa"      "Myonuclei_IIb"      "Myonuclei_IIx"      "FAPs"               "Endothelial"        "Pericyte"           "Immune"            
#[8] "Smooth muscle cell"

#### study cell types passing QC in both strains
both.celltype.qc   <- sort(intersect(Female.celltype.qc, Male.celltype.qc))

# extract pseudobulk information for samples that pass the cell number cutoff
counts.pb.f <- counts.pb.tmp.f[both.celltype.qc]

# extract pseudobulk information for samples that pass the cell number cutoff
counts.pb.m <- counts.pb.tmp.m[both.celltype.qc]

#### save counts
save(counts.pb.f, counts.pb.m, 
     Female.celltype.qc, Male.celltype.qc,
     file = paste0(Sys.Date(),"_muscat_PB_Female_Male_objects_QC_Clean.RData"))
###############################################################################################


##############################################################################################
# 2. Use SVA to clean up batch effects on expression and DEseq2 for DE analysis

# clean up memory and reload only muscat PBs
load('2024-12-13_muscat_PB_Female_Male_objects_QC_Clean.RData')

# run for the cell types with at least 20 cells from every each sample

##################################################
#######   DEG analysis   ++++   Female     #######
##################################################

# import metadata and order it
my.Female.meta           <- read.csv("2024-12-13_Female_sample_metadata_table.txt", sep = "\t")
my.Female.meta           <- setorder(my.Female.meta, Treat.gp, Sample)
my.Female.meta           <- my.Female.meta[c(3:4,1:2,5:6),]
rownames(my.Female.meta) <- my.Female.meta$Sample

# reorder count tables in sensical order
for  (i in 1:length(counts.pb.f)) {
  counts.pb.f[[i]] <- counts.pb.f[[i]][,my.Female.meta$Sample]
}

# Create list object to receive VST normalized counts
vst.cts.f        <- vector(mode = "list", length = length(counts.pb.f))
names(vst.cts.f) <- names(counts.pb.f)

# Create list object to receive DESeq2 results
deseq.res.list.f.run        <- vector(mode = "list", length = length(counts.pb.f))
names(deseq.res.list.f.run) <- names(counts.pb.f)
deseq.res.list.f.tfeb        <- vector(mode = "list", length = length(counts.pb.f))
names(deseq.res.list.f.tfeb) <- names(counts.pb.f)


# loop over pseudobulk data
for  (i in 1:length(counts.pb.f)) {
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Pseudobulk_Female_",names(counts.pb.f)[[i]])
  
  ###################################
  #######       Run DE      #######
  
  # legend
  my.cols  <- rep("",nrow(my.Female.meta))
  my.cols[my.Female.meta$Treat.gp %in% "SED"]  <- "deeppink"
  my.cols[my.Female.meta$Treat.gp %in% "RUN"]  <- "lightpink2"
  my.cols[my.Female.meta$Treat.gp %in% "TFEB"] <- "mediumpurple4"
  
  # get matrix using age as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = counts.pb.f[[i]],
                                colData   = my.Female.meta,
                                design    = ~ Treat.gp)
  
  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds)
  
  # plot dispersion
  my.disp.out <- paste(my.outprefix,"_dispersion_plot.pdf")
  
  pdf(my.disp.out)
  plotDispEsts(dds.deseq)
  dev.off()
  
  # get DESeq2 normalized expression value
  vst.cts.f[[i]] <- getVarianceStabilizedData(dds.deseq)
  
  # MDS analysis
  mds.result <- cmdscale(1-cor(vst.cts.f[[i]],method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  x <- mds.result[, 1]
  y <- mds.result[, 2]
  
  pdf(paste0(my.outprefix,"_MDS_plot.pdf"))
  plot(x, y,
       xlab = "MDS dimension 1", ylab = "MDS dimension 2",
       main= paste0(names(counts.pb.f)[[i]]," (MDS)"),
       cex=3, col = my.cols, pch = 16,
       cex.lab = 1.25,
       cex.axis = 1.25, las = 1)
  legend("bottomright", c("Sedentary","Runner","TFEb"),col = c("deeppink","lightpink2","mediumpurple4"), pch = 16, bty = 'n')
  dev.off()
  
  # extract gene significance by DEseq2
  res.run  <- results(dds.deseq, contrast = c("Treat.gp", "RUN" , "SED")) 
  res.tfeb <- results(dds.deseq, contrast = c("Treat.gp", "TFEB", "SED")) 
  
  # exclude genes with NA FDR value
  res.run  <- res.run [!is.na(res.run $padj),]
  res.tfeb <- res.tfeb[!is.na(res.tfeb$padj),]
  
  # store results
  deseq.res.list.f.run[[i]]        <- data.frame(res.run  )
  deseq.res.list.f.tfeb[[i]]       <- data.frame(res.tfeb )
  
  ### get sex dimorphic changes at FDR5
  genes.run  <- rownames(res.run)[res.run$padj < 0.05]
  genes.tfeb <- rownames(res.tfeb)[res.tfeb$padj < 0.05]
  
  if (length(genes.run) > 2) {
    # heatmap drawing - only if there is at least 2 gene
    my.heatmap.out <- paste0(my.outprefix,"_RUNNING_Heatmap_FDR5_GENES.pdf")
    
    pdf(my.heatmap.out, onefile = F, height = 10, width = 10)
    my.heatmap.title <- paste0(names(counts.pb.f)[[i]], " running significant (FDR<5%), ", length(genes.run), " genes")
    pheatmap::pheatmap(vst.cts.f[[i]][genes.run,],
                       cluster_cols = F,
                       cluster_rows = T,
                       colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                       show_rownames = F, scale="row",
                       main = my.heatmap.title,
                       cellwidth = 15,
                       border    = NA,
                       cellheight = 0.5)
    dev.off()
  }
  
  if (length(genes.tfeb) > 2) {
    # heatmap drawing - only if there is at least 2 gene
    my.heatmap.out <- paste0(my.outprefix,"_TFEB_Heatmap_FDR5_GENES.pdf")
    
    pdf(my.heatmap.out, onefile = F, height = 10, width = 10)
    my.heatmap.title <- paste0(names(counts.pb.f)[[i]], " TFEB significant (FDR<5%), ", length(genes.tfeb), " genes")
    pheatmap::pheatmap(vst.cts.f[[i]][genes.tfeb,],
                       cluster_cols = F,
                       cluster_rows = T,
                       colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                       show_rownames = F, scale="row",
                       main = my.heatmap.title,
                       cellwidth = 15,
                       border    = NA,
                       cellheight = 0.5)
    dev.off()
  }
  
  
  ## compare
  my.merged <- merge(data.frame(res.run), data.frame(res.tfeb), by = "row.names", suffixes = c(".run",".tfeb"))
  
  my.test <- cor.test(my.merged$log2FoldChange.run, my.merged$log2FoldChange.tfeb, method = "spearman")
  my.rho  <- signif(my.test$estimate,3)
  
  pdf(paste0(Sys.Date(),"_", names(counts.pb.f)[[i]], "_running_vs_TFEB_FC_scatterplot_FEMALE_FDR5.pdf"), height = 10, width = 10)
  smoothScatter(my.merged$log2FoldChange.run, my.merged$log2FoldChange.tfeb,
                xlim = c(-8, 8), ylim = c(-8, 8),
                xlab = "Log2FC (Running/Sedentary)",
                ylab = "Log2FC (TFEB/Sedentary)",
                main = paste0("FEMALE ", names(counts.pb.f)[[i]]) )
  abline(0, 1,  col = "red", lty = "dashed")
  abline(h = 0, col = "grey", lty = "dashed")
  abline(v = 0, col = "grey", lty = "dashed")
  text(-7, 7.5, paste("Rho = ", my.rho), pos = 4)
  text(-7, 6.5, paste("p = ",signif(my.test$p.value,2)), pos = 4)
  dev.off()
  
}

# save R object with all DEseq2 results
my.rdata <- paste0(Sys.Date(),"_pseudobulk_Muscle_cell_types_RUN_TFEB_Female_DEseq2_objects.RData")
save(deseq.res.list.f.run, deseq.res.list.f.tfeb, file = my.rdata)

my.vst <- paste0(Sys.Date(),"_pseudobulk_Muscle_cell_types_RUN_TFEB_Female_VST_data_objects.RData")
save(vst.cts.f, file = my.vst)
###############################################



###############################################
#######   DEG analysis   ++++   Male     #######
###############################################

# import metadata and order it
my.Male.meta           <- read.csv("2024-12-13_Male_sample_metadata_table.txt", sep = "\t")
my.Male.meta           <- setorder(my.Male.meta, Treat.gp, Sample)
my.Male.meta           <- my.Male.meta[c(3:4,1:2,5:6),]
rownames(my.Male.meta) <- my.Male.meta$Sample

# reorder count tables in sensical order
for  (i in 1:length(counts.pb.m)) {
  counts.pb.m[[i]] <- counts.pb.m[[i]][,my.Male.meta$Sample]
}

# Create list object to receive VST normalized counts
vst.cts.m        <- vector(mode = "list", length = length(counts.pb.m))
names(vst.cts.m) <- names(counts.pb.m)

# Create list object to receive DESeq2 results
deseq.res.list.m.run        <- vector(mode = "list", length = length(counts.pb.m))
names(deseq.res.list.m.run) <- names(counts.pb.m)
deseq.res.list.m.tfeb        <- vector(mode = "list", length = length(counts.pb.m))
names(deseq.res.list.m.tfeb) <- names(counts.pb.m)


# loop over pseudobulk data
for  (i in 1:length(counts.pb.m)) {
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Pseudobulk_Male_",names(counts.pb.m)[[i]])
  
  ###################################
  #######       Run DE      #######
  
  # legend
  my.cols  <- rep("",nrow(my.Male.meta))
  my.cols[my.Male.meta$Treat.gp %in% "SED"]  <- "deepskyblue"
  my.cols[my.Male.meta$Treat.gp %in% "RUN"]  <- "mediumturquoise"
  my.cols[my.Male.meta$Treat.gp %in% "TFEB"] <- "royalblue3"
  
  # get matrix using age as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = counts.pb.m[[i]],
                                colData   = my.Male.meta,
                                design    = ~ Treat.gp)
  
  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds)
  
  # plot dispersion
  my.disp.out <- paste(my.outprefix,"_dispersion_plot.pdf")
  
  pdf(my.disp.out)
  plotDispEsts(dds.deseq)
  dev.off()
  
  # get DESeq2 normalized expression value
  vst.cts.m[[i]] <- getVarianceStabilizedData(dds.deseq)
  
  # MDS analysis
  mds.result <- cmdscale(1-cor(vst.cts.m[[i]],method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  x <- mds.result[, 1]
  y <- mds.result[, 2]
  
  pdf(paste0(my.outprefix,"_MDS_plot.pdf"))
  plot(x, y,
       xlab = "MDS dimension 1", ylab = "MDS dimension 2",
       main= paste0(names(counts.pb.m)[[i]]," (MDS)"),
       cex=3, col = my.cols, pch = 16,
       cex.lab = 1.25,
       cex.axis = 1.25, las = 1)
  legend("bottomright", c("Sedentary","Runner","TFEb"),col = c("deepskyblue","mediumturquoise","royalblue3"), pch = 16, bty = 'n')
  dev.off()
  
  # extract gene significance by DEseq2
  res.run  <- results(dds.deseq, contrast = c("Treat.gp", "RUN" , "SED")) 
  res.tfeb <- results(dds.deseq, contrast = c("Treat.gp", "TFEB", "SED")) 
  
  # exclude genes with NA FDR value
  res.run  <- res.run [!is.na(res.run $padj),]
  res.tfeb <- res.tfeb[!is.na(res.tfeb$padj),]
  
  # store results
  deseq.res.list.m.run[[i]]        <- data.frame(res.run  )
  deseq.res.list.m.tfeb[[i]]       <- data.frame(res.tfeb )
  
  ### get sex dimorphic changes at FDR5
  genes.run  <- rownames(res.run)[res.run$padj < 0.05]
  genes.tfeb <- rownames(res.tfeb)[res.tfeb$padj < 0.05]
  
  if (length(genes.run) > 2) {
    # heatmap drawing - only if there is at least 2 gene
    my.heatmap.out <- paste0(my.outprefix,"_RUNNING_Heatmap_FDR5_GENES.pdf")
    
    pdf(my.heatmap.out, onefile = F, height = 10, width = 10)
    my.heatmap.title <- paste0(names(counts.pb.m)[[i]], " running significant (FDR<5%), ", length(genes.run), " genes")
    pheatmap::pheatmap(vst.cts.m[[i]][genes.run,],
                       cluster_cols = F,
                       cluster_rows = T,
                       colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                       show_rownames = F, scale="row",
                       main = my.heatmap.title,
                       cellwidth = 15,
                       border    = NA,
                       cellheight = 0.5)
    dev.off()
  }
  
  if (length(genes.tfeb) > 2) {
    # heatmap drawing - only if there is at least 2 gene
    my.heatmap.out <- paste0(my.outprefix,"_TFEB_Heatmap_FDR5_GENES.pdf")
    
    pdf(my.heatmap.out, onefile = F, height = 10, width = 10)
    my.heatmap.title <- paste0(names(counts.pb.m)[[i]], " TFEB significant (FDR<5%), ", length(genes.tfeb), " genes")
    pheatmap::pheatmap(vst.cts.m[[i]][genes.tfeb,],
                       cluster_cols = F,
                       cluster_rows = T,
                       colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                       show_rownames = F, scale="row",
                       main = my.heatmap.title,
                       cellwidth = 15,
                       border    = NA,
                       cellheight = 0.5)
    dev.off()
  }
  
  
  ## compare
  my.merged <- merge(data.frame(res.run), data.frame(res.tfeb), by = "row.names", suffixes = c(".run",".tfeb"))
  
  my.test <- cor.test(my.merged$log2FoldChange.run, my.merged$log2FoldChange.tfeb, method = "spearman")
  my.rho  <- signif(my.test$estimate,3)
  
  pdf(paste0(Sys.Date(),"_", names(counts.pb.m)[[i]], "_running_vs_TFEB_FC_scatterplot_Male_FDR5.pdf"), height = 10, width = 10)
  smoothScatter(my.merged$log2FoldChange.run, my.merged$log2FoldChange.tfeb,
                xlim = c(-8, 8), ylim = c(-8, 8),
                xlab = "Log2FC (Running/Sedentary)",
                ylab = "Log2FC (TFEB/Sedentary)",
                main = paste0("Male ", names(counts.pb.m)[[i]]) )
  abline(0, 1,  col = "red", lty = "dashed")
  abline(h = 0, col = "grey", lty = "dashed")
  abline(v = 0, col = "grey", lty = "dashed")
  text(-7, 7.5, paste("Rho = ", my.rho), pos = 4)
  text(-7, 6.5, paste("p = ",signif(my.test$p.value,2)), pos = 4)
  dev.off()
  
}

# save R object with all DEseq2 results
my.rdata <- paste0(Sys.Date(),"_pseudobulk_Muscle_cell_types_RUN_TFEB_Male_DEseq2_objects.RData")
save(deseq.res.list.m.run, deseq.res.list.m.tfeb, file = my.rdata)

my.vst <- paste0(Sys.Date(),"_pseudobulk_Muscle_cell_types_RUN_TFEB_Male_VST_data_objects.RData")
save(vst.cts.m, file = my.vst)
##########################################################################################################################################

##########################################################################################################################################
## 3. Make excel workbooks with results
load('2024-12-13_pseudobulk_Muscle_cell_types_RUN_TFEB_Female_DEseq2_objects.RData')
load('2024-12-13_pseudobulk_Muscle_cell_types_RUN_TFEB_Male_DEseq2_objects.RData')

options(java.parameters = "-Xmx16g" )
require(openxlsx)

# "deseq.res.list.f.run"  "deseq.res.list.f.tfeb" "deseq.res.list.m.run"  "deseq.res.list.m.tfeb" 

names(deseq.res.list.f.run )  <- paste0("F_running_",names(deseq.res.list.f.run ))
names(deseq.res.list.f.tfeb)  <- paste0("F_TFEb_"   ,names(deseq.res.list.f.tfeb))
names(deseq.res.list.m.run )  <- paste0("M_running_",names(deseq.res.list.m.run ))
names(deseq.res.list.m.tfeb)  <- paste0("M_TFEb_"   ,names(deseq.res.list.m.tfeb))

deseq.out <-  paste0(Sys.Date(),"_Running_TFeb_DESeq2_Muscle_Results.xlsx")

write.xlsx(c(deseq.res.list.f.run,deseq.res.list.f.tfeb,deseq.res.list.m.run,deseq.res.list.m.tfeb), rowNames = TRUE, file = deseq.out)
##########################################################################################################################################

##########################################################################################################################################
## 4.get some summary plots
load('2024-12-13_pseudobulk_Muscle_cell_types_RUN_TFEB_Female_DEseq2_objects.RData')
load('2024-12-13_pseudobulk_Muscle_cell_types_RUN_TFEB_Female_VST_data_objects.RData')
load('2024-12-13_pseudobulk_Muscle_cell_types_RUN_TFEB_Male_DEseq2_objects.RData')
load('2024-12-13_pseudobulk_Muscle_cell_types_RUN_TFEB_Male_VST_data_objects.RData')

# "deseq.res.list.f.run"  "deseq.res.list.f.tfeb" "deseq.res.list.m.run"  "deseq.res.list.m.tfeb" "vst.cts.f"             "vst.cts.m"            


######### A. Calculate correlations
# get cell types
cell.types <- names(deseq.res.list.f.run)

# prepare lists for merged data
deseq.res.F.merged <- vector(mode = "list", length = length(cell.types))
deseq.res.M.merged <- vector(mode = "list", length = length(cell.types))

# get table of spearman correlation
corr.tab              <- data.frame("Cell_Type" = rep(cell.types,2),
                                    "Sex"       = c(rep("Females",length(cell.types)), rep("Males",length(cell.types))))
corr.tab$Spearman_Rho <- NA
corr.tab$pval         <- NA


for (i in 1:length(cell.types)) {
  
  ##### a. female data
  deseq.res.F.merged[[i]] <- merge(data.frame(deseq.res.list.f.run[[i]]), data.frame(deseq.res.list.f.tfeb[[i]]), by = "row.names", suffixes = c(".run",".tfeb"))
  my.test <- cor.test(deseq.res.F.merged[[i]]$log2FoldChange.run, deseq.res.F.merged[[i]]$log2FoldChange.tfeb, method = "spearman")
  
  # get values
  corr.tab$Spearman_Rho[i]  <- signif(my.test$estimate,3)
  corr.tab$pval[i]          <- signif(my.test$p.value,3)
  
  ##### b. male data
  deseq.res.M.merged[[i]] <- merge(data.frame(deseq.res.list.m.run[[i]]), data.frame(deseq.res.list.m.tfeb[[i]]), by = "row.names", suffixes = c(".run",".tfeb"))
  my.test <- cor.test(deseq.res.M.merged[[i]]$log2FoldChange.run, deseq.res.M.merged[[i]]$log2FoldChange.tfeb, method = "spearman")
  
  # get values
  corr.tab$Spearman_Rho[i+length(cell.types)]  <- signif(my.test$estimate,3)
  corr.tab$pval[i+length(cell.types)]          <- signif(my.test$p.value,3)
  
}

write.table(corr.tab, file = paste0(Sys.Date(),"_Running_TFeb_correlation_analyses.txt"), sep = "\t", row.names = F, quote = F)


######### B. Plot correlations
corr.mat.plot <- data.frame(row.names = cell.types,
                            "Female"  = corr.tab$Spearman_Rho[1:6],
                            "Male"    = corr.tab$Spearman_Rho[7:12])

col_fun = colorRamp2(c(0, 0.5, 1), c("white", "firebrick3", "firebrick4"))

pdf(paste0(Sys.Date(),"_Running_TFEB_Spearman_Rho_Heatmap.pdf"), onefile = F, height = 7, width = 5.5)
Heatmap(corr.mat.plot, 
        name = "Rho", 
        col = col_fun, 
        column_title = "Rho - Running vs TFeb",
        rect_gp = gpar(col = "grey", lwd = 0.5),
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(Sex = c("Female","Male"), col = list(Sex = c("Male" = "deepskyblue", "Female" = "deeppink"))) )
dev.off()
##########################################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_MuscatDEseq2_PB_DESeq2_scRNAseq_Muscle_session_Info.txt", sep =""))
sessionInfo()
sink()
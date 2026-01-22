setwd('/Volumes/BB_HQ_4/Collaboration/Cortes_lab/Muscle_snRNAseq/Preprocessing')
options(stringsAsFactors = F)
options (future.globals.maxSize = 32000 * 1024^2)

# General use packages
library('Seurat')
library(bitops)
library(sctransform)

# removal of ambient RNA with DecontX
# https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html#seurat
library('celda')

# Doublet identification packages
# https://github.com/chris-mcginnis-ucsf/DoubletFinder
library(DoubletFinder)

# cxds_bcds_hybrid related packages
# https://www.bioconductor.org/packages/release/bioc/vignettes/scds/inst/doc/scds.html
library(scds)
library(scater)
library(bitops)


################################################################################################
# 2024-11-29
# Processing of all snRNAseq muscle data for Cortes lab
#      - use decontX to diminish the impact of ambient RNA
#      - clean doublets, poor QC cells, etc in each library
#      - use intersection of "generous" run of scds and doubletfinder (since we know nuclei are more susceptible)
#      - decontX filter to remove any cell with predicted contamination level ≥ 0.25
################################################################################################

################################################################################################
#### 0. Assumed doublet information/to calculate %age for prediction
# Targeted Cell Recovery  # of Cells Loaded	Barcodes Detected	Singlets Multiplets	Multiplet Rate
# 3,000	                         4,950           ~3,000        ~2,900	   ~80	    ~2.4%
# 4,000	                         6,600           ~3,900	       ~3,800	   ~140     ~3.2%
# 5,000	                         8,250           ~4,800        ~4,600	   ~210     ~4.0%
# 6,000                          9,900           ~5,700	       ~5,400	   ~300	    ~4.8%
# 7,000                          11,550	         ~6,600	       ~6,200	   ~400	    ~5.6%
# 8,000                          13,200	         ~7,500	       ~7,000	   ~510	    ~6.4%
# 9,000                          14,850	         ~8,400	       ~7,700	   ~640	    ~7.2%
#10,000                          16,500	         ~9,200	       ~8,400	   ~780	    ~8.0%
#12,000                          19,800	         ~10,900	     ~9,800	   ~1,100   ~9.6%

pred.10x.dblt <- data.frame( "cell_number" = c(3000,4000,5000,6000,7000,8000,9000, 10000, 12000),
                             "dblt_rate"   = c(2.4 ,3.2 ,4.0 ,4.8 ,5.6 ,6.4 ,7.2 , 8.0  , 9.6))

pred_dblt_lm <- lm(dblt_rate ~ cell_number, data = pred.10x.dblt)

pdf(paste0(Sys.Date(), "_10x_cel_number_vs_doublet_rate.pdf"))
plot(dblt_rate ~ cell_number, data = pred.10x.dblt)
abline(pred_dblt_lm, col = "red", lty = "dashed")
dev.off()

# predict(pred_dblt_lm, data.frame("cell_number" = 8500))
################################################################################################


################################################################################################
#### 1a. Read data from CellRanger, perform background removal

# Calculate and clean the contribution of ambient RNA with DecontX
# read 10x libraries cell ranger gene barcode matrices for DecontX
cts.F_RUN_1  <- Read10X('../cellranger/Female_6mo_runner_1/outs/filtered_feature_bc_matrix/')
cts.F_RUN_2  <- Read10X('../cellranger/Female_6mo_runner_2/outs/filtered_feature_bc_matrix/')
cts.F_SED_1  <- Read10X('../cellranger/Female_6mo_sed_1/outs/filtered_feature_bc_matrix/')
cts.F_SED_2  <- Read10X('../cellranger/Female_6mo_sed_2/outs/filtered_feature_bc_matrix/')
cts.F_TFB_1  <- Read10X('../cellranger/Female_6mo_TFEB_1/outs/filtered_feature_bc_matrix/')
cts.F_TFB_2  <- Read10X('../cellranger/Female_6mo_TFEB_2/outs/filtered_feature_bc_matrix/')
cts.M_RUN_1  <- Read10X('../cellranger/Male_6mo_runner_1/outs/filtered_feature_bc_matrix/')
cts.M_RUN_2  <- Read10X('../cellranger/Male_6mo_runner_2/outs/filtered_feature_bc_matrix/')
cts.M_SED_1  <- Read10X('../cellranger/Male_6mo_sed_1/outs/filtered_feature_bc_matrix/')
cts.M_SED_2  <- Read10X('../cellranger/Male_6mo_sed_2/outs/filtered_feature_bc_matrix/')
cts.M_TFB_1  <- Read10X('../cellranger/Male_6mo_TFEB_1/outs/filtered_feature_bc_matrix/')
cts.M_TFB_2  <- Read10X('../cellranger/Male_6mo_TFEB_2/outs/filtered_feature_bc_matrix/')

cts_raw.F_RUN_1  <- Read10X('../cellranger/Female_6mo_runner_1/outs/raw_feature_bc_matrix/')
cts_raw.F_RUN_2  <- Read10X('../cellranger/Female_6mo_runner_2/outs/raw_feature_bc_matrix/')
cts_raw.F_SED_1  <- Read10X('../cellranger/Female_6mo_sed_1/outs/raw_feature_bc_matrix/')
cts_raw.F_SED_2  <- Read10X('../cellranger/Female_6mo_sed_2/outs/raw_feature_bc_matrix/')
cts_raw.F_TFB_1  <- Read10X('../cellranger/Female_6mo_TFEB_1/outs/raw_feature_bc_matrix/')
cts_raw.F_TFB_2  <- Read10X('../cellranger/Female_6mo_TFEB_2/outs/raw_feature_bc_matrix/')
cts_raw.M_RUN_1  <- Read10X('../cellranger/Male_6mo_runner_1/outs/raw_feature_bc_matrix/')
cts_raw.M_RUN_2  <- Read10X('../cellranger/Male_6mo_runner_2/outs/raw_feature_bc_matrix/')
cts_raw.M_SED_1  <- Read10X('../cellranger/Male_6mo_sed_1/outs/raw_feature_bc_matrix/')
cts_raw.M_SED_2  <- Read10X('../cellranger/Male_6mo_sed_2/outs/raw_feature_bc_matrix/')
cts_raw.M_TFB_1  <- Read10X('../cellranger/Male_6mo_TFEB_1/outs/raw_feature_bc_matrix/')
cts_raw.M_TFB_2  <- Read10X('../cellranger/Male_6mo_TFEB_2/outs/raw_feature_bc_matrix/')

# Create SingleCellExperiment objects
sce.F_RUN_1     <- SingleCellExperiment(list(counts = cts.F_RUN_1    ))
sce.F_RUN_2     <- SingleCellExperiment(list(counts = cts.F_RUN_2    ))
sce.F_SED_1     <- SingleCellExperiment(list(counts = cts.F_SED_1    ))
sce.F_SED_2     <- SingleCellExperiment(list(counts = cts.F_SED_2    ))
sce.F_TFB_1     <- SingleCellExperiment(list(counts = cts.F_TFB_1    ))
sce.F_TFB_2     <- SingleCellExperiment(list(counts = cts.F_TFB_2    ))
sce.M_RUN_1     <- SingleCellExperiment(list(counts = cts.M_RUN_1    ))
sce.M_RUN_2     <- SingleCellExperiment(list(counts = cts.M_RUN_2    ))
sce.M_SED_1     <- SingleCellExperiment(list(counts = cts.M_SED_1    ))
sce.M_SED_2     <- SingleCellExperiment(list(counts = cts.M_SED_2    ))
sce.M_TFB_1     <- SingleCellExperiment(list(counts = cts.M_TFB_1    ))
sce.M_TFB_2     <- SingleCellExperiment(list(counts = cts.M_TFB_2    ))

sce_raw.F_RUN_1 <- SingleCellExperiment(list(counts = cts_raw.F_RUN_1    ))
sce_raw.F_RUN_2 <- SingleCellExperiment(list(counts = cts_raw.F_RUN_2    ))
sce_raw.F_SED_1 <- SingleCellExperiment(list(counts = cts_raw.F_SED_1    ))
sce_raw.F_SED_2 <- SingleCellExperiment(list(counts = cts_raw.F_SED_2    ))
sce_raw.F_TFB_1 <- SingleCellExperiment(list(counts = cts_raw.F_TFB_1    ))
sce_raw.F_TFB_2 <- SingleCellExperiment(list(counts = cts_raw.F_TFB_2    ))
sce_raw.M_RUN_1 <- SingleCellExperiment(list(counts = cts_raw.M_RUN_1    ))
sce_raw.M_RUN_2 <- SingleCellExperiment(list(counts = cts_raw.M_RUN_2    ))
sce_raw.M_SED_1 <- SingleCellExperiment(list(counts = cts_raw.M_SED_1    ))
sce_raw.M_SED_2 <- SingleCellExperiment(list(counts = cts_raw.M_SED_2    ))
sce_raw.M_TFB_1 <- SingleCellExperiment(list(counts = cts_raw.M_TFB_1    ))
sce_raw.M_TFB_2 <- SingleCellExperiment(list(counts = cts_raw.M_TFB_2    ))

# Run decontX
sce.F_RUN_1  <- decontX(sce.F_RUN_1, background = sce_raw.F_RUN_1 )
sce.F_RUN_2  <- decontX(sce.F_RUN_2, background = sce_raw.F_RUN_2 )
sce.F_SED_1  <- decontX(sce.F_SED_1, background = sce_raw.F_SED_1 )
sce.F_SED_2  <- decontX(sce.F_SED_2, background = sce_raw.F_SED_2 )
sce.F_TFB_1  <- decontX(sce.F_TFB_1, background = sce_raw.F_TFB_1 )
sce.F_TFB_2  <- decontX(sce.F_TFB_2, background = sce_raw.F_TFB_2 )
sce.M_RUN_1  <- decontX(sce.M_RUN_1, background = sce_raw.M_RUN_1 )
sce.M_RUN_2  <- decontX(sce.M_RUN_2, background = sce_raw.M_RUN_2 )
sce.M_SED_1  <- decontX(sce.M_SED_1, background = sce_raw.M_SED_1 )
sce.M_SED_2  <- decontX(sce.M_SED_2, background = sce_raw.M_SED_2 )
sce.M_TFB_1  <- decontX(sce.M_TFB_1, background = sce_raw.M_TFB_1 )
sce.M_TFB_2  <- decontX(sce.M_TFB_2, background = sce_raw.M_TFB_2 )

# get seurat objects
seurat.F_RUN_1 <- CreateSeuratObject( round(decontXcounts(sce.F_RUN_1)) )
seurat.F_RUN_2 <- CreateSeuratObject( round(decontXcounts(sce.F_RUN_2)) )
seurat.F_SED_1 <- CreateSeuratObject( round(decontXcounts(sce.F_SED_1)) )
seurat.F_SED_2 <- CreateSeuratObject( round(decontXcounts(sce.F_SED_2)) )
seurat.F_TFB_1 <- CreateSeuratObject( round(decontXcounts(sce.F_TFB_1)) )
seurat.F_TFB_2 <- CreateSeuratObject( round(decontXcounts(sce.F_TFB_2)) )
seurat.M_RUN_1 <- CreateSeuratObject( round(decontXcounts(sce.M_RUN_1)) )
seurat.M_RUN_2 <- CreateSeuratObject( round(decontXcounts(sce.M_RUN_2)) )
seurat.M_SED_1 <- CreateSeuratObject( round(decontXcounts(sce.M_SED_1)) )
seurat.M_SED_2 <- CreateSeuratObject( round(decontXcounts(sce.M_SED_2)) )
seurat.M_TFB_1 <- CreateSeuratObject( round(decontXcounts(sce.M_TFB_1)) )
seurat.M_TFB_2 <- CreateSeuratObject( round(decontXcounts(sce.M_TFB_2)) )

seurat.F_RUN_1 <- AddMetaData(object = seurat.F_RUN_1, colData(sce.F_RUN_1)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.F_RUN_2 <- AddMetaData(object = seurat.F_RUN_2, colData(sce.F_RUN_2)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.F_SED_1 <- AddMetaData(object = seurat.F_SED_1, colData(sce.F_SED_1)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.F_SED_2 <- AddMetaData(object = seurat.F_SED_2, colData(sce.F_SED_2)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.F_TFB_1 <- AddMetaData(object = seurat.F_TFB_1, colData(sce.F_TFB_1)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.F_TFB_2 <- AddMetaData(object = seurat.F_TFB_2, colData(sce.F_TFB_2)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.M_RUN_1 <- AddMetaData(object = seurat.M_RUN_1, colData(sce.M_RUN_1)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.M_RUN_2 <- AddMetaData(object = seurat.M_RUN_2, colData(sce.M_RUN_2)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.M_SED_1 <- AddMetaData(object = seurat.M_SED_1, colData(sce.M_SED_1)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.M_SED_2 <- AddMetaData(object = seurat.M_SED_2, colData(sce.M_SED_2)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.M_TFB_1 <- AddMetaData(object = seurat.M_TFB_1, colData(sce.M_TFB_1)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.M_TFB_2 <- AddMetaData(object = seurat.M_TFB_2, colData(sce.M_TFB_2)$decontX_contamination  , col.name = "decontX_contamination"  )

# Merge objects for the cohort
muscle.combined <- merge(seurat.F_RUN_1,
                           y =  c(seurat.F_RUN_2,
                                  seurat.F_SED_1,
                                  seurat.F_SED_2,
                                  seurat.F_TFB_1,
                                  seurat.F_TFB_2,
                                  seurat.M_RUN_1,
                                  seurat.M_RUN_2,
                                  seurat.M_SED_1,
                                  seurat.M_SED_2,
                                  seurat.M_TFB_1,
                                  seurat.M_TFB_2),
                           add.cell.ids = c("F_RUN_1"  ,
                                            "F_RUN_2"  ,
                                            "F_SED_1"  ,
                                            "F_SED_2"  ,
                                            "F_TFB_1"  ,
                                            "F_TFB_2"  ,
                                            "M_RUN_1"  ,
                                            "M_RUN_2"  ,
                                            "M_SED_1"  ,
                                            "M_SED_2"  ,
                                            "M_TFB_1"  ,
                                            "M_TFB_2"  ),
                           project = "10x_Muscle_exercise")
muscle.combined
# An object of class Seurat 
# 33696 features across 70699 samples within 1 assay 
# Active assay: RNA (33696 features, 0 variable features)
# 2 layers present: counts, data

#clean memory
rm(cts.F_RUN_1, cts.F_RUN_2, cts.F_SED_1, cts.F_SED_2, cts.F_TFB_1, cts.F_TFB_2, cts.M_RUN_1, cts.M_RUN_2, cts.M_SED_1, cts.M_SED_2, cts.M_TFB_1, cts.M_TFB_2, cts_raw.F_RUN_1, cts_raw.F_RUN_2, cts_raw.F_SED_1, cts_raw.F_SED_2, cts_raw.F_TFB_1, cts_raw.F_TFB_2, cts_raw.M_RUN_1, cts_raw.M_RUN_2, cts_raw.M_SED_1, cts_raw.M_SED_2, cts_raw.M_TFB_1, cts_raw.M_TFB_2, sce.F_RUN_1, sce.F_RUN_2, sce.F_SED_1, sce.F_SED_2, sce.F_TFB_1, sce.F_TFB_2, sce.M_RUN_1, sce.M_RUN_2, sce.M_SED_1, sce.M_SED_2, sce.M_TFB_1, sce.M_TFB_2, sce_raw.F_RUN_1, sce_raw.F_RUN_2, sce_raw.F_SED_1, sce_raw.F_SED_2, sce_raw.F_TFB_1, sce_raw.F_TFB_2, sce_raw.M_RUN_1, sce_raw.M_RUN_2, sce_raw.M_SED_1, sce_raw.M_SED_2, sce_raw.M_TFB_1, sce_raw.M_TFB_2, seurat.F_RUN_1, seurat.F_RUN_2, seurat.F_SED_1, seurat.F_SED_2, seurat.F_TFB_1, seurat.F_TFB_2, seurat.M_RUN_1, seurat.M_RUN_2, seurat.M_SED_1, seurat.M_SED_2, seurat.M_TFB_1, seurat.M_TFB_2)
################################################################################################


################################################################################################
#### 1b. Add key metadata to Seurat object

# create Group label
my.F_RUN_1   <- grep("F_RUN_1"  , colnames(muscle.combined@assays$RNA))
my.F_RUN_2   <- grep("F_RUN_2"  , colnames(muscle.combined@assays$RNA))
my.F_SED_1   <- grep("F_SED_1"  , colnames(muscle.combined@assays$RNA))
my.F_SED_2   <- grep("F_SED_2"  , colnames(muscle.combined@assays$RNA))
my.F_TFB_1   <- grep("F_TFB_1"  , colnames(muscle.combined@assays$RNA))
my.F_TFB_2   <- grep("F_TFB_2"  , colnames(muscle.combined@assays$RNA))
my.M_RUN_1   <- grep("M_RUN_1"  , colnames(muscle.combined@assays$RNA))
my.M_RUN_2   <- grep("M_RUN_2"  , colnames(muscle.combined@assays$RNA))
my.M_SED_1   <- grep("M_SED_1"  , colnames(muscle.combined@assays$RNA))
my.M_SED_2   <- grep("M_SED_2"  , colnames(muscle.combined@assays$RNA))
my.M_TFB_1   <- grep("M_TFB_1"  , colnames(muscle.combined@assays$RNA))
my.M_TFB_2   <- grep("M_TFB_2"  , colnames(muscle.combined@assays$RNA))

#####
Group <- rep("NA", length(colnames(muscle.combined@assays$RNA)))
Group[ my.F_RUN_1  ]   <- "F_runner"
Group[ my.F_RUN_2  ]   <- "F_runner"
Group[ my.F_SED_1  ]   <- "F_sedentary"
Group[ my.F_SED_2  ]   <- "F_sedentary"
Group[ my.F_TFB_1  ]   <- "F_TFEB"
Group[ my.F_TFB_2  ]   <- "F_TFEB"
Group[ my.M_RUN_1  ]   <- "M_runner"
Group[ my.M_RUN_2  ]   <- "M_runner"
Group[ my.M_SED_1  ]   <- "M_sedentary"
Group[ my.M_SED_2  ]   <- "M_sedentary"
Group[ my.M_TFB_1  ]   <- "M_TFEB"
Group[ my.M_TFB_2  ]   <- "M_TFEB"
Group <- data.frame(Group)
rownames(Group) <- colnames(muscle.combined@assays$RNA)

#####
Sex <- rep("NA", length(colnames(muscle.combined@assays$RNA)))
Sex[ my.F_RUN_1 ]   <- "F"  
Sex[ my.F_RUN_2 ]   <- "F" 
Sex[ my.F_SED_1 ]   <- "F" 
Sex[ my.F_SED_2 ]   <- "F"
Sex[ my.F_TFB_1 ]   <- "F" 
Sex[ my.F_TFB_2 ]   <- "F"  
Sex[ my.M_RUN_1 ]   <- "M"   
Sex[ my.M_RUN_2 ]   <- "M" 
Sex[ my.M_SED_1 ]   <- "M" 
Sex[ my.M_SED_2 ]   <- "M"  
Sex[ my.M_TFB_1 ]   <- "M"
Sex[ my.M_TFB_2 ]   <- "M" 
Sex <- data.frame(Sex)
rownames(Sex) <- colnames(muscle.combined@assays$RNA)

##### Treatment
Treat.gp <- rep("NA", length(colnames(muscle.combined@assays$RNA)))
Treat.gp[ my.F_RUN_1 ]   <- "RUN"
Treat.gp[ my.F_RUN_2 ]   <- "RUN"
Treat.gp[ my.F_SED_1 ]   <- "SED"
Treat.gp[ my.F_SED_2 ]   <- "SED"
Treat.gp[ my.F_TFB_1 ]   <- "TFEB"
Treat.gp[ my.F_TFB_2 ]   <- "TFEB"
Treat.gp[ my.M_RUN_1 ]   <- "RUN"
Treat.gp[ my.M_RUN_2 ]   <- "RUN"
Treat.gp[ my.M_SED_1 ]   <- "SED"
Treat.gp[ my.M_SED_2 ]   <- "SED"
Treat.gp[ my.M_TFB_1 ]   <- "TFEB"
Treat.gp[ my.M_TFB_2 ]   <- "TFEB"
Treat.gp <- data.frame(Treat.gp)
rownames(Treat.gp) <- colnames(muscle.combined@assays$RNA)

#####
Sample <- rep("NA", length(colnames(muscle.combined@assays$RNA)))
Sample[ my.F_RUN_1 ]   <- "F_RUN_1"  
Sample[ my.F_RUN_2 ]   <- "F_RUN_2" 
Sample[ my.F_SED_1 ]   <- "F_SED_1" 
Sample[ my.F_SED_2 ]   <- "F_SED_2"
Sample[ my.F_TFB_1 ]   <- "F_TFB_1" 
Sample[ my.F_TFB_2 ]   <- "F_TFB_2"  
Sample[ my.M_RUN_1 ]   <- "M_RUN_1"   
Sample[ my.M_RUN_2 ]   <- "M_RUN_2" 
Sample[ my.M_SED_1 ]   <- "M_SED_1" 
Sample[ my.M_SED_2 ]   <- "M_SED_2"  
Sample[ my.M_TFB_1 ]   <- "M_TFB_1"
Sample[ my.M_TFB_2 ]   <- "M_TFB_2" 
Sample <- data.frame(Sample)
rownames(Sample) <- colnames(muscle.combined@assays$RNA)


# update Seurat with metadata
muscle.combined <- AddMetaData(object = muscle.combined, metadata = as.vector(Group)       , col.name = "Group"       )
muscle.combined <- AddMetaData(object = muscle.combined, metadata = as.vector(Sex)         , col.name = "Sex"         )
muscle.combined <- AddMetaData(object = muscle.combined, metadata = as.vector(Treat.gp)    , col.name = "Treat.gp"   )
muscle.combined <- AddMetaData(object = muscle.combined, metadata = as.vector(Sample)      , col.name = "Sample"       )
################################################################################################


################################################################################################
#### 2. Basic QC and filtering with Seurat and decontX

### No filtering on genes at this stage - only after all cohorts merged for fairness
muscle.combined <- SetIdent(muscle.combined, value = "Group")

# DecontX contamination levels for filtration
pdf(paste(Sys.Date(),"Muscle_data_Cortes_lab_violinPlots_QC_DecontX.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = muscle.combined, features = c("decontX_contamination"), pt.size = 0)
dev.off()

# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
muscle.combined[["percent.mito"]] <- PercentageFeatureSet(muscle.combined, pattern = "^mt-")
head(muscle.combined@meta.data)

pdf(paste(Sys.Date(),"Muscle_data_Cortes_lab_violinPlots_QC_gene_UMI_mito.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = muscle.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(muscle.combined, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(muscle.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste(Sys.Date(),"Muscle_data_Cortes_lab_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
plot1 + plot2
dev.off()

# filter dead/low Q cells
muscle.combined <- subset(muscle.combined, subset = nFeature_RNA > 250 & nFeature_RNA < 5000 & percent.mito < 15 & nCount_RNA < 10000 & decontX_contamination < 0.25 )
muscle.combined
# An object of class Seurat 
# 33696 features across 66339 samples within 1 assay 
# Active assay: RNA (33696 features, 0 variable features)
# 2 layers present: counts, data

### Check data after cell filtering
head(muscle.combined@meta.data)

table(muscle.combined@meta.data$Group)
# F_runner F_sedentary      F_TFEB    M_runner M_sedentary      M_TFEB 
# 11280       10696        9536       14273       12794        7760 

#### Normalize the data for doublet analysis, etc
# global-scaling normalization method 'LogNormalize' normalizes gene expression measurements for each cell 
# by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
muscle.combined <- NormalizeData(object = muscle.combined, normalization.method = "LogNormalize",  scale.factor = 10000)
################################################################################################

################################################################################################
#### 2b. Cell cycle prediction and storage
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "../cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")

# make into mouse gene names
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

cc.genes.mouse <- firstup(tolower(cc.genes))

# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes   <- cc.genes.mouse[1:43]
g2m.genes <- cc.genes.mouse[44:97]

# Assign Cell-Cycle Scores
muscle.combined <- CellCycleScoring(object = muscle.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# write predictions to file
write.table(muscle.combined@meta.data, file = paste0(Sys.Date(),"_Muscle_data_Cortes_lab_CellCycle_predictions.txt"), sep = "\t", quote = F)

table(muscle.combined@meta.data$Group, muscle.combined$Phase)
#               G1  G2M    S
# F_runner    4647 3088 3545
# F_sedentary 4365 2850 3481
# F_TFEB      4019 2450 3067
# M_runner    5811 3943 4519
# M_sedentary 4939 3473 4382
# M_TFEB      2809 1997 2954
################################################################################################

################################################################################################
#### 3. Find and remove doublets using doublet finder & scds workflow

# https://github.com/chris-mcginnis-ucsf/DoubletFinder
muscle.combined <- SCTransform(object = muscle.combined, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mito", "Phase"))
save(muscle.combined, file = paste0(Sys.Date(),"_Muscle_data_Cortes_lab_Seurat_object_postSCT.RData"))

# Run first pass analysis just for doublet identification (not final clustering)
muscle.combined <- RunPCA(muscle.combined, npcs = 30)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_Muscle_data_Cortes_lab_ElbowPlot.pdf"))
ElbowPlot(muscle.combined, ndims = 30)
dev.off()

# run dimensionality reduction
# Keep all PCs here, we'll do the clean clustering analysis on the merged object across all cohorts
muscle.combined <- RunUMAP(muscle.combined, dims = 1:30)
muscle.combined <- FindNeighbors(muscle.combined, dims = 1:30)
muscle.combined <- FindClusters(object = muscle.combined)

#### need to split by 10x sample to make sure to identify real doublets
# will run on one object at a time
cohort.list <- SplitObject(muscle.combined, split.by = "Sample")

## Assume doublet rate based on 10x information (add a 15% fudge factor due to nuclei being more sticky)
pred.dblt.rate <- 1.15 * predict(pred_dblt_lm, data.frame("cell_number" = unlist(lapply(cohort.list, ncol))))/100


## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& A. Run DoubletFinder
# loop over samples
for (i in 1:length(cohort.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list_musc <- paramSweep_v3(cohort.list[[i]], PCs = 1:30, sct = TRUE, num.cores	 = 4)
  sweep.stats_musc    <- summarizeSweep(sweep.res.list_musc, GT = FALSE)
  bcmvn_musc          <- find.pK(sweep.stats_musc)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yield row number
  pk.killi <- as.numeric(as.character(bcmvn_musc[as.numeric(bcmvn_musc$pK[bcmvn_musc$BCmetric == max(bcmvn_musc$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(cohort.list[[i]]@meta.data$seurat_clusters)             ## ex: annotations
  nExp_poi       <- round((pred.dblt.rate[i]) *length(cohort.list[[i]]@meta.data$Group))      ## Assume doublets based on nuclei isolation protocol performance
  
  ## Run DoubletFinder with varying classification stringencies
  cohort.list[[i]] <- doubletFinder_v3(cohort.list[[i]], PCs = 1:30, pN = 0.25, pK = pk.killi, nExp = nExp_poi,     reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(cohort.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(cohort.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(cohort.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(cohort.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(cohort.list)) {
  
  pdf(paste(Sys.Date(),"Muscle",names(cohort.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(cohort.list[[i]], reduction = "umap", group.by = "DoubletFinder"), raster = T)
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
singlets.annot <- merge(cohort.list[[1]],
                                 y = c(cohort.list[[2]],
                                       cohort.list[[3]],
                                       cohort.list[[4]],
                                       cohort.list[[5]],
                                       cohort.list[[6]],
                                       cohort.list[[7]],
                                       cohort.list[[8]],
                                       cohort.list[[9]],
                                       cohort.list[[10]],
                                       cohort.list[[11]],
                                       cohort.list[[12]]
                                       ),
                                 project = "Muscle_data_Cortes_lab")
singlets.annot


# remove pANN columns that are 10xGenomics library lane specific
singlets.annot@meta.data <- singlets.annot@meta.data[,-grep("pANN",colnames(singlets.annot@meta.data))]


## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& B. Run scds:single cell doublet scoring (hybrid method)
# cxds is based on co-expression of gene pairs and works with absence/presence calls only, 
# bcds uses the full count information and a binary classification approach using artificially generated doublets. 
# cxds_bcds_hybrid combines both approaches

# create scds working object - convert list to SingleCellExperiment
cohort.list.scds        <- lapply(cohort.list, as.SingleCellExperiment)

# loop over sample
for (i in 1:length(cohort.list.scds)) {
  
  # Annotate doublets using co-expression based doublet scoring:
  cohort.list.scds[[i]] <- cxds_bcds_hybrid(cohort.list.scds[[i]])
  
  # predicted doublet rate
  n.db <- round((pred.dblt.rate[i])*ncol(cohort.list.scds[[i]]))                         ## Assume doublets based on nuclei isolation protocol performance
  
  # sort prediction, get top n.db cells
  srt.db.score <- sort(cohort.list.scds[[i]]$hybrid_score, index.return = T, decreasing = T)
  cohort.list.scds[[i]]$scds <- "Singlet"
  cohort.list.scds[[i]]$scds[srt.db.score$ix[1:n.db]] <- "Doublet"
  
}

# run UMAP plots
for (i in 1:length(cohort.list.scds)) {
  
  p <- plotReducedDim(cohort.list.scds[[i]], dimred = "UMAP", colour_by = "scds")
  
  pdf(paste(Sys.Date(),"Muscle",names(cohort.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to doubletFinder annotated Seurat object
singlets.annot@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(cohort.list.scds)) {
  
  # for each object compare and move doublet annotations over
  singlets.annot@meta.data[colnames(cohort.list.scds[[i]]), ]$scds_hybrid <- cohort.list.scds[[i]]$scds
  
}

## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& C. Merge and summarize doublet findings

table(singlets.annot@meta.data$DoubletFinder, singlets.annot@meta.data$scds_hybrid)
#          Doublet Singlet
# Doublet     458    3156
# Singlet    3156   59569

# Union (more conservative)
singlets.annot@meta.data$DoubletCall <- ifelse( bitOr(singlets.annot@meta.data$DoubletFinder == "Doublet", singlets.annot@meta.data$scds_hybrid == "Doublet") > 0, 
                                                         "Doublet", "Singlet")
table(singlets.annot@meta.data$DoubletCall)
# Doublet Singlet 
#   6770   59569 

# re-run dimensionality reduction for plotting purposes
singlets.annot <- SCTransform(object = singlets.annot, vars.to.regress =  c("nFeature_RNA", "nCount_RNA", "percent.mito", "Phase"))
singlets.annot <- RunPCA(singlets.annot, npcs = 30)
singlets.annot <- RunUMAP(singlets.annot, dims = 1:30)

pdf(paste0(Sys.Date(),"_Muscle_data_Cortes_lab_UMAP_Singlets_labelled_UNION.pdf"), width = 6, height = 5)
DimPlot(singlets.annot, reduction = "umap", group.by = "DoubletCall", raster = T)
dev.off()

pdf(paste0(Sys.Date(),"_Muscle_data_Cortes_lab_Key_Marker_Gene_expression_plots.pdf"), width = 3.5, height = 3)
FeaturePlot(singlets.annot, features = "Myod1", raster = T)
FeaturePlot(singlets.annot, features = "Acta1", raster = T)
FeaturePlot(singlets.annot, features = "Tgfbi", raster = T)
FeaturePlot(singlets.annot, features = "Pax7", raster = T)
FeaturePlot(singlets.annot, features = "Ptprc"  , raster = T)
FeaturePlot(singlets.annot, features = "Adgre1", raster = T)
dev.off()


# save annotated object
save(singlets.annot, file = paste0(Sys.Date(),"_Muscle_data_Cortes_lab_Seurat_object_with_AnnotatedDoublets.RData"))


### extract/subset only singlets
# save data for singlets df
muscle.singlets   <- subset(singlets.annot, subset = DoubletCall %in% "Singlet")  # only keep singlets
muscle.singlets
# An object of class Seurat 
# 57541 features across 59569 samples within 2 assays 
# Active assay: SCT (23845 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

pdf(paste0(Sys.Date(),"_Muscle_data_Cortes_lab_UMAP_Singlets_ONLY_UNION.pdf"), width = 6, height = 5)
DimPlot(muscle.singlets, reduction = "umap", group.by = "DoubletCall", raster = T)
dev.off()

table(muscle.singlets@meta.data$Group)
#    F_runner F_sedentary      F_TFEB    M_runner M_sedentary      M_TFEB 
#      10189        9678        8733       12399       11367        7203 

# save filtered/annotated object
save(muscle.singlets, file = paste0(Sys.Date(),"_Muscle_data_Cortes_lab_Seurat_object_SINGLETS_ONLY.RData"))
################################################################################################

#######################
sink(file = paste(Sys.Date(),"_Muscle_data_Cortes_lab_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()
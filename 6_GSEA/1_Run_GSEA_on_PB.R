setwd('/Volumes/BB_HQ_4/Collaboration/Cortes_lab/Muscle_snRNAseq/GSEA')
options(stringsAsFactors = F)

#### Packages
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

library(phenoTest)        #
library(qusage)           #

theme_set(theme_bw())   

# 2025-03-07
# Process scRNAseq exercise/TFEb cohorts for differential gene analysis

###############################################################################################
# 1. PhenoTest GO and REACTOME analysis

# Clean workspace and reload necessary objects
load('../DEG_Analysis/2024-12-13_pseudobulk_Muscle_cell_types_RUN_TFEB_Male_DEseq2_objects.RData'     )
load('../DEG_Analysis/2024-12-13_pseudobulk_Muscle_cell_types_RUN_TFEB_Male_VST_data_objects.RData'   )
load('../DEG_Analysis/2024-12-13_pseudobulk_Muscle_cell_types_RUN_TFEB_Female_DEseq2_objects.RData'   )
load('../DEG_Analysis/2024-12-13_pseudobulk_Muscle_cell_types_RUN_TFEB_Female_VST_data_objects.RData' )

ls()
# [1] "deseq.res.list.f.run"  "deseq.res.list.f.tfeb" "deseq.res.list.m.run"  "deseq.res.list.m.tfeb"

# load prepped gene sets
Sym.m5.GO              <- read.gmt("../GeneSets/m5.go.v2024.1.Mm.symbols.gmt")
Sym.m2.reactome        <- read.gmt("../GeneSets/m2.cp.reactome.v2024.1.Mm.symbols.gmt")
Sym.m3.gtrd            <- read.gmt("../GeneSets/m3.gtrd.v2024.1.Mm.symbols.gmt")

# load function
source("Run_PhenoTest_GSEA.R")

##########################################################
#######   GSEA analysis   ++++   Females Running   #######
##########################################################

# Create list object to receive GO ALL results
goall.results.F_run        <- vector(mode = "list", length = length(deseq.res.list.f.run))
names(goall.results.F_run) <- names(deseq.res.list.f.run)

# Create list object to receive REACTOME results
react.results.F_run        <- vector(mode = "list", length = length(deseq.res.list.f.run))
names(react.results.F_run) <- names(deseq.res.list.f.run)

# Create list object to receive GTRD results
gtrd.results.F_run        <- vector(mode = "list", length = length(deseq.res.list.f.run))
names(gtrd.results.F_run) <- names(deseq.res.list.f.run)

# Loop over DEseq2 results
for (i in 1:length(deseq.res.list.f.run)) {
  
  # i. Prepare GeneLists using DEseq2 log2FoldChange to rank genes
  cur.geneList         <- deseq.res.list.f.run[[i]]$log2FoldChange
  names(cur.geneList)  <- rownames(deseq.res.list.f.run[[i]])
  cur.geneList         <- sort(cur.geneList , decreasing = TRUE)
  
  # ii. Run Gene Set Enrichment Analysis 
  goall.results.F_run[[i]]    <- run_enrich(cur.geneList, names(deseq.res.list.f.run)[[i]], "Female_Running", my.fdr = 0.05, my.ontology = Sym.m5.GO        , my.ontology.name = "MsigDB_GOALL"   )
  react.results.F_run[[i]]    <- run_enrich(cur.geneList, names(deseq.res.list.f.run)[[i]], "Female_Running", my.fdr = 0.05, my.ontology = Sym.m2.reactome  , my.ontology.name = "MsigDB_REACTOME")
  gtrd.results.F_run[[i]]     <- run_enrich(cur.geneList, names(deseq.res.list.f.run)[[i]], "Female_Running", my.fdr = 0.05, my.ontology = Sym.m3.gtrd      , my.ontology.name = "MsigDB_GTRD"    )
}


# save R object with all DEseq2 results
save(goall.results.F_run , file = paste0(Sys.Date(),"_pseudobulk_Muscle_F_running_cell_types_GOALL_GSEA.RData"))
save(react.results.F_run , file = paste0(Sys.Date(),"_pseudobulk_Muscle_F_running_cell_types_REACTOME_GSEA.RData"))
save(gtrd.results.F_run  , file = paste0(Sys.Date(),"_pseudobulk_Muscle_F_running_cell_types_GTRD_GSEA.RData"))


##########################################################
#######   GSEA analysis   ++++   Females TFEb      #######
##########################################################

# Create list object to receive GO ALL results
goall.results.F_TFEb        <- vector(mode = "list", length = length(deseq.res.list.f.tfeb))
names(goall.results.F_TFEb) <- names(deseq.res.list.f.tfeb)

# Create list object to receive REACTOME results
react.results.F_TFEb        <- vector(mode = "list", length = length(deseq.res.list.f.tfeb))
names(react.results.F_TFEb) <- names(deseq.res.list.f.tfeb)

# Create list object to receive GTRD results
gtrd.results.F_TFEb        <- vector(mode = "list", length = length(deseq.res.list.f.tfeb))
names(gtrd.results.F_TFEb) <- names(deseq.res.list.f.tfeb)

# Loop over DEseq2 results
for (i in 1:length(deseq.res.list.f.tfeb)) {
  
  # i. Prepare GeneLists using DEseq2 log2FoldChange to rank genes
  cur.geneList         <- deseq.res.list.f.tfeb[[i]]$log2FoldChange
  names(cur.geneList)  <- rownames(deseq.res.list.f.tfeb[[i]])
  cur.geneList         <- sort(cur.geneList , decreasing = TRUE)
  
  # ii. Run Gene Set Enrichment Analysis 
  goall.results.F_TFEb[[i]]    <- run_enrich(cur.geneList, names(deseq.res.list.f.tfeb)[[i]], "Female_TFEb", my.fdr = 0.05, my.ontology = Sym.m5.GO        , my.ontology.name = "MsigDB_GOALL"   )
  react.results.F_TFEb[[i]]    <- run_enrich(cur.geneList, names(deseq.res.list.f.tfeb)[[i]], "Female_TFEb", my.fdr = 0.05, my.ontology = Sym.m2.reactome  , my.ontology.name = "MsigDB_REACTOME")
  gtrd.results.F_TFEb[[i]]     <- run_enrich(cur.geneList, names(deseq.res.list.f.tfeb)[[i]], "Female_TFEb", my.fdr = 0.05, my.ontology = Sym.m3.gtrd      , my.ontology.name = "MsigDB_GTRD"    )
}


# save R object with all DEseq2 results
save(goall.results.F_TFEb , file = paste0(Sys.Date(),"_pseudobulk_Muscle_F_TFEb_cell_types_GOALL_GSEA.RData"))
save(react.results.F_TFEb , file = paste0(Sys.Date(),"_pseudobulk_Muscle_F_TFEb_cell_types_REACTOME_GSEA.RData"))
save(gtrd.results.F_TFEb  , file = paste0(Sys.Date(),"_pseudobulk_Muscle_F_TFEb_cell_types_GTRD_GSEA.RData"))


##########################################################
#######   GSEA analysis   ++++   Males Running     #######
##########################################################

# Create list object to receive GO ALL results
goall.results.M_run        <- vector(mode = "list", length = length(deseq.res.list.m.run))
names(goall.results.M_run) <- names(deseq.res.list.m.run)

# Create list object to receive REACTOME results
react.results.M_run        <- vector(mode = "list", length = length(deseq.res.list.m.run))
names(react.results.M_run) <- names(deseq.res.list.m.run)

# Create list object to receive GTRD results
gtrd.results.M_run        <- vector(mode = "list", length = length(deseq.res.list.m.run))
names(gtrd.results.M_run) <- names(deseq.res.list.m.run)

# Loop over DEseq2 results
for (i in 1:length(deseq.res.list.m.run)) {
  
  # i. Prepare GeneLists using DEseq2 log2FoldChange to rank genes
  cur.geneList         <- deseq.res.list.m.run[[i]]$log2FoldChange
  names(cur.geneList)  <- rownames(deseq.res.list.m.run[[i]])
  cur.geneList         <- sort(cur.geneList , decreasing = TRUE)
  
  # ii. Run Gene Set Enrichment Analysis 
  goall.results.M_run[[i]]    <- run_enrich(cur.geneList, names(deseq.res.list.m.run)[[i]], "Male_Running", my.fdr = 0.05, my.ontology = Sym.m5.GO        , my.ontology.name = "MsigDB_GOALL"   )
  react.results.M_run[[i]]    <- run_enrich(cur.geneList, names(deseq.res.list.m.run)[[i]], "Male_Running", my.fdr = 0.05, my.ontology = Sym.m2.reactome  , my.ontology.name = "MsigDB_REACTOME")
  gtrd.results.M_run[[i]]     <- run_enrich(cur.geneList, names(deseq.res.list.m.run)[[i]], "Male_Running", my.fdr = 0.05, my.ontology = Sym.m3.gtrd      , my.ontology.name = "MsigDB_GTRD"    )
}


# save R object with all DEseq2 results
save(goall.results.M_run , file = paste0(Sys.Date(),"_pseudobulk_Muscle_M_running_cell_types_GOALL_GSEA.RData"))
save(react.results.M_run , file = paste0(Sys.Date(),"_pseudobulk_Muscle_M_running_cell_types_REACTOME_GSEA.RData"))
save(gtrd.results.M_run  , file = paste0(Sys.Date(),"_pseudobulk_Muscle_M_running_cell_types_GTRD_GSEA.RData"))


##########################################################
#######   GSEA analysis   ++++   Males TFEb        #######
##########################################################


# Create list object to receive GO ALL results
goall.results.M_TFEb        <- vector(mode = "list", length = length(deseq.res.list.m.tfeb))
names(goall.results.M_TFEb) <- names(deseq.res.list.m.tfeb)

# Create list object to receive REACTOME results
react.results.M_TFEb        <- vector(mode = "list", length = length(deseq.res.list.m.tfeb))
names(react.results.M_TFEb) <- names(deseq.res.list.m.tfeb)

# Create list object to receive GTRD results
gtrd.results.M_TFEb        <- vector(mode = "list", length = length(deseq.res.list.m.tfeb))
names(gtrd.results.M_TFEb) <- names(deseq.res.list.m.tfeb)

# Loop over DEseq2 results
for (i in 1:length(deseq.res.list.m.tfeb)) {
  
  # i. Prepare GeneLists using DEseq2 log2FoldChange to rank genes
  cur.geneList         <- deseq.res.list.m.tfeb[[i]]$log2FoldChange
  names(cur.geneList)  <- rownames(deseq.res.list.m.tfeb[[i]])
  cur.geneList         <- sort(cur.geneList , decreasing = TRUE)
  
  # ii. Run Gene Set Enrichment Analysis 
  goall.results.M_TFEb[[i]]    <- run_enrich(cur.geneList, names(deseq.res.list.m.tfeb)[[i]], "Male_TFEb", my.fdr = 0.05, my.ontology = Sym.m5.GO        , my.ontology.name = "MsigDB_GOALL"   )
  react.results.M_TFEb[[i]]    <- run_enrich(cur.geneList, names(deseq.res.list.m.tfeb)[[i]], "Male_TFEb", my.fdr = 0.05, my.ontology = Sym.m2.reactome  , my.ontology.name = "MsigDB_REACTOME")
  gtrd.results.M_TFEb[[i]]     <- run_enrich(cur.geneList, names(deseq.res.list.m.tfeb)[[i]], "Male_TFEb", my.fdr = 0.05, my.ontology = Sym.m3.gtrd      , my.ontology.name = "MsigDB_GTRD"    )
}


# save R object with all DEseq2 results
save(goall.results.M_TFEb , file = paste0(Sys.Date(),"_pseudobulk_Muscle_M_TFEb_cell_types_GOALL_GSEA.RData"))
save(react.results.M_TFEb , file = paste0(Sys.Date(),"_pseudobulk_Muscle_M_TFEb_cell_types_REACTOME_GSEA.RData"))
save(gtrd.results.M_TFEb  , file = paste0(Sys.Date(),"_pseudobulk_Muscle_M_TFEb_cell_types_GTRD_GSEA.RData"))
###############################################################################################

###############################################################################################
# 2. Plot PhenoTest GO, REACTOME and GTRD analysis

source('Plot_bubble_chart_function.R')

for(i in 1:length(gtrd.results.F_run)) {
  if(nrow(gtrd.results.F_run[[i]])>0) { # only if something passed significance
    get_age_bubble_plot(names(gtrd.results.F_run)[i], gtrd.results.F_run[[i]], max.path.plot = 10, my.onto = "F_RUN_GTRD")
  }
  if(nrow(goall.results.F_run[[i]])>0) { # only if something passed significance
    get_age_bubble_plot(names(goall.results.F_run)[i], goall.results.F_run[[i]], max.path.plot = 10, my.onto = "F_RUN_GO_ALL")
  }
  if(nrow(react.results.F_run[[i]])>0) { # only if something passed significance
    get_age_bubble_plot(names(react.results.F_run)[i], react.results.F_run[[i]], max.path.plot = 10, my.onto = "F_RUN_REACTOME")
  }
}


for(i in 1:length(gtrd.results.F_TFEb)) {
  if(nrow(gtrd.results.F_TFEb[[i]])>0) { # only if something passed significance
    get_age_bubble_plot(names(gtrd.results.F_TFEb)[i], gtrd.results.F_TFEb[[i]], max.path.plot = 10, my.onto = "F_TFEb_GTRD")
  }
  if(nrow(goall.results.F_TFEb[[i]])>0) { # only if something passed significance
    get_age_bubble_plot(names(goall.results.F_TFEb)[i], goall.results.F_TFEb[[i]], max.path.plot = 10, my.onto = "F_TFEb_GO_ALL")
  }
  if(nrow(react.results.F_TFEb[[i]])>0) { # only if something passed significance
    get_age_bubble_plot(names(react.results.F_TFEb)[i], react.results.F_TFEb[[i]], max.path.plot = 10, my.onto = "F_TFEb_REACTOME")
  }
}

for(i in 1:length(gtrd.results.M_run)) {
  if(nrow(gtrd.results.M_run[[i]])>0) { # only if something passed significance
    get_age_bubble_plot(names(gtrd.results.M_run)[i], gtrd.results.M_run[[i]], max.path.plot = 10, my.onto = "M_RUN_GTRD")
  }
  if(nrow(goall.results.M_run[[i]])>0) { # only if something passed significance
    get_age_bubble_plot(names(goall.results.M_run)[i], goall.results.M_run[[i]], max.path.plot = 10, my.onto = "M_RUN_GO_ALL")
  }
  if(nrow(react.results.M_run[[i]])>0) { # only if something passed significance
    get_age_bubble_plot(names(react.results.M_run)[i], react.results.M_run[[i]], max.path.plot = 10, my.onto = "M_RUN_REACTOME")
  }
}


for(i in 1:length(gtrd.results.M_TFEb)) {
  if(nrow(gtrd.results.M_TFEb[[i]])>0) { # only if something passed significance
    get_age_bubble_plot(names(gtrd.results.M_TFEb)[i], gtrd.results.M_TFEb[[i]], max.path.plot = 10, my.onto = "M_TFEb_GTRD")
  }
  if(nrow(goall.results.M_TFEb[[i]])>0) { # only if something passed significance
    get_age_bubble_plot(names(goall.results.M_TFEb)[i], goall.results.M_TFEb[[i]], max.path.plot = 10, my.onto = "M_TFEb_GO_ALL")
  }
  if(nrow(react.results.M_TFEb[[i]])>0) { # only if something passed significance
    get_age_bubble_plot(names(react.results.M_TFEb)[i], react.results.M_TFEb[[i]], max.path.plot = 10, my.onto = "M_TFEb_REACTOME")
  }
}


# #######   Summary plot across cell types   #######
source('Plot_cross_cell_type_summary_bubble_Plot_v2.R')

get_age_summary_plot(goall.results.F_run , my.onto = "F_RUN_GO_ALL"        , n.thresh = 4, max.path.plot = 10)
get_age_summary_plot(react.results.F_run , my.onto = "F_RUN_REACTOME"      , n.thresh = 4, max.path.plot = 10)

get_age_summary_plot(goall.results.F_TFEb , my.onto = "F_TFEb_GO_ALL"        , n.thresh = 4, max.path.plot = 10)
get_age_summary_plot(react.results.F_TFEb , my.onto = "F_TFEb_REACTOME"      , n.thresh = 4, max.path.plot = 10)

get_age_summary_plot(goall.results.M_run , my.onto = "M_RUN_GO_ALL"        , n.thresh = 4, max.path.plot = 10)
get_age_summary_plot(react.results.M_run , my.onto = "M_RUN_REACTOME"      , n.thresh = 4, max.path.plot = 10)

get_age_summary_plot(goall.results.M_TFEb , my.onto = "M_TFEb_GO_ALL"        , n.thresh = 4, max.path.plot = 10)
get_age_summary_plot(react.results.M_TFEb , my.onto = "M_TFEb_REACTOME"      , n.thresh = 4, max.path.plot = 10)
#################################################################################################################################################################



#######################
sink(file = paste(Sys.Date(),"_MuscatDEseq2_PB_DESeq2_GSEA_snRNAseq_Muscle_Atlas_session_Info.txt", sep =""))
sessionInfo()
sink()
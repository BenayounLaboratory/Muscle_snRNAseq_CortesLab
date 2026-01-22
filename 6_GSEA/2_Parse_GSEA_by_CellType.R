setwd('/Volumes/BB_HQ_4/Collaboration/Cortes_lab/Muscle_snRNAseq/GSEA')
options(stringsAsFactors = F)

#### Packages
library(ggplot2)          # 
library(scales)           # 
library("bitops")         # 
library(Vennerable)       # 
library(data.table)       #


theme_set(theme_bw())   

# 2025-03-11
# Process scRNAseq exercise/TFEb cohorts for differential gene analysis

###############################################################################################
# 1. Load PhenoTest GO and REACTOME results

load('2025-03-07_pseudobulk_Muscle_F_running_cell_types_REACTOME_GSEA.RData' )
load('2025-03-07_pseudobulk_Muscle_F_running_cell_types_GOALL_GSEA.RData' )
load('2025-03-07_pseudobulk_Muscle_F_TFEb_cell_types_GOALL_GSEA.RData' )
load('2025-03-07_pseudobulk_Muscle_F_TFEb_cell_types_REACTOME_GSEA.RData' )
load('2025-03-07_pseudobulk_Muscle_M_running_cell_types_GOALL_GSEA.RData' )
load('2025-03-07_pseudobulk_Muscle_M_running_cell_types_REACTOME_GSEA.RData' )
load('2025-03-07_pseudobulk_Muscle_M_TFEb_cell_types_GOALL_GSEA.RData' )
load('2025-03-07_pseudobulk_Muscle_M_TFEb_cell_types_REACTOME_GSEA.RData' )

ls()
# "goall.results.F_run"  "goall.results.F_TFEb" "goall.results.M_run"  "goall.results.M_TFEb" 
# "react.results.F_run" "react.results.F_TFEb" "react.results.M_run"  "react.results.M_TFEb"
###############################################################################################

###############################################################################################
# 2. merge results for each cell type

celltypes <- gsub(" ", "", names(goall.results.F_run))

merged.go.list    <- vector(mode = "list", length = length(celltypes))
merged.react.list <- vector(mode = "list", length = length(celltypes))

names(merged.go.list) <- celltypes
names(merged.react.list) <- celltypes

for (i in 1:length(celltypes)) {
  
  ### GO ALL
  go.tmp.f  <- merge(goall.results.F_run[[i]], goall.results.F_TFEb[[i]], by = "row.names", suffixes = c(".run",".tfeb"), all = T)
  go.tmp.m  <- merge(goall.results.M_run[[i]], goall.results.M_TFEb[[i]], by = "row.names", suffixes = c(".run",".tfeb"), all = T)
  merged.go.list[[i]] <- merge(go.tmp.f, go.tmp.m, by = "Row.names", suffixes = c(".F",".M"), all = T)
  colnames(merged.go.list[[i]])[1] <- "GeneSet"
  
  ## REACTOME
  react.tmp.f <- merge(react.results.F_run[[i]], react.results.F_TFEb[[i]], by = "row.names", suffixes = c(".run",".tfeb"), all = T)
  react.tmp.m <- merge(react.results.M_run[[i]], react.results.M_TFEb[[i]], by = "row.names", suffixes = c(".run",".tfeb"), all = T)
  merged.react.list[[i]] <- merge(react.tmp.f, react.tmp.m, by = "Row.names", suffixes = c(".F",".M"), all = T)
  colnames(merged.react.list[[i]])[1] <- "GeneSet"
}

save(merged.go.list,merged.react.list, file = paste0(Sys.Date(),"_merged_enrichments_GO_REACTOME.RData") )
###############################################################################################

###############################################################################################
# 3. plot top recurrent for each cell type

load('2025-03-11_merged_enrichments_GO_REACTOME.RData')

celltypes <- names(merged.go.list)


for(i in 1:length(celltypes)) {
  # get those observed in all conditions
  tmp.mat <- merged.go.list[[i]]
  
  all.sig  <- tmp.mat[rowSums(!is.na(tmp.mat[,grep("fdr."   , colnames(tmp.mat))])>0) == 4,]
  
  all.sig$comb_fdr <- as.numeric(apply(all.sig[,grep("fdr."   , colnames(all.sig))]+ 1e-20,1,"prod"))
  
  ## get top 10
  sort.ix <- sort(all.sig$comb_fdr, decreasing = F, index.return = T) # top 10 most sig
  top.sig <- all.sig[sort.ix$ix[1:10],]
  
  # tabulate for ggplot
  all.sig.tab <- data.frame(rbindlist(list(top.sig[,c(1,grep("run.F", colnames(top.sig)))] ,
                                           top.sig[,c(1,grep("tfeb.F",colnames(top.sig)))],
                                           top.sig[,c(1,grep("run.M", colnames(top.sig)))] ,
                                           top.sig[,c(1,grep("tfeb.M",colnames(top.sig)))] )  ))
  all.sig.tab$Group <- factor(c(rep("F_run",nrow(top.sig)),
                                rep("F_TFEb",nrow(top.sig)),
                                rep("M_run",nrow(top.sig)),
                                rep("M_TFEb",nrow(top.sig))))
  
  # clean
  all.sig.tab$minlog10fdr  <- -log10(all.sig.tab$fdr.run.F + 1e-20)
  colnames(all.sig.tab) <- c("Path", "n", "es", "NES","pval","fdr","Group","minlog10fdr")
  
  all.sig.tab$Path <- tolower(gsub("_"," ", all.sig.tab$Path))
  all.sig.tab$Path <- factor(all.sig.tab$Path)
  
  # Up/Down color scale
  my.max <-  3.5 # add a value in case there are not terms biased in one direction
  my.min <- -3.5 # add a value in case there are not terms biased in one direction
  my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
  my.scaled <- rescale(my.values, to = c(0, 1))
  
  my.color.vector.age <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")
  
  my.plot <- ggplot(all.sig.tab, aes(x=Group, y=Path, colour=NES, size=minlog10fdr))+ theme_bw()+ geom_point(shape = 16)
  my.plot <- my.plot + ggtitle(paste0(celltypes[i], " GO gsea Analysis (Top 10)")) + labs(x = "-log10(FDR)", y = "")
  my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.age, 
                                              na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min, my.max) )
  my.plot <- my.plot + scale_size_area(limits = c(1,20)) + scale_y_discrete(labels = wrap_format(40))
  my.plot
  
  my.pdfname <- paste(Sys.Date(),"GSEA_BALLOON_plot",celltypes[i], "top", nrow(top.sig),"significant_GO_terms.pdf", sep="_")
  
  pdf(my.pdfname, onefile=T, height = 5, width= 6)
  print(my.plot)
  dev.off()  
  
  
}

####

for(i in 1:length(celltypes)) {
  # get those observed in all conditions
  tmp.mat <- merged.react.list[[i]]
  
  all.sig  <- tmp.mat[rowSums(!is.na(tmp.mat[,grep("fdr."   , colnames(tmp.mat))])>0) == 4,]
  
  all.sig$comb_fdr <- as.numeric(apply(all.sig[,grep("fdr."   , colnames(all.sig))]+ 1e-20,1,"prod"))
  
  ## get top 10
  sort.ix <- sort(all.sig$comb_fdr, decreasing = F, index.return = T) # top 10 most sig
  top.sig <- all.sig[sort.ix$ix[1:10],]
  
  # tabulate for ggplot
  all.sig.tab <- data.frame(rbindlist(list(top.sig[,c(1,grep("run.F", colnames(top.sig)))] ,
                                           top.sig[,c(1,grep("tfeb.F",colnames(top.sig)))],
                                           top.sig[,c(1,grep("run.M", colnames(top.sig)))] ,
                                           top.sig[,c(1,grep("tfeb.M",colnames(top.sig)))] )  ))
  all.sig.tab$Group <- factor(c(rep("F_run",nrow(top.sig)),
                                rep("F_TFEb",nrow(top.sig)),
                                rep("M_run",nrow(top.sig)),
                                rep("M_TFEb",nrow(top.sig))))
  
  # clean
  all.sig.tab$minlog10fdr  <- -log10(all.sig.tab$fdr.run.F + 1e-20)
  colnames(all.sig.tab) <- c("Path", "n", "es", "NES","pval","fdr","Group","minlog10fdr")
  
  all.sig.tab$Path <- tolower(gsub("_"," ", all.sig.tab$Path))
  all.sig.tab$Path <- factor(all.sig.tab$Path)
  
  # Up/Down color scale
  my.max <-  3.5 # add a value in case there are not terms biased in one direction
  my.min <- -3.5 # add a value in case there are not terms biased in one direction
  my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
  my.scaled <- rescale(my.values, to = c(0, 1))
  
  my.color.vector.age <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")
  
  my.plot <- ggplot(all.sig.tab, aes(x=Group, y=Path, colour=NES, size=minlog10fdr))+ theme_bw()+ geom_point(shape = 16)
  my.plot <- my.plot + ggtitle(paste0(celltypes[i], " GO gsea Analysis (Top 10)")) + labs(x = "-log10(FDR)", y = "")
  my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.age, 
                                              na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min, my.max) )
  my.plot <- my.plot + scale_size_area(limits = c(1,20)) + scale_y_discrete(labels = wrap_format(40))
  my.plot
  
  my.pdfname <- paste(Sys.Date(),"GSEA_BALLOON_plot",celltypes[i], "top", nrow(top.sig),"significant_REACTOME_terms.pdf", sep="_")
  
  pdf(my.pdfname, onefile=T, height = 5, width= 6)
  print(my.plot)
  dev.off()  
  
  
}
# #################################################################################################################################################################
# 


#######################
sink(file = paste(Sys.Date(),"_MuscatDEseq2_PB_DESeq2_GSEA_snRNAseq_Muscle_Atlas_session_Info.txt", sep =""))
sessionInfo()
sink()
setwd('/Volumes/BB_HQ_4/Collaboration/Cortes_lab/Muscle_snRNAseq/DEG_analysis/')
options(stringsAsFactors = F)

#### Packages
library(ggplot2)          # 
library(scales)           # 
library("bitops")         # 
library(Vennerable)       # 
library(data.table)       #

library(ComplexHeatmap)   #
library(circlize)         #


theme_set(theme_bw())   


# 2025-03-13
# get jitter plots

###############################################################################################
# 0. load DESeq2 objects

load('2024-12-13_pseudobulk_Muscle_cell_types_RUN_TFEB_Male_DEseq2_objects.RData')
load('2024-12-13_pseudobulk_Muscle_cell_types_RUN_TFEB_Female_DEseq2_objects.RData')

ls()
# [1] "deseq.res.list.f.run"  "deseq.res.list.f.tfeb" "deseq.res.list.m.run"  "deseq.res.list.m.tfeb"
###############################################################################################


###############################################################################################
# 1.  Make jitter plot of DE genes

make_jitter <- function(deseq2.list, condition, neg_col =   "#333399", pos_col = "#CC3333", my.fdr = 0.05) {
  ######## 
  ## Order by pvalue:
  comp.res <- lapply(deseq2.list,function(x) {x[order(x$padj),]})
  n        <- sapply(comp.res, nrow)
  names(n) <- names(comp.res)
  
  cols <- list()
  xlab <- character(length = length(comp.res))
  for(i in seq(along = comp.res)){
    cols[[i]]      <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), n[i]) # grey60
    ind.sig.i      <- comp.res[[i]]$padj < my.fdr
    ind.sig.i.up   <- bitAnd(comp.res[[i]]$padj < my.fdr, comp.res[[i]]$log2FoldChange >0)>0
    ind.sig.i.down <- bitAnd(comp.res[[i]]$padj < my.fdr, comp.res[[i]]$log2FoldChange <0)>0
    
    cols[[i]][ind.sig.i.up]   <- pos_col
    cols[[i]][ind.sig.i.down] <- neg_col
    xlab[i] <- paste(names(comp.res)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
  }
  names(cols) <- names(comp.res)
  
  pdf(paste0(Sys.Date(),condition,"_stripplot_DESeq2_with_reg_colors_FDR", my.fdr, ".pdf"), width = 6, height = 5.5)
  par(mar = c(3.1, 4.1, 1, 1))
  par(oma = c(6, 2, 1, 1))
  plot(x = 1,
       y = 1,
       type = "n",
       xlim = c(0.5, 6.5),
       ylim = c(-10, 10),
       axes = FALSE,
       xlab = "",
       ylab = "Log2 fold change",
       main = condition
  )
  abline(h = 0)
  abline(h = seq(-10, 10, by = 2)[-6],
         lty = "dotted",
         col = "grey")
  for(i in 1:length(comp.res)){
    set.seed(1234)
    points(x = jitter(rep(i, nrow(comp.res[[i]])), amount = 0.2),
           y = rev(comp.res[[i]]$log2FoldChange),
           pch = 16,
           col = rev(cols[[i]]),
           bg = rev(cols[[i]]),
           cex = 0.75)
  }
  axis(1,
       at = 1:6,
       tick = FALSE,
       las = 2,
       lwd = 0,
       labels = xlab,
       cex.axis = 0.7)
  axis(2,
       las = 1,
       at = seq(-10, 10, by = 2))
  box()
  dev.off()
  
  png(paste0(Sys.Date(),condition,"_stripplot_DESeq2_with_reg_colors_FDR", my.fdr, ".png"), width = 2400, height = 2200, res = 300, units = "px")
  par(mar = c(3.1, 4.1, 1, 1))
  par(oma = c(6, 2, 1, 1))
  plot(x = 1,
       y = 1,
       type = "n",
       xlim = c(0.5, 6.5),
       ylim = c(-10, 10),
       axes = FALSE,
       xlab = "",
       ylab = "Log2 fold change",
       main = condition
  )
  abline(h = 0)
  abline(h = seq(-10, 10, by = 2)[-6],
         lty = "dotted",
         col = "grey")
  for(i in 1:length(comp.res)){
    set.seed(1234)
    points(x = jitter(rep(i, nrow(comp.res[[i]])), amount = 0.2),
           y = rev(comp.res[[i]]$log2FoldChange),
           pch = 16,
           col = rev(cols[[i]]),
           bg = rev(cols[[i]]),
           cex = 0.75)
  }
  axis(1,
       at = 1:6,
       tick = FALSE,
       las = 2,
       lwd = 0,
       labels = xlab,
       cex.axis = 0.7)
  axis(2,
       las = 1,
       at = seq(-10, 10, by = 2))
  box()
  dev.off()
  
  
  pdf(paste0(Sys.Date(),condition,"_stripplot_DESeq2_with_reg_colors_FDR", my.fdr, "NO_DOTS.pdf"), width = 6, height = 5.5)
  par(mar = c(3.1, 4.1, 1, 1))
  par(oma = c(6, 2, 1, 1))
  plot(x = 1,
       y = 1,
       type = "n",
       xlim = c(0.5, 6.5),
       ylim = c(-10, 10),
       axes = FALSE,
       xlab = "",
       ylab = "Log2 fold change",
       main = condition
  )
  abline(h = 0)
  abline(h = seq(-10, 10, by = 2)[-6],
         lty = "dotted",
         col = "grey")
  # for(i in 1:length(comp.res)){
  #   set.seed(1234)
  #   points(x = jitter(rep(i, nrow(comp.res[[i]])), amount = 0.2),
  #          y = rev(comp.res[[i]]$log2FoldChange),
  #          pch = 16,
  #          col = rev(cols[[i]]),
  #          bg = rev(cols[[i]]),
  #          cex = 0.75)
  # }
  axis(1,
       at = 1:6,
       tick = FALSE,
       las = 2,
       lwd = 0,
       labels = xlab,
       cex.axis = 0.7)
  axis(2,
       las = 1,
       at = seq(-10, 10, by = 2))
  box()
  dev.off()
}


# FDR5
make_jitter(deseq.res.list.f.run , "F_running", neg_col =   "deeppink"   , pos_col = "lightpink2"     , my.fdr = 0.05)
make_jitter(deseq.res.list.f.tfeb, "F_TFEb"   , neg_col =   "deeppink"   , pos_col = "mediumpurple4"  , my.fdr = 0.05)
make_jitter(deseq.res.list.m.run , "M_running", neg_col =   "deepskyblue", pos_col = "mediumturquoise", my.fdr = 0.05)
make_jitter(deseq.res.list.m.tfeb, "M_TFEb"   , neg_col =   "deepskyblue", pos_col = "royalblue3"     , my.fdr = 0.05)

# FDR10
make_jitter(deseq.res.list.f.run , "F_running", neg_col =   "deeppink"   , pos_col = "lightpink2"     , my.fdr = 0.1)
make_jitter(deseq.res.list.f.tfeb, "F_TFEb"   , neg_col =   "deeppink"   , pos_col = "mediumpurple4"  , my.fdr = 0.1)
make_jitter(deseq.res.list.m.run , "M_running", neg_col =   "deepskyblue", pos_col = "mediumturquoise", my.fdr = 0.1)
make_jitter(deseq.res.list.m.tfeb, "M_TFEb"   , neg_col =   "deepskyblue", pos_col = "royalblue3"     , my.fdr = 0.1)


##########################################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_MuscatDEseq2_PB_DESeq2_scRNAseq_Muscle_session_Info.txt", sep =""))
sessionInfo()
sink()
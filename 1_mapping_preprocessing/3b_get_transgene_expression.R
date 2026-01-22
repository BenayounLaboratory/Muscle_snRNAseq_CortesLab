setwd('/Volumes/BB_HQ_4/Collaboration/Cortes_lab/Muscle_snRNAseq/Trasngene_Counting')

library(beeswarm)

### read kallisto results,
my.ct.files <- c('Female_6mo_sed_1.kallisto_res/abundance.tsv',
                 'Female_6mo_sed_2.kallisto_res/abundance.tsv',
                 'Female_6mo_runner_1.kallisto_res/abundance.tsv',
                 'Female_6mo_runner_2.kallisto_res/abundance.tsv',
                 'Female_6mo_TFEB_1.kallisto_res/abundance.tsv',
                 'Female_6mo_TFEB_2.kallisto_res/abundance.tsv',
                 'Male_6mo_sed_1.kallisto_res/abundance.tsv',
                 'Male_6mo_sed_2.kallisto_res/abundance.tsv',
                 'Male_6mo_runner_1.kallisto_res/abundance.tsv',
                 'Male_6mo_runner_2.kallisto_res/abundance.tsv',
                 'Male_6mo_TFEB_1.kallisto_res/abundance.tsv',
                 'Male_6mo_TFEB_2.kallisto_res/abundance.tsv')

# get sample name
names(my.ct.files) <- unlist(lapply(strsplit(my.ct.files, ".kallisto"), '[[',1))

# read 1st file
first.mat <- read.table(my.ct.files[1], header = T)

# start result matrix
tpm.mat <- first.mat[,c(1,5)]
colnames(tpm.mat)[2] <- names(my.ct.files) [1]

# parse into result matrix
for (i in 2:length(my.ct.files)) {
  my.name <- names(my.ct.files)[i]
  my.tmp  <- data.frame(read.table(my.ct.files[i], header = T)[,5])
  colnames(my.tmp) <- my.name
  
  tpm.mat <- cbind(tpm.mat,my.tmp)
  
}

tpm.mat <- data.frame(tpm.mat)

tpm.mat[2,-1]

barplot(as.numeric(tpm.mat[2,-1]))


my.res <- list("F_sed"   = as.numeric(tpm.mat[2,2:3]),
               "F_run"   = as.numeric(tpm.mat[2,4:5]),
               "F_TFEb"  = as.numeric(tpm.mat[2,6:7]),
               "M_sed"   = as.numeric(tpm.mat[2,8:9]),
               "M_run"   = as.numeric(tpm.mat[2,10:11]),
               "M_TFEb"  = as.numeric(tpm.mat[2,12:13]))

pdf(paste0(Sys.Date(),"_human_TFEb_expression_per_library.pdf"), height = 5, width = 3.5)
beeswarm(my.res, pch = 16, 
         pwcol = c("deeppink","deeppink",
                   "lightpink2","lightpink2",
                   "mediumpurple4","mediumpurple4",
                   "deepskyblue","deepskyblue",
                   "mediumturquoise","mediumturquoise",
                   "royalblue3","royalblue3"),
         ylab = "Transcript expression per 10x library (Kallisto TPM)",
         ylim = c(0,60),
         cex = 1.5,
         main = "Human TFEb expression",
         las = 3)
abline(v = 3.5, col = "grey", lty = "dashed")
abline(h = 0  , col = "red", lty = "dashed")
dev.off()


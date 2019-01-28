# Additional individual analyses
rm(list=ls())
setwd("/home/boyangzhao/NMDcyt/")

#libraries
if(!require("pacman")) install.packages("pacman")
pacman::p_load(caret, car, MASS, ROCR, randomForest, rfPermute)
pacman::p_load(gplots, ggplot2, ggrepel, GGally, RColorBrewer, ggpubr, factoextra)
pacman::p_load(Hmisc)
pacman::p_load(plyr, dplyr)
pacman::p_load(future)
source('src/utils.R')
source('src/utils_surv.R')

#------- parameters --------
dir.out <- './outputs/18.0604/' #output directory
dir.db <- './datasets/' #location of the datasets
dirout <- file.path(dir.out, '/summaries_aggr')

load(file.path(dir.out, 'pancancer.rda'))

#there is one outlier with value greater than 15000, removed
df.merged <- df.merged[df.merged$silent < 15000,]

#------- scatter plots --------
ggplot(df.merged, aes(silent, missense)) + geom_point(size=3, alpha=0.7) + geom_smooth(method='lm', level=0.95)
ggsave(file.path(dirout, 'scatter_silent_missense.pdf'))

ggplot(df.merged, aes(silent, nonsense)) + geom_point(size=3, alpha=0.7) + geom_smooth(method='lm', level=0.95)
ggsave(file.path(dirout, 'scatter_silent_nonsense.pdf'))

ggplot(df.merged, aes(silent, frameshift)) + geom_point(size=3, alpha=0.7)
ggsave(file.path(dirout, 'scatter_silent_frameshift.pdf'))

#------- heatmap --------
varnames <- colnames(df.merged)
selnames <- c('silent', 'missense', 'nonstop', 'nonsense', 'frameshift', 'total', varnames[grepl('^nmd', varnames)],
              'GZMA', 'PRF1', 'cyt', 'cyt.bool')
df <- df.merged[, selnames]
c <- rcorr(sapply(df, as.numeric), type="spearman")

pdf(file.path(dirout, paste0('corr_matrix_mutNMDonly.pdf')), width=12, height=12)
heatmap.2(as.matrix(c$r), 
          col=colorRampPalette(brewer.pal(9, "RdBu"))(100), scale="none", 
          key.xlab='Spearman corr', key.title=F,
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5, cexCol=0.5,
          margins = c(8,8))
dev.off()

#------- forest plot --------
library(forestplot)
dirout <- file.path(dir.out, '/summaries_aggr')

#SKCM uni.cox.withctrl2/univariate_cox_withctrl_TMBPDL1cat0
#for nmdptc_med
ds.forest <- 
  structure(list(
    mean  = c(NA, 1.2, 1.8, 1.4, 1.2, 0.4, 1.4, 2.1, 3.4, 0.6, 0.5, 0.4), 
    lower = c(NA, 1.0, 1.2, 1.0, 0.4, 0.1, 0.5, 0.7, 1.0, 0.5, 0.4, 0.2),
    upper = c(NA, 1.4, 2.5, 2.0, 3.4, 1.8, 4.0, 6.0, 11.7, 0.9, 0.8, 0.6)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -11L), 
    class = "data.frame")

hr <- paste0(ds.forest$mean, " (", ds.forest$lower, "-", ds.forest$upper, ")")
hr[1] <- "Hazard ratio (95% CI)"

tabletext<-cbind(
  c("Covariate", 
    "nmdptc_med", "Age (>=65 vs <65)", "Gender (male vs female)",
    "TNM stage I", "TNM stage I/II NOS", "TNM stage II", "TNM stage III", "TNM stage IV", 
    "TMB (>median vs <= median)", "PDL1 med", "PDL1 high"),
  hr,
  c("p-value",
    "0.0364", "0.00156",  "0.0482", 
    "0.797", "0.217", "0.522",  "0.165", "0.0534",
    "0.011", "0.00154", "0.000148") )

pdf(file.path(dirout, paste0('forest_ind_SKCM_multi_nmdptc_med.pdf')), width=12, height=12)
forestplot(tabletext, 
           ds.forest,
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,11)),
           xlog=TRUE, 
           col=fpColors(box="black",line="black", summary="black"))
dev.off()

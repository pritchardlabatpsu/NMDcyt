# Analyze models - aggregate summaries
#  - generate combined ROC curves
# Model results are saved in as .rda - to be used as input here
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
pval.ptc.excludeNoMut <- TRUE #exclude non mut patients, for plotting/mann-whitney test on NMD against discretized CYT
param.parallel <- TRUE #use futures, for parallel computing

cancertypes <- c('BLCA','BRCA','CESC','COADREAD','GBMLGG','HNSC','KIPAN',
                 'LUAD','LUSC','OV','PRAD','SKCM','STAD','THCA','UCEC')

#------- combine ROC curves
fnames <- list('bin_rf_mut_fl.rda', 'bin_rf_NMD_fl.rda', 'bin_rf_mut+NMD_fl.rda')
mod.legnames.default <- c('Mut', 'NMD', 'Mut+NMD')

if(!file.exists(file.path(dir.out, 'summaries'))) dir.create(file.path(dir.out, 'summaries'))

for(cancertype in cancertypes){
  dirout <- file.path(dir.out, paste0('_',cancertype))
  
  a.name <- 'train.val_set'
  mod.legnames <- mod.legnames.default
  pdf(file.path(dir.out, 'summaries', paste0('ROC_', cancertype, '_', a.name, '.pdf')))
  n = 1
  for(fname in fnames){
    load(file.path(dirout, fname))
    
    mod <- modelobj.save[[1]]
    mod.cv <- modelobj.save[[2]]
    yname <- modelobj.save[[3]]
    mod.data.tval <- modelobj.save[[4]]
    mod.data.test <- modelobj.save[[5]]
    mod.name <- modelobj.save[[6]]
    mod.type <- modelobj.save[[7]]
    mod.method <- modelobj.save[[8]]
    
    a.mod <- mod.cv$finalModel
    a.data <- mod.data.tval
    
    #predict
    p.test <- predict(a.mod, type='response')
    
    #ROC curve
    p.test <- predict(a.mod, type='prob')
    pred <- prediction(p.test[,'high'], a.data$cyt.bool)
    perf <- performance(pred, "tpr", "fpr")
    perf.auc <- performance(pred, "auc")
    auc.val <- round(perf.auc@y.values[[1]],2)
    if(n==1)
      plot(perf, col=n, lwd=2)#, main=paste0(mod.name,'; AUC: ', auc.val, '; ', a.name))
    else
      plot(perf, col=n, lwd=2, add=TRUE)
    
    #update legend with the AUC values
    mod.legnames[n] <- paste0(mod.legnames[n], '; AUC=', auc.val)
    
    n = n+1
  }
  abline(0,1,col='gray')
  legend('bottomright', legend=mod.legnames, col=seq(1, n-1), lty=1, bty='n')
  title(paste0(cancertype, ' - ', a.name))
  dev.off()
}

#------- aggregate model summaries
res <- summarize.combineResults(cancertypes)
res.model <- res[[2]]
all.res <- lapply(cancertypes, analyze.res.model, res.model) #loop through each cancer and analyze data

#------- aggregate important variables
df.pval <- NULL
df.val <- NULL
for(cancertype in cancertypes){
  dirout <- file.path(dir.out, paste0('_',cancertype))
  
  df.tmp <- read.csv(file.path(dirout, 'bin_rf_mut+NMD_fl_RFimp_scaled.csv'), stringsAsFactors = FALSE)
  df.tmp.pval <- df.tmp %>% select(c(X, MeanDecreaseAccuracy.pval))
  df.tmp.val <- df.tmp %>% select(c(X, MeanDecreaseAccuracy))
  colnames(df.tmp.pval) <- c('X', cancertype)
  colnames(df.tmp.val) <- c('X', cancertype)
  
  if(is.null(df.pval))
    df.pval <- df.tmp.pval
  else
    df.pval <- df.pval %>% inner_join(df.tmp.pval, by='X')
  
  if(is.null(df.val))
    df.val <- df.tmp.val
  else
    df.val <- df.val %>% inner_join(df.tmp.val, by='X')
}

df.pval <- df.pval[match(c("silent", "missense", "nonsense", "nonstop", "frameshift", "total",
                           "nmdfs_frac_decayed", "nmdfs_max", "nmdfs_mean",
                           "nmdfs_mean.wt", "nmdfs_med", "nmdfs_med.wt", "nmdfs_mut.bool", "nmdfs_n_decayed",     
                           "nmdns_frac_decayed", "nmdns_max", "nmdns_mean", "nmdns_mean.wt", "nmdns_med", 
                           "nmdns_med.wt", "nmdns_mut.bool", "nmdns_n_decayed", "nmdptc_frac_decayed", "nmdptc_max", 
                           "nmdptc_mean", "nmdptc_mean.wt", "nmdptc_med", "nmdptc_med.wt", "nmdptc_n_decayed"), df.pval$X),]

df.val <- df.val[match(c("silent", "missense", "nonsense", "nonstop", "frameshift", "total",
                           "nmdfs_frac_decayed", "nmdfs_max", "nmdfs_mean",
                           "nmdfs_mean.wt", "nmdfs_med", "nmdfs_med.wt", "nmdfs_mut.bool", "nmdfs_n_decayed",     
                           "nmdns_frac_decayed", "nmdns_max", "nmdns_mean", "nmdns_mean.wt", "nmdns_med", 
                           "nmdns_med.wt", "nmdns_mut.bool", "nmdns_n_decayed", "nmdptc_frac_decayed", "nmdptc_max", 
                           "nmdptc_mean", "nmdptc_mean.wt", "nmdptc_med", "nmdptc_med.wt", "nmdptc_n_decayed"), df.val$X),]

write.csv(df.pval, file=file.path(dir.out, 'summaries', 'mut+NMD_pval.csv'))
write.csv(df.val, file=file.path(dir.out, 'summaries', 'mut+NMD_val.csv'))

df2 <- df.val %>% select(-X)
df2.std <- as.data.frame(apply(df2, 2, scale))
rownames(df2.std) <- df.val$X
pdf(file.path(dir.out, 'summaries', paste0('heatmap_accdesc_std.pdf')))
heatmap.2(as.matrix(df2.std), 
          col=rev(colorRampPalette(brewer.pal(9, "RdBu"))(100)), scale="none", 
          key.xlab='mean desc in accuracy; std', key.title=F,
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5, cexCol=0.5,
          margins = c(8,8))
dev.off()

rownames(df2) <- df.val$X
pdf(file.path(dir.out, 'summaries', paste0('heatmap_accdesc.pdf')))
heatmap.2(as.matrix(df2), 
          col=rev(colorRampPalette(brewer.pal(9, "RdBu"))(100)), scale="none", 
          key.xlab='mean desc in accuracy', key.title=F,
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5, cexCol=0.5,
          margins = c(8,8))
  #Rowv=FALSE, dendrogram='column'
dev.off()

#pval with cut-off
pval_cutoff <- 0.1
df2 <- df.pval %>% select(-X)
df2 <- df2 < pval_cutoff
df2 <- as.data.frame(apply(df2, 2, as.numeric))
rownames(df2) <- df.pval$X
pdf(file.path(dir.out, 'summaries', paste0('heatmap_accdesc_pval', pval_cutoff, '.pdf')))
heatmap.2(as.matrix(df2), 
          col=c("#ffffff", "#000000"), scale="none", 
          key.xlab=paste0('Boolean of p-value < ', pval_cutoff), key.title=F,
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5, cexCol=0.5,
          margins = c(8,8))
dev.off()


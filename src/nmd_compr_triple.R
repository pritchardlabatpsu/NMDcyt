#compare three categories (cooccur, 1var, and ctrl), the sig test is based on trend (three cateogries)
rm(list=ls())
setwd("/home/boyangzhao/NMDcyt/")

#libraries
library(Hmisc)
library(gplots)
library(RColorBrewer)
library(ggpubr)
library(reshape2)
library(dplyr)
source('./src/utils.R')

library(PMCMR)

#load dataset manually
load("./outputs/19.0110_cnacutoff_1_NMDpos/results.rda")

#settings
dirout.p <- './outputs/19.0110_cnacutoff_1_NMDpos_CooccurVs1variantVsCtrl/' #the old and new analyses used dir.out and dirout slightly differently, this is added to avoid mixing up the usage
dirout <- dirout.p
if(!file.exists(dirout)) dir.create(dirout)
df.NMD.mutcna2 <- df.NMD.any[df.NMD.any$patientid %in% df.mutampN[df.mutampN$n>1,'patientid'],]
df.NMD.mutcna1 <- df.NMD.any[df.NMD.any$patientid %in% df.mutampN[df.mutampN$n<2,'patientid'],]
NMD.ref.name <- 'NMD_ctrl'
df.NMD.ref <- df.NMD.ctrl #reference dataset to be used for the comparison, e.g. background control
toAnalyze <- list(NMD.any.cooccur=df.NMD.mutcna2) #NMD with gain and/or loss of function, only co-occurred cna/mut
cancertypes <- c('BLCA','BRCA','CESC','COADREAD','GBMLGG','HNSC','KIPAN',
                 'LUAD','LUSC','OV','PRAD','SKCM','STAD','THCA','UCEC')

for(anlyzname in names(toAnalyze)){
  df.NMD <- toAnalyze[[anlyzname]]
  if(nrow(df.NMD) < 1) next
  
  #correlation between NMD and NMD ctrl metrics, per cancer type
  nmdstat <- data.frame()
  nmdstat.log <- data.frame()
  for(cancertype in cancertypes){
    for(nmdcolname in nmdcolnames){
      
      if(cancertype == 'pancancer'){
        nmdval <- df.NMD[nmdcolname]
        nmdval.1var <- df.NMD.mutcna1[nmdcolname]
        nmdval.ctrl <- df.NMD.ref[nmdcolname]
      } else {
        nmdval <- df.NMD[df.NMD$cancertype==cancertype,nmdcolname]
        nmdval.1var <- df.NMD.mutcna1[df.NMD.mutcna1$cancertype==cancertype,nmdcolname]
        nmdval.ctrl <- df.NMD.ref[df.NMD.ref$cancertype==cancertype,nmdcolname]
      }
      if(!all(is.na(nmdval)) && length(nmdval) > 0 && 
         !all(is.na(nmdval.1var)) && length(nmdval.1var) > 0 && 
         !all(is.na(nmdval.ctrl)) && length(nmdval.ctrl) > 0){
        df1 <- data.frame(val=nmdval, type=anlyzname, lab=2)
        colnames(df1) <- c('val','type','lab')
        df2 <- data.frame(val=nmdval.1var, type='NMD_1variant', lab=1)
        colnames(df2) <- c('val','type','lab')
        df3 <- data.frame(val=nmdval.ctrl, type='NMD_ctrl', lab=0)
        colnames(df3) <- c('val','type','lab')
        df <- rbind(df1, df2, df3)
        
        val <- wilcox.test(df1$val, df3$val)
        diff <- median(df1$val, na.rm=TRUE) - median(df3$val, na.rm=TRUE)
        jt.test <- jonckheere.test(df$val, df$lab, alternative="monotonic") #Jonckheere-Terpstrata Test
        tmp <- data.frame(cancertype=cancertype, nmdmetric=nmdcolname, 
                          diff=diff, diff.pval=val$p.value, jt.pval=jt.test$p.value)
        if(nrow(nmdstat) < 1) { nmdstat <- tmp } else { nmdstat <- rbind(nmdstat, tmp) }
        
        r <- tryCatch({
          val <- wilcox.test(log10(df1$val), log10(df3$val)) #Wilcox test of difference of NMD metric
          diff <- median(log10(df1$val), na.rm=TRUE) - median(log10(df3$val), na.rm=TRUE) #median difference of NMD metric
          df.calc <- df
          df.calc$val <- log10(df$val)
          df.calc <- df.calc[!is.nan(df.calc$val) & !is.na(df.calc$val) & df.calc$val != Inf & df.calc$val != -Inf,]
          jt.test <- jonckheere.test(df.calc$val, df.calc$lab, alternative="monotonic") #Jonckheere-Terpstrata Test
          tmp <- data.frame(cancertype=cancertype, nmdmetric=nmdcolname, 
                            diff=diff, diff.pval=val$p.value, jt.pval=jt.test$p.value)
          if(nrow(nmdstat.log) < 1) { nmdstat.log <- tmp } else { nmdstat.log <- rbind(nmdstat.log, tmp) }
        }, error=catch.errorfunc, finally={})
      }
      
    }
  }
  
  nmdstat$diff.qval <- p.adjust(nmdstat$diff.pval, method="BH")
  nmdstat$jt.qval <- p.adjust(nmdstat$jt.pval, method="BH")
  nmdstat.log$diff.qval <- p.adjust(nmdstat.log$diff.pval, method="BH")
  nmdstat.log$jt.qval <- p.adjust(nmdstat.log$jt.pval, method="BH")
  
  #Keep the statistically significant ones and plots those
  nmdstat.log.sig <- nmdstat.log[nmdstat.log$jt.qval < qval_cutoff,]
  
  plotBarsNMD <- function(x){
    dirout <- paste0(dirout.p, '/sig', qval_cutoff,'_', anlyzname, '/')
    if(!file.exists(dirout)) dir.create(dirout)
    cancertype <- x['cancertype']
    nmdcolname <- x['nmdmetric']
    
    if(cancertype == 'pancancer'){
      nmdval <- df.NMD[nmdcolname]
      nmdval.1var <- df.NMD.mutcna1[nmdcolname]
      nmdval.ctrl <- df.NMD.ref[nmdcolname]
    } else {
      nmdval <- df.NMD[df.NMD$cancertype==cancertype,nmdcolname]
      nmdval.1var <- df.NMD.mutcna1[df.NMD.mutcna1$cancertype==cancertype,nmdcolname]
      nmdval.ctrl <- df.NMD.ref[df.NMD.ref$cancertype==cancertype,nmdcolname]
    }
    if(length(nmdval) == 1 && is.na(nmdval)) return()
    if(length(nmdval.ctrl) == 1 && is.na(nmdval.ctrl)) return()
    if(length(nmdval) > 0 && length(nmdval.ctrl) > 0){
      df1 <- data.frame(val=nmdval, type=anlyzname)
      colnames(df1) <- c('val','type')
      df2 <- data.frame(val=nmdval.1var, type='NMD_1variant')
      colnames(df2) <- c('val','type')
      df3 <- data.frame(val=nmdval.ctrl, type='NMD_ctrl')
      colnames(df3) <- c('val','type')
      df <- rbind(df1, df2, df3)
      
      r <- tryCatch({
        ggplot(df, aes(as.factor(df$type), val)) + geom_boxplot() + geom_jitter(width=0.2) +
          xlab(cancertype) + ylab(nmdcolname) + scale_y_log10() +
          stat_compare_means(method='wilcox.test',comparisons=list(c(anlyzname,'NMD_ctrl'), c(anlyzname,'NMD_1variant')))
        ggsave(file.path(dirout, paste0('nmdlog10_', cancertype,'_',nmdcolname,'.pdf')))
        #note since using scale transform and not coord transform, the statistics for stat_compare is test on the
        #transformed data and not the raw data values
      }, error=catch.errorfunc, finally={})
    }
  }
  
  if(nrow(nmdstat.log.sig)>0) apply(nmdstat.log.sig, 1, plotBarsNMD)
  write.csv(nmdstat.log, file=paste0(dirout, '/nmdstatlog_', anlyzname,'.csv'))
}


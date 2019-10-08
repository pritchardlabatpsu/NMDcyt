# Analyze null model results
#  - generate averaged ROC curves and bar plots
# Model results are saved in as .rda - to be used as input here
rm(list=ls())
setwd("/home/boyangzhao/NMDcyt/")

#libraries
if(!require("pacman")) install.packages("pacman")
pacman::p_load(caret, car, MASS, ROCR, randomForest, rfPermute)
pacman::p_load(gplots, ggplot2, ggrepel, GGally, RColorBrewer, ggpubr, factoextra)
pacman::p_load(Hmisc)
pacman::p_load(plyr, dplyr)

#------- parameters --------
dir.out.base <- './out/19.0817_nullmodel' #the base output directory name from analyze_nullmodel.R
dir.db <- '../datasets/' #location of the datasets
pval.ptc.excludeNoMut <- TRUE #exclude non mut patients, for plotting/mann-whitney test on NMD against discretized CYT
param.parallel <- TRUE #use futures, for parallel computing

cancertypes <- c('BLCA','BRCA','CESC','COADREAD','GBMLGG','HNSC','KIPAN',
                 'LUAD','LUSC','OV','PRAD','SKCM','STAD','THCA','UCEC')

if(!file.exists(dir.out.base)) dir.create(dir.out.base)

#------- aggregate model summaries (bar plots of AUCs and OOB errors)
seedsN <- 10 #number of randomization datasets used
seedsSeq <- seq(1,seedsN)
for(cancertype in cancertypes){
  res.model <- data.frame() #multivariate model results
  for(n in seedsSeq){
    #aggregate all the randomized results per given cancer type
    dir.out <- paste0(dir.out.base, n)
    f <- file.path(dir.out, paste0('_',cancertype), 'summary_multivariate.csv')
    if(file.exists(f)){
      #append res.model
      res.model.s <- read.csv(f)
      res.model.s <- res.model.s[,1:3]
      colnames(res.model.s)[2:3] <- paste0(colnames(res.model.s)[2:3],'_',n)
      if(nrow(res.model) < 1) res.model <- res.model.s
      else res.model <- cbind(res.model, res.model.s[3])
    }
  }
  res.model$model_name <- sub('.*_','',res.model$model_name)
  res.model$model_name <- factor(res.model$model_name, res.model$model_name)
  colnames(res.model)[2] <- 'model'
  
  getAUCCV <- function(x){
    y <-  sub(';auc\\.cv\\.sd:.*','',x)
    y <- sub('.*auc\\.cv:', '', y)
    y <- as.numeric(y)
    return(y)
  }
  getOOB <- function(x){
    y <- sub('.*oob:','',x)
    y <- as.numeric(y)
    return(y)
  }
  end_idx <- 3+length(seedsSeq)-1
  vals <- res.model[,3:end_idx]
  vals.auc <- sapply(vals,getAUCCV)
  colnames(vals.auc) <- paste0('AUC_',seedsSeq)
  vals.oob <- sapply(vals,getOOB)
  colnames(vals.oob) <- paste0('OOB_',seedsSeq)
  
  auc.mean <- apply(vals.auc,1,function(x) mean(x, na.rm = T))
  auc.sd <- apply(vals.auc,1,function(x) sd(x, na.rm = T))
  auc.se <- auc.sd/sqrt(length(seedsSeq))
  auc.stats <- data.frame(auc.mean=auc.mean, auc.sd=auc.sd, auc.se=auc.se)
  
  oob.mean <- apply(vals.oob,1,function(x) mean(x, na.rm = T))
  oob.sd <- apply(vals.oob,1,function(x) sd(x, na.rm = T))
  oob.se <- oob.sd/sqrt(length(seedsSeq))
  oob.stats <- data.frame(oob.mean=oob.mean, oob.sd=oob.sd, oob.se=oob.se)
  
  df <- cbind(res.model[,1:2], vals.auc, auc.stats, vals.oob, oob.stats)
  
  #deal with MSI, if MSI is empty, don't need to include models that are also MSI, redundant
  if(!any(grepl('msi',df$model))){
    df <- df[!grepl("+MSI$",df$model_name),]
  }
  
  #AUROC CV +/- SE(AUROC CV) of the n number of randomizations
  ggplot(data=df, aes(x=model_name, y=auc.mean)) + 
    geom_bar(stat='identity') + geom_errorbar(aes(ymin=auc.mean-auc.se, ymax=auc.mean+auc.se), width=.1) +
    ylab('AUROC (null model)') + xlab('Models') + 
    coord_cartesian(ylim=c(0.4,0.9))
  ggsave(file.path(dir.out.base, paste0('auroc_bar_', cancertype, '.pdf')))
  
  #OOB+/- SE(OOB) of the n number of randomizations
  ggplot(data=df, aes(x=model_name, y=oob.mean)) + 
    geom_bar(stat='identity') + geom_errorbar(aes(ymin=oob.mean-oob.se, ymax=oob.mean+oob.se), width=.1) +
    ylab('OOB (null model)') + xlab('Models') + 
    coord_cartesian(ylim=c(0.0,0.65))
  ggsave(file.path(dir.out.base, paste0('oob_bar_', cancertype, '.pdf')))
}

#------- combined ROC curve (averaged across randomized datasets)
x.fpr <- seq(0,1,0.0001) #fixed false positive rate (x values for ROC curve)

getROC <- function(f, x.fpr=seq(0,1,0.0001)){
  #f: filename
  load(f)
  a.mod <- modelobj.save[[2]]$finalModel #cv best model
  a.data <- modelobj.save[[5]] #test data, for the analysis done, this is the whole dataset
  
  p.test <- predict(a.mod, type='prob')
  pROC_obj <- roc(a.data$cyt.bool, p.test[,'high'])
  x.sp <- 1 - x.fpr
  y.tpr <- coords(pROC_obj, x = x.sp, input = "specificity", ret = "se")[1,]
  df <- data.frame(fpr=x.fpr, tpr=y.tpr)
  return(df)
}

for(cancertype in cancertypes){
  df.mut <- data.frame(fpr=x.fpr,tpr=rep(0,length(x.fpr)), mod='Mut')
  df.NMD <- data.frame(fpr=x.fpr,tpr=rep(0,length(x.fpr)), mod='NMD')
  df.mutNMD <- data.frame(fpr=x.fpr,tpr=rep(0,length(x.fpr)), mod='Mut+NMD')
  for(n in seedsSeq){
    #aggregate all the randomized results per given cancer type
    dir.out <- paste0(dir.out.base, n)
    df.mut.s <- getROC(file.path(dir.out, paste0('_',cancertype), 'bin_rf_mut_fl.rda'), x.fpr)
    df.NMD.s <- getROC(file.path(dir.out, paste0('_',cancertype), 'bin_rf_NMD_fl.rda'), x.fpr)
    df.mutNMD.s <- getROC(file.path(dir.out, paste0('_',cancertype), 'bin_rf_mut+NMD_fl.rda'), x.fpr)
    
    df.mut$tpr = df.mut$tpr+df.mut.s$tpr
    df.NMD$tpr = df.NMD$tpr+df.NMD.s$tpr
    df.mutNMD$tpr = df.mutNMD$tpr+df.mutNMD.s$tpr
  }
  df.mut$tpr = df.mut$tp/length(seedsSeq)
  df.NMD$tpr = df.NMD$tpr/length(seedsSeq)
  df.mutNMD$tpr = df.mutNMD$tpr/length(seedsSeq)
  
  df <- rbind(df.mut, df.NMD, df.mutNMD)
  
  ggplot(df, aes(x=fpr, y=tpr, group=mod, color=mod)) + geom_line(size=1) + 
    xlab('False positive rate') + ylab('True positive rate') + theme_bw() +
    geom_abline(intercept = 0, color='gray') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_color_manual(values=c('black', 'red', 'green'))
  ggsave(file.path(dir.out.base, paste0('ROC_', cancertype, '.pdf')))
}


# Aggregate survival analysis results
rm(list=ls())
setwd("/home/boyangzhao/NMDcyt/")

#libraries
if(!require("pacman")) install.packages("pacman")
pacman::p_load(caret, car, MASS, ROCR, randomForest, rfPermute)
pacman::p_load(gplots, ggplot2, ggrepel, GGally, RColorBrewer, ggpubr, factoextra)
pacman::p_load(plyr, Hmisc)
pacman::p_load(dplyr)
pacman::p_load(future)
source('src/utils.R')
source('src/utils_surv.R')

#------- parameters --------
dir.out <- './outputs/18.1212_survival/' #output directory
dir.db <- './datasets/' #location of the datasets
pval.ptc.excludeNoMut <- TRUE #exclude non mut patients, for plotting/mann-whitney test on NMD against discretized CYT
param.parallel <- TRUE #use futures, for parallel computing

#------- retrieve/parse datasets --------
cancertypes <- c('BLCA','BRCA','CESC','COADREAD','GBMLGG','HNSC','KIPAN',
                 'LUAD','LUSC','OV','PRAD','SKCM','STAD','THCA','UCEC')

#pan-cancer table
createPanCancer(cancertypes, reanalyze=TRUE, appendIndication=TRUE, fname='pancancer_m')
load(file.path(dir.out, 'pancancer_m.rda'))
xvars <- c('silent','missense','nonstop','nonsense','frameshift','total', colnames(df.merged)[grep('^nmd' ,colnames(df.merged))], 'cyt')

#pan-cancer statistics
df.stats <- df.merged %>%
  select(c(xvars, 'cancertype')) %>%
  group_by(cancertype) %>% summarize_all(median) %>% collect

#get median survival per cancer type
clin.stats <- data.frame()
for(cancertype in cancertypes){
  sf <- survfit(Surv(duration, event)~1, data=df.merged.clin[df.merged.clin$cancertype==cancertype,])
  sfit.tbl <- summary(sf)$table
  if('strata' %in% names(sf)){
    res <- data.frame(strata=rownames(sfit.tbl),
                      n=sfit.tbl[,'n.start'], n_event=sfit.tbl[,'events'],
                      medsurv=sfit.tbl[,'median'],
                      ci95L=sfit.tbl[,'0.95LCL'], ci95H=sfit.tbl[,'0.95UCL'])
  } else {
    res <- data.frame(strata='-',
                      n=sfit.tbl['n.start'], n_event=sfit.tbl['events'],
                      medsurv=sfit.tbl['median'],
                      ci95L=sfit.tbl['0.95LCL'], ci95H=sfit.tbl['0.95UCL'])
  }
  
  df <- data.frame(cancertype=cancertype, medsurv=res$medsurv, n=res$n, stringsAsFactors=FALSE)
  clin.stats <- rbind(clin.stats, df)
}
clin.stats$cancertype <- as.character(clin.stats$cancertype)

#merge
df.stats <- df.stats %>% left_join(clin.stats, by='cancertype') %>% collect
df.stats <- as.data.frame(df.stats)

#------- analyses --------
dirout <- file.path(dir.out, '/summaries_aggr')
if(!file.exists(dirout)) dir.create(dirout)
for(xvar in xvars){
  if(pval.ptc.excludeNoMut && grepl("^nmd.*", xvar) && !grepl(".*_mut\\.bool$", xvar)){
    #if the xvar is a NMD metric, and pval.ptc.excludeNoMut is TRUE, then will exclude patients with no corresponding
    #mutation in the plotting and calculation of p-val (based on Mann-Whitney test)
    xvar.nmd <- sub("_.*$",'',xvar) #get which PTC to look at (fs, ns, or ptc)
    df.tmp <- df.merged[df.merged[paste0(xvar.nmd, '_mut.bool')]==1,] #df with just patients with mutations for the PTC type
  } else {
    df.tmp <- df.merged
  }
  
  #box-whisker plot
  ggplot(df.tmp, aes(x=reorder(cancertype, eval(parse(text=xvar)), function(x)1-median(log10(x[x>0]))), y=eval(parse(text=xvar)) )) + 
    geom_boxplot() + scale_y_log10() + theme(axis.title.x=element_blank()) + ylab(xvar)
  ggsave(file.path(dirout, paste0('bw_', xvar, '.pdf')))
  
  #scatter OS vs parameter
  # ggplot(df.stats, aes(eval(parse(text=xvar)), medsurv, label=cancertype)) + geom_point(aes(size=n)) + 
  #   geom_text_repel() + stat_smooth(method=lm, aes(weight=n)) +
  #   ylab('Median survival (years)') + xlab(paste0('Median ', xvar)) + scale_size_continuous("No. of patients") + theme_classic(base_size = 15)
  # ggsave(file.path(dirout, paste0('corr_OS_', xvar, '.pdf')))
  
  #scatter cyt vs parameter
  if(xvar=='cyt') next
  ggplot(df.stats, aes(eval(parse(text=xvar)), cyt, label=cancertype)) + geom_point(aes(size=n)) + 
    geom_text_repel() + stat_smooth(method=lm, aes(weight=n)) +
    ylab('Cytolytic activity (median)') + xlab(paste0(xvar, ' (median)')) + scale_size_continuous("No. of patients") + theme_classic(base_size = 15)
  ggsave(file.path(dirout, paste0('corr_cyt_', xvar, '.pdf')))
}

#correlation matrix - in aggregate by cancer type
df <- df.stats %>% select(-medsurv,-n, -nmdfs_mut.bool, -nmdns_mut.bool, -nmdptc_mut.bool, -cancertype, -nonstop)
c <- rcorr(sapply(df, as.numeric), type="spearman")
r <- tryCatch({
  pdf(file.path(dirout, paste0('corr_matrix_aggcancer.pdf')))
  heatmap.2(as.matrix(c$r), 
            col=colorRampPalette(brewer.pal(9, "RdBu"))(100), scale="none", 
            key.xlab='Spearman corr', key.title=F,
            key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.8, cexCol=0.8,
            margins = c(8,8))
  dev.off()
}, error=catch.errorfunc, finally={})

#------- combine the univariate Cox regression results (table)
library(rJava)
library(xlsx)

wb1 <- createWorkbook()
for(cancertype in cancertypes){
  dirout <- file.path(dir.out, paste0('_',cancertype))
  df.tmp <- read.csv(file.path(dirout, 'uni.cox2', 'univariate_summary_cox.csv'), stringsAsFactors = FALSE)
  
  #include only the relevant rows
  selbool <- !(grepl('_mut$',df.tmp$var) | grepl('_cna$',df.tmp$var))
  df.tmp <- df.tmp[selbool,]
  colnames(df.tmp)[1] <- 'var'

  df.tmp$pval<-as.numeric(format(df.tmp$pval, scientific=TRUE)) #convert
  df.tmp$qval <- p.adjust(df.tmp$pval, 'fdr')
  
  sheet1 <- createSheet(wb1, sheetName=cancertype)
  addDataFrame(df.tmp, sheet1, row.names=FALSE)
  
}

saveWorkbook(wb1, file.path(dir.out, 'summaries_aggr', "surv_uni_cox_combined.xlsx"))

#------- combine multivariate Cox regression results into Excel
library(rJava)
library(xlsx)

wb1 <- createWorkbook()
for(cancertype in cancertypes){
  dirout <- file.path(dir.out, paste0('_',cancertype))
  sheet1 <- createSheet(wb1, sheetName=cancertype)
  addDataFrame(read.csv(file.path(dirout, 'uni.cox.withctrl2', 'univariate_cox_withctrl_TMBPDL1cat.csv')), sheet1, row.names=FALSE)
}

#saving the workbook
saveWorkbook(wb1, file.path(dir.out, 'summaries_aggr', "surv_cox_withctrl_TMBPDL1cat.xlsx"))


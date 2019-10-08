# Additional individual analyses
# for Fig 2C and Supp Fig S6
# redid to simplify the graphics, to only group by cancer type
rm(list=ls())
setwd("/home/boyangzhao/NMDcyt/")

if(!require("pacman")) install.packages("pacman")
pacman::p_load(caret, car, MASS, ROCR, randomForest, rfPermute)
pacman::p_load(gplots, ggplot2, GGally, RColorBrewer, ggpubr, factoextra)
pacman::p_load(plyr, Hmisc)
pacman::p_load(dplyr)

dir.in <- './out/19.0110_cnacutoff_1_NMDpos_CooccurVsCtrl_copy/'
dir.out <- './out/19.0820_NMDglobal/'
qval_cutoff <- 0.05
if(!file.exists(dir.out)) dir.create(dir.out)

genGlobalNMDplot <- function(df, filename){
  c1 <- df[df$qval < qval_cutoff,'cancertype'] %>% unique()
  ckSame <- function(x){all(c(1,3)> 0) | all(c(1,3)<0)} #check for the same directionality
  c2 <- df[df$qval < qval_cutoff,] %>% group_by(cancertype) %>% select(cancertype,diff) %>% summarize(consistent=ckSame(diff)) %>% as.data.frame()
  c2 <- c2[c2$consistent,'cancertype']
  cancertype.sel <- union(c1,c2)
  df.sel <- df[df$cancertype %in% cancertype.sel,]
  df.sel$sig <- df.sel$qval<qval_cutoff
  ggplot(df.sel,aes(x=cancertype,diff, color=sig)) + geom_point(size=3, position=position_jitter(width=0.2, height=0.2)) + 
    geom_hline(yintercept = 0) +
    labs(x="Indication_NMDmetric", y="delta median (NMD - NMD ctrl)", color=paste0("qval < ",qval_cutoff)) +
    scale_color_manual(values=c("gray","red"))
  ggsave(file.path(dir.out, filename), width=12, height=6)
  
}

df <- read.csv(file.path(dir.in, 'nmdstatlog_NMD.any.cooccur.csv'))
genGlobalNMDplot(df, 'plot_nmdstatlog_NMD.any.cooccur_sigonly.pdf')

df <- read.csv(file.path(dir.in, 'nmdstatlog_NMD.amp.csv'))
genGlobalNMDplot(df, 'plot_nmdstatlog_NMD.amp_sigonly.pdf')

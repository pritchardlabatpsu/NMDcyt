#Cox regression (SKCM nmdptc_med)
#plot Schoenfeld residuals

rm(list=ls())
setwd("/home/boyangzhao/NMDcyt/")

if(!require("pacman")) install.packages("pacman")
pacman::p_load(caret, car, MASS, ROCR, randomForest, rfPermute)
pacman::p_load(gplots, ggplot2, GGally, RColorBrewer, ggpubr, factoextra)
pacman::p_load(plyr, Hmisc)
pacman::p_load(dplyr)
source('src/utils.R')
source('src/utils_surv.R')
library(gridExtra)

dir.in <- './out/18.1212_survival_copy/' #input directory of previously analyzed survival results, *_prep_base.rda
dir.output <- './out/19.0819_cox/' #output directory
cancertype <- 'SKCM'

#------- portions code reused from utils and utils_surv, for ad-hoc analyses here --------
load(file.path(dir.in, paste0(cancertype, '_prep_base.rda')))

dir.out <- dir.in #a hack, for retrieveRNASeqGene function below, which uses dir.out
df.rnaseq.PDL1 <- retrieveRNASeqGene(cancertype,'CD274')
df.merged.clin <- append.clinVars(df.merged.clin, df.rnaseq.PDL1)

#corrections - age
if('age' %in% colnames(df.merged.clin)){
  df.merged.clin$age <- df.merged.clin$age*365.2 #the age was divided by 365.2, correct it here
  
  #categorize age
  df.merged.clin$age.cat <- '<65'
  df.merged.clin$age.cat[df.merged.clin$age >= 65] <- '>=65'
  df.merged.clin$age.cat <- factor(df.merged.clin$age.cat, levels=c('<65','>=65'))
}

#define variables to look at
varnames <- c(colnames(df.merged), 'tnm_stage', 'gender','race','age.cat', 'TMB', 'TMB.cat', 'PDL1', 'PDL1.cat', 'PDL1.cat0')

#determine which variables are factor/numeric
var.cat <- retrieveVarTypeNames(df.merged.clin[, colnames(df.merged.clin) %in% varnames], is.factor)
var.cont <- retrieveVarTypeNames(df.merged.clin[, colnames(df.merged.clin) %in% varnames], is.numeric)

var.cat <- c(var.cat, var.cont[grepl('\\.bool$',var.cont)]) #take care of .bool labels, these are categorical, not continuous
var.cont <- var.cont[!grepl('\\.bool$',var.cont)] #and remove them from the continuous list here

if('nmdns_mut.bool' %in% colnames(df.merged.clin)) df.merged.clin$nmdns_mut.bool <- factor(df.merged.clin$nmdns_mut.bool)
if('nmdfs_mut.bool' %in% colnames(df.merged.clin)) df.merged.clin$nmdfs_mut.bool <- factor(df.merged.clin$nmdfs_mut.bool)
if('nmdptc_mut.bool' %in% colnames(df.merged.clin)) df.merged.clin$nmdptc_mut.bool <- factor(df.merged.clin$nmdptc_mut.bool)

var.cont.km <- c('silent','missense','nonstop','nonsense','frameshift','total', 'GZMA', 'PRF1', 'cyt',
                 var.cont[grep('^nmd' , var.cont)])
for(varname in var.cont.km){
  if(!varname %in% colnames(df.merged.clin)) next
  varname.cat <- paste0(varname, ".cat")
  q1 <- quantile(df.merged.clin[[varname]])['25%']
  q3 <- quantile(df.merged.clin[[varname]])['75%']
  df.merged.clin[varname.cat] <- 'med'
  df.merged.clin[df.merged.clin[varname] <= q1, varname.cat] <- 'low'
  df.merged.clin[df.merged.clin[varname] >= q3, varname.cat] <- 'high'
  df.merged.clin[varname.cat] <- factor(df.merged.clin[[varname.cat]], levels=c('low','med','high'))
}
var.cont.km <- paste0(var.cont.km, '.cat')

#multivariate Cox regression - set up
outpref.plot <<- FALSE
varnames <- union(union(var.cat, var.cont),var.cont.km)
varnames <- varnames[!varnames %in% c('age.cat','gender','tnm_stage','race')] #remove the variables to control

#------- Cox regression --------
#multivariate Cox regression - controlling for baselines: age, gender, and tnm_stage, TMB, PDL1
#plot the Schoenfeld resdiuals
dir.out <- dir.output
res.cox <- coxph(Surv(duration, event) ~ nmdptc_med+age.cat+gender+tnm_stage+TMB.cat+PDL1.cat, data = df.merged.clin)
summary(res.cox)
test.ph <- cox.zph(res.cox)
test.ph
plot_zph = ggcoxzph(test.ph)
ggsave(file.path(dir.out, paste0('residual_Schoenfeld_SKCM_nmdptc_med_TMBPBL1cat.pdf')), arrangeGrob(grobs = plot_zph))


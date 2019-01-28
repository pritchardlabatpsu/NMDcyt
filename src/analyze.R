# Mutations to Tumor Infiltration analysis
# Cytolytic activity assessed by GZMA and PRF1
rm(list=ls())
setwd("/home/boyangzhao/NMDcyt/")

#libraries
if(!require("pacman")) install.packages("pacman")
pacman::p_load(caret, car, MASS, ROCR, randomForest, rfPermute)
pacman::p_load(gplots, ggplot2, GGally, RColorBrewer, ggpubr, factoextra)
pacman::p_load(plyr, Hmisc)
pacman::p_load(dplyr)
pacman::p_load(future)
source('src/utils.R')
source('src/utils_surv.R')

#------- parameters --------
dir.out <- './outputs/18.0604/' #output directory
dir.db <- './datasets/' #location of the datasets
dbset.trainval <- 'all' #train/CV using a train.val dataset (values: train.val, all)
dbset.test <- 'all' #test using a test dataset (values: test, all)
dbset.preprocess.std <- FALSE #standardize data during preprocessing
dbset.preprocess.pca <- FALSE #perform PCA on data during preprocessing, if TRUE, then by default, data is standardized before PCA
dbset.preprocess.MSIcomplete <- FALSE #if true, only include cases with MSI-H and MSS, exclude any samples with missing info
pval.ptc.excludeNoMut <- TRUE #exclude non mut patients, for plotting/mann-whitney test on NMD against discretized CYT
anlyz.pancancer <- TRUE #generate and include pan-cancer in analysis
param.parallel <- FALSE #use futures, for parallel computing
genesList <- c(#list
  #'CD8A', 'CD4', 'ZBTB7B',  'TNFRSF18', 
  #tumor cell intrinsic: T cell exclusion
  'CD274', #PD-L1
  #tumor cell intrinsic: absense of antigen presentation
  'TAP1', 'B2M', 'HLA-A', 'HLA-B', 'HLA-C',
  #tumor cell intrinsic: insensibility to T cells
  'IFNG', 
  #tumor cell intrinsic: INF gamma signaling
  'IFNGR1', 'IFNGR2', 'JAK2', 'IRF1',
  'IDO1', 'CXCL9', #IDO1/CXCL9 based on Herbst et al, Nature
  'CXCL10', 'CXCL11',
  #tumor cell extrinsic: cytokines
  'IL10', 'TGFB1', 'IL2', 
  #tumor cell extrinsic: immunosuppressive cells
  'FOXP3',
  #tumor cell extrinsic: inhibitory immune checkpoints
  'HAVCR2','VSIR', 'LAG3',
  #tumor cell extrinsic:
  'PDCD1', 'CTLA4', 'TNFRSF4',
  #drivers MutSigCV, Rooney et al
  'CASP8', 'TP53', 'PIK3CA', 'NF1', 'MET', 'SOS1', 'SMC1A', 'TET2', 'ARID5B', 'ALPK2',
  'LPAR2', 'COL5A1', 'NCOR1', 'SSX5', 'DNER', 'MORC4', 'IRF6', 'MYOCD', 'CIC', 'SLC22A14',
  'CNKSR1', 'CUL4B', 'DDX3X', 'FUBP1', 'ARID2', 'TCP11L2', 'CSNK2A1', 'ASXL1', 'TMEM88',
  'DNMT3A', 'EP300', 'MUC17', 'OVOL1'
)

cancertypes <- c('BLCA','BRCA','CESC','COADREAD','GBMLGG','HNSC','KIPAN',
                 'LUAD','LUSC','OV','PRAD','SKCM','STAD','THCA','UCEC')

#------- main --------
if(!file.exists(dir.out)) dir.create(dir.out)
if(!file.exists(dir.db)) stop('Datasets directory not found...')

set.seed(1833)

#setting up parallel computing
if(param.parallel){
  plan(multiprocess) #multi-core
} else {
  plan(sequential) #single-core
}

#------- retrieve and parse datasets
cluster.retrievefunc <- function(cancertype){
  if(length(cancertype) > 1){
    #if input cancertype is a vector of multiple cancer types, then create a pancancer dataset based on this
    createPanCancer(cancertype, dirout=file.path(dir.out, '_pancancer')) #create pan-cancer
  } else {
    dirout <- file.path(dir.out, paste0('_',cancertype))
    if(!file.exists(dirout)) dir.create(dirout)
    print(sprintf("Retrieving datasets for %s...", cancertype))
    sink(file.path(dirout, paste0('out_get_', cancertype, '.txt')), append=FALSE)
    retrieve.dsets(cancertype, dirout) #retrieve data, results saved as rda file
    sink()
  }
}

all.res <- future_lapply(cancertypes, cluster.retrievefunc) #loop through each cancer and analyze data
if(anlyz.pancancer){
  cluster.retrievefunc(cancertypes) #create and analyze pancancer
  cancertypes <- c(cancertypes, 'pancancer') #append pancancer to the list of 'cancer types', for analysis etc below
}
if(!all(resolved(all.res))) stop('Parallel compute somehow not finished... terminating...')

#------- analyze each cancer type
cluster.analyzefunc <- function(cancertype){
  dirout <- file.path(dir.out, paste0('_',cancertype))
  if(!file.exists(dirout)) dir.create(dirout)
  
  print(sprintf("Analyzing %s...", cancertype))
  sink(file.path(dirout, paste0('out_anlyz', cancertype, '.txt')), append=FALSE)
  
  #preprocess datasets, results saved as rda file
  preprocess.dsets(cancertype, dirout)
  
  #analyze univariate/multivariate models
  analyze.dsets.uni(cancertype, dirout=dirout) #univariate
  analyze.dsets.multi(cancertype, dirout=dirout) #multivariate, model building
  
  #analyze clinical
# analyze.dsets.clin(cancertype, dirout=dirout) #clinical - kaplan-meier, cox regression, etc
  analyze.dsets.clin2(cancertype, dirout=dirout) #clinical - kaplan-meier, cox regression, etc with PDL1 and TMB
  
  sink()
  graphics.off()
}

all.res <- future_lapply(cancertypes, cluster.analyzefunc) #loop through each cancer and analyze data
if(!all(resolved(all.res))) stop('Parallel compute somehow not finished... terminating...')

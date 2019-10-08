# Mutations to Tumor Infiltration analysis
# Cytolytic activity assessed by GZMA and PRF1
# Null model based on randomized data
rm(list=ls())
setwd("/home/boyangzhao/NMDcyt/")

#libraries
if(!require("pacman")) install.packages("pacman")
pacman::p_load(caret, car, MASS, ROCR, randomForest, rfPermute)
pacman::p_load(gplots, ggplot2, GGally, RColorBrewer, ggpubr, factoextra)
pacman::p_load(plyr, Hmisc)
pacman::p_load(dplyr)
source('src/utils.R')

#------- parameters --------
dir.out <- './out/19.0817_nullmodel' #output directory
dir.in <- './out/18.0604_copy/' #input directory of previously analyzed _prep_base.rda
dir.db <- '../datasets/' #location of the datasets
dbset.trainval <- 'all' #train/CV using a train.val dataset (values: train.val, all)
dbset.test <- 'all' #test using a test dataset (values: test, all)
dbset.preprocess.std <- FALSE #standardize data during preprocessing
dbset.preprocess.pca <- FALSE #perform PCA on data during preprocessing, if TRUE, then by default, data is standardized before PCA
dbset.preprocess.MSIcomplete <- FALSE #if true, only include cases with MSI-H and MSS, exclude any samples with missing info
pval.ptc.excludeNoMut <- TRUE #exclude non mut patients, for plotting/mann-whitney test on NMD against discretized CYT
anlyz.pancancer <- TRUE #generate and include pan-cancer in analysis
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
if(!file.exists(dir.db)) stop('Datasets directory not found...')
dir.out.base <- dir.out
seeds <- c(12,593,3059,23,85, 212,2593,23059,223,285)
for(cancertype in cancertypes){
  for(n in seq(5,length(seeds))){
    set.seed(seeds[n])
    dir.out <- paste0(dir.out.base, n)
    if(!file.exists(dir.out)) dir.create(dir.out)
    
    #------- randomize dependent variables in data
    #get prep_base from dir.in and save a new prep_base in dir.out
    df_randomize <- function(df){
      varsY <- c('cyt','GZMA','PRF1','cyt.bool','cyt.cat')
      varsX <- colnames(df)
      varsX <- varsX[!(varsX %in% varsY)] #exclude the Y vars
      
      dfX <- df[varsX]
      dfY <- df[sample(nrow(df)),varsY]
      row.names(dfY) <- row.names(dfX)
      df.null <- cbind(dfX,dfY)
      return(df.null)
    }
    
    #load the original unrandomized data
    load(file.path(dir.in, paste0(cancertype, '_prep_base.rda')))
    
    df.merged <- df_randomize(df.merged)
    df.merged.bool <- df_randomize(df.merged.bool)
    df.merged.clin <- df_randomize(df.merged.clin)
    save(df.merged, df.merged.bool, df.merged.clin, file=file.path(dir.out, paste0(cancertype, '_prep_base.rda')))
    
    #------- analyze each cancer type
    dirout <- file.path(dir.out, paste0('_',cancertype))
    if(!file.exists(dirout)) dir.create(dirout)
    
    print(sprintf("Analyzing %s...", cancertype))
    sink(file.path(dirout, paste0('out_anlyz', cancertype, '.txt')), append=FALSE)
    
    #analyze univariate/multivariate models
    # analyze.dsets.uni(cancertype, dirout=dirout) #univariate
    analyze.dsets.multi(cancertype, dirout=dirout) #multivariate, model building
    
    sink()
    graphics.off()
  }
}

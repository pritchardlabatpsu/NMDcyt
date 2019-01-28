#compare co-occur with cytolytic activity
setwd("/home/boyangzhao/NMDcyt/")

#load dataset manually
load('./outputs/18.0604/pancancer.rda')
load("./outputs/19.0110_cnacutoff_1_NMDpos/results.rda")

nmdcolnames <- 'cyt'
df.cyt <- data.frame(patientid=row.names(df.merged), cyt=df.merged$cyt, stringsAsFactors=FALSE)
df.NMD.any <- df.NMD.any.all %>% select(cancertype, patientid) %>% unique() %>% left_join(df.cyt, by="patientid")

#get the control - no any NMD abberations (mutation nor CNA)
df.NMD.ctrl <- df.merged[!row.names(df.merged) %in% df.NMD.any$patientid,]
df.NMD.ctrl <- data.frame(patientid=row.names(df.NMD.ctrl),
                          cancertype=df.NMD.ctrl$cancertype,
                          df.NMD.ctrl[nmdcolnames], stringsAsFactors=FALSE)

#settings
df.NMD.mutcna1 <- df.NMD.any[df.NMD.any$patientid %in% df.mutampN[df.mutampN$n<2,'patientid'],]
df.NMD.mutcna2 <- df.NMD.any[df.NMD.any$patientid %in% df.mutampN[df.mutampN$n>1,'patientid'],]
toAnalyze <- list(NMD.any.cooccur=df.NMD.mutcna2) #NMD with gain and/or loss of function, only co-occurred cna/mut
NMD.ref.name <- 'NMD_ctrl'
df.NMD.ref <- df.NMD.ctrl #reference dataset to be used for the comparison, e.g. background control
cancertypes <- c('BLCA','BRCA','CESC','COADREAD','GBMLGG','HNSC','KIPAN',
                 'LUAD','LUSC','OV','PRAD','SKCM','STAD','THCA','UCEC')

dirout <- dirout.p <- './out/19.0110_cnacutoff_1_NMDpos_CooccurCytTrend/'
if(!file.exists(dirout)) dir.create(dirout)
#then run code for(anlyzname in names(toAnalyze)) in nmd_compr_triple.R

dirout <- dirout.p <- './out/19.0110_cnacutoff_1_NMDpos_CooccurCyt/'
if(!file.exists(dirout)) dir.create(dirout)
#then run code for(anlyzname in names(toAnalyze)) in nmd_compr.R
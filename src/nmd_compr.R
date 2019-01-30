#compare NMD metrics (of patients with NMD-gene alternations vs WT)
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

#------- parameters --------
dir.db <- '../datasets/' #location of the datasets
dir.out <- './outputs/'
dirout <- paste0(dir.out, '19.0110_cnacutoff_1_NMDpos/')
#NMD_genes <- c('SMG1', 'SMG5', 'SMG6', 'SMG7', 'SMG8', 'SMG9', 'UPF1', 'UPF2', 'UPF3A', 'UPF3B') #all
NMD_genes <- c('SMG1', 'SMG5', 'SMG6', 'SMG7', 'UPF1', 'UPF2', 'UPF3B') #positive regulators
cna_cutoff <- 1 #the CNA cut-off for analysis, values <(-1*val) or >val are considered CNA and used for analysis
cna_ctrl_cutoff <- 0 #the CNA cut-off with which to be used as no CNA, values not <(-1*val) or >val are considered as controls
qval_cutoff <- 0.05
rm.gainlossboth <- TRUE #remove patients with both gain and loss of function of the NMD genes
cancertypes <- c('BLCA','BRCA','CESC','COADREAD','GBMLGG','HNSC','KIPAN',
                 'LUAD','LUSC','OV','PRAD','SKCM','STAD','THCA','UCEC')

if(!file.exists(dirout)) dir.create(dirout)

#------- retrieve data --------
#Retrieve the mutations for all cancer types
df.NMD.lof <- data.frame() #for NMD genes loss of function (mutations or deletions)
df.NMD.amp <- data.frame() #for NMD genes amplifications (without overlapping mutations)
df.NMD.del <- data.frame() #for NMD genes deletions (without overlapping mutations)
df.NMD.any <- data.frame() #for NMD genes any alterations (gain or loss of function)
df.mutcna.n <- data.frame() #for all genes, with mut and/or cna
df.NMD.mut <- data.frame() #for NMD genes with mutations
df.NMD.cna <- data.frame() #for NMD genes with cna (!=0)
for(cancertype in cancertypes){
  print(paste0("Reading ", cancertype, '...'))
  
  #retrieve CNA
  #CNA values are reported as +/-2, +/-1, and 0
  fname <- paste0(cancertype,'_cna/')
  fname.txt <- 'all_thresholded.by_genes.txt'
  df.cna <- read.table(file.path(dir.db, fname, fname.txt),
                       header=T, stringsAsFactors=F, fill=T, sep='\t')
  pids <- sub('(\\.[[:alnum:]]*){3}$', '', colnames(df.cna)[-(1:3)]) #to format patient id to TCGA-XX-XXXX
  pids <- gsub('\\.','-',pids)
  colnames(df.cna) <- c('gene', 'locusid', 'cytoband', sub('(-[[:alnum:]]*){1}$', '', pids))
  
  #remove duplicates, so only to keep one unique patient
  toremoveidx <- getIdxDupPatientIDs(c('gene', 'locusid', 'cytoband', pids))
  if(length(toremoveidx) > 0){
    print(sprintf('There are %.0f duplicates due to tumor type and/or vial, removing them...', length(toremoveidx)))
    df.cna <- df.cna[, !(1:ncol(df.cna) %in% toremoveidx)] #ignore the first three columns, which are gene, locusid, cytoband
  }
  upids.cna <- sub('(-[[:alnum:]]*){1}$', '', pids) #unique patients with CNA
  
  #retrieve MAF
  fname <- paste0(cancertype,'_mut/')
  flist <- list.files(path=file.path(dir.db,fname), pattern="*.maf.txt")
  pids <- sub(".maf.txt",'',flist) #patient ids - tumortype
  toremoveidx <- getIdxDupPatientIDs(pids) #only keep one unique patient
  if(length(toremoveidx) >  0){
    print(sprintf('There are %.0f duplicates due to tumor type and/or vial, removing them...', length(toremoveidx)))
    pids <- pids[-toremoveidx]
  }
  upids.mut <- sub('(-[[:alnum:]]*){1}$', '', pids) #unique patients with mutation
  
  #loop through each patient and get the maf file
  df.mut <- data.frame()
  for(f in pids){
    patientid <- sub('(-[[:alnum:]]*){1}$', '', f) #patient id (remove the tumortype tag)
    
    #SNV data for this patient
    df.snv.p <- read.table(file.path(dir.db, fname, paste0(f, '.maf.txt')), header=T, stringsAsFactors=F, fill=T, sep='\t', quote='')
    
    df.mut.p <- df.snv.p[df.snv.p$Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins',
                                                                'Nonstop_Mutation', 'In_Frame_Del', 'Nonsense_Mutation',
                                                                'Translation_Start_Site', 'Splice_Site'),]
    if(nrow(df.mut.p) > 0){
      df.mut.p <- data.frame(patientid=patientid, gene=df.mut.p$Hugo_Symbol, stringsAsFactors = F)
      if(nrow(df.mut) < 1) {
        df.mut <- df.mut.p
      } else {
        df.mut <- rbind(df.mut, df.mut.p)
      }
    }
  }
  
  #tally patient numbers, by mutation, cna, or both
  mutcna.n.c <- data.frame(n.mut=length(upids.mut),
                           n.cna=length(upids.cna),
                           n.mutcna=length(union(upids.cna, upids.mut)),
                           cancertype=cancertype, stringsAsFactors=FALSE)
  
  if(nrow(df.mutcna.n) < 1) {
    df.mutcna.n <- mutcna.n.c
  } else {
    df.mutcna.n <- rbind(df.mutcna.n, mutcna.n.c)
  }
  
  #go through each patient and get the CNA or Mut info. The data has be retrieved above for all patients for this indication
  for(patientid in union(upids.cna, upids.mut)){
    #CNA data for this patient
    df.cna.p.all <- data.frame(genes=df.cna$gene, cna=0, stringsAsFactors=FALSE)
    if(patientid %in% colnames(df.cna))
      df.cna.p.all$cna <- df.cna[,colnames(df.cna)==patientid]
    df.cna.p <- df.cna.p.all[df.cna.p.all$cna < (-1*cna_ctrl_cutoff) | df.cna.p.all$cna > cna_ctrl_cutoff,] #only keep the genes that have CNA
    df.cna.p.del <- data.frame()
    df.cna.p.amp <- data.frame()
    if(nrow(df.cna.p) > 0){
      df.cna.p.del <- df.cna.p[df.cna.p$cna < (-1*cna_cutoff),] #loss of function
      df.cna.p.amp <- df.cna.p[df.cna.p$cna > cna_cutoff,] #gain of function, amplification
    }
    
    #mutations for this patient
    df.mut.p <- df.mut[df.mut$patientid==patientid,]
    
    #1. retrieve NMD genes (as likely loss of function) in this patient (mutations or loss)
    df.tmp <- data.frame()
    NMDgene.p.lof <- unique(c(df.mut.p$gene[df.mut.p$gene %in% NMD_genes],
                              df.cna.p.del$genes[df.cna.p.del$genes %in% NMD_genes]))
    
    if(length(NMDgene.p.lof) > 0){
      df.tmp <- data.frame(gene=NMDgene.p.lof, cancertype=cancertype, patientid=patientid, stringsAsFactors=FALSE)
      if(nrow(df.NMD.lof) < 1) {
        df.NMD.lof <- df.tmp
      } else {
        df.NMD.lof <- rbind(df.NMD.lof, df.tmp)
      }
    }
    
    #2. retrieve NMD genes (as likely gain of function) in this patient
    df.tmp <- data.frame()
    NMDgene.p.amp <- df.cna.p.amp$genes[df.cna.p.amp$genes %in% NMD_genes]
    if(length(NMDgene.p.amp) > 0){
      df.tmp <- data.frame(gene=NMDgene.p.amp, cancertype=cancertype, patientid=patientid, stringsAsFactors=FALSE)
      if(nrow(df.NMD.amp) < 1) {
        df.NMD.amp <- df.tmp
      } else {
        df.NMD.amp <- rbind(df.NMD.amp, df.tmp)
      }
    }
    
    df.tmp <- data.frame()
    NMDgene.p.del <- df.cna.p.del$genes[df.cna.p.del$genes %in% NMD_genes]
    if(length(NMDgene.p.del) > 0){
      df.tmp <- data.frame(gene=NMDgene.p.del, cancertype=cancertype, patientid=patientid, stringsAsFactors=FALSE)
      if(nrow(df.NMD.del) < 1) {
        df.NMD.del <- df.tmp
      } else {
        df.NMD.del <- rbind(df.NMD.del, df.tmp)
      }
    }
    
    #3. retrieve NMD genes (as likely gain of function and/or loss of function) for this patient
    df.tmp <- data.frame()
    NMDgene.p.any <- unique(c(df.mut.p$gene[df.mut.p$gene %in% NMD_genes],
                              df.cna.p$genes[df.cna.p$genes %in% NMD_genes]))
    
    if(length(NMDgene.p.any) > 0){
      df.tmp <- data.frame(gene=NMDgene.p.any, cancertype=cancertype, patientid=patientid, stringsAsFactors=FALSE)
      if(nrow(df.NMD.any) < 1) {
        df.NMD.any <- df.tmp
      } else {
        df.NMD.any <- rbind(df.NMD.any, df.tmp)
      }
    }
    
    #4. retrieve NMD genes, by mut, and by cna
    df.tmp <- data.frame()
    if(nrow(df.mut.p) > 0){
      NMDgene.p.mut <- unique(df.mut.p$gene[df.mut.p$gene %in% NMD_genes])
      
      if(length(NMDgene.p.mut) > 0){
        df.tmp <- data.frame(gene=NMDgene.p.mut, cancertype=cancertype, patientid=patientid, stringsAsFactors=FALSE)
        if(nrow(df.NMD.mut) < 1) { df.NMD.mut <- df.tmp } else { df.NMD.mut <- rbind(df.NMD.mut, df.tmp) }
      }
    }
    
    df.tmp <- data.frame()
    if(nrow(df.cna.p) > 0){
      NMDgene.p.cna <- unique(df.cna.p$genes[df.cna.p$genes %in% NMD_genes])
      
      if(length(NMDgene.p.cna) > 0){
        df.tmp <- data.frame(gene=NMDgene.p.cna, cancertype=cancertype, patientid=patientid, stringsAsFactors=FALSE)
        if(nrow(df.NMD.cna) < 1) { df.NMD.cna <- df.tmp } else { df.NMD.cna <- rbind(df.NMD.cna, df.tmp) }
      }
    }
  }
  
}

#Retrieve NMD metrics from previous analyses and merge
load('./outputs/18.0604/pancancer.rda')

nmdcolnames <- colnames(df.merged)[grep('nmd',colnames(df.merged))]
nmdcolnames <- nmdcolnames[!grepl('bool',nmdcolnames) & 
                             !grepl('_n_decayed',nmdcolnames) &
                             !grepl('_max',nmdcolnames)] #don't need the .bool, _n_decayed, _max metrics
df.nmd <- data.frame(patientid=row.names(df.merged), df.merged[nmdcolnames], stringsAsFactors=FALSE)

df.NMD.lof.all <- data.table::copy(df.NMD.lof)
df.NMD.amp.all <- data.table::copy(df.NMD.amp)
df.NMD.del.all <- data.table::copy(df.NMD.del)
df.NMD.any.all <- data.table::copy(df.NMD.any)
df.NMD.mut.all <- data.table::copy(df.NMD.mut)
rm(df.NMD.mut)
df.NMD.cna.all <- data.table::copy(df.NMD.cna)
rm(df.NMD.cna)

df.NMD.lof <- df.NMD.lof %>% select(cancertype, patientid) %>% unique() %>% left_join(df.nmd, by="patientid")
df.NMD.amp <- df.NMD.amp %>% select(cancertype, patientid) %>% unique() %>% left_join(df.nmd, by="patientid")
df.NMD.del <- df.NMD.del %>% select(cancertype, patientid) %>% unique() %>% left_join(df.nmd, by="patientid")
df.NMD.any <- df.NMD.any %>% select(cancertype, patientid) %>% unique() %>% left_join(df.nmd, by="patientid")

#get the control - no any NMD abberations (mutation nor CNA)
df.NMD.ctrl <- df.merged[!row.names(df.merged) %in% df.NMD.any$patientid,]
df.NMD.ctrl <- data.frame(patientid=row.names(df.NMD.ctrl),
                          cancertype=df.NMD.ctrl$cancertype,
                          df.NMD.ctrl[nmdcolnames], stringsAsFactors=FALSE)

#remove any patients where there are both loss and gain of NMD genes
patients.ambi <- c()
if(rm.gainlossboth){
  patients.ambi <- intersect(df.NMD.lof$patientid, df.NMD.amp$patientid)
  df.NMD.lof <- df.NMD.lof[!df.NMD.lof$patientid %in% patients.ambi,]
  df.NMD.amp <- df.NMD.amp[!df.NMD.amp$patientid %in% patients.ambi,]
  df.NMD.del <- df.NMD.del[!df.NMD.del$patientid %in% patients.ambi,]
}

#get patients with co-occurrence mut/amp
df.mutampN <- df.NMD.any.all %>% group_by(patientid) %>% mutate(n=n()) %>% select(cancertype,patientid,n) %>% as.data.frame()

#internal: double check the logics above and the ctrl, amp, and del results make sense
# ck.cna.ctrl <- df.cna[colnames(df.cna) %in% c('gene',df.NMD.ctrl$patientid)]
# ck.cna.ctrl <- ck.cna.ctrl[ck.cna.ctrl$gene %in% NMD_genes,] %>% select(-gene)
# sum(apply(ck.cna.ctrl,2,sum))
# 
# ck.cna.amp <- df.cna[colnames(df.cna) %in% c('gene',df.NMD.amp$patientid)]
# ck.cna.amp <- ck.cna.amp[ck.cna.amp$gene %in% NMD_genes,] %>% select(-gene)
# apply(ck.cna.amp,2,sum) > 0
# 
# ck.cna.del <- df.cna[colnames(df.cna) %in% c('gene',df.NMD.lof$patientid)]
# ck.cna.del <- ck.cna.del[ck.cna.del$gene %in% NMD_genes,] %>% select(-gene)
# apply(ck.cna.del,2,sum) <= 0

#------- analyses --------
df.NMD.mutlossamp <- unique(rbind(df.NMD.lof, df.NMD.amp))

#any alterations
toAnalyze <- list(NMD.any=df.NMD.any) #NMD with any alterations

#co-alterations
# toAnalyze <- list(NMD.amp=df.NMD.amp[df.NMD.amp$patientid %in% df.mutampN[df.mutampN$n>1,'patientid'],], #NMD with co-amp
#                   NMD.any.cooccur=df.NMD.any[df.NMD.any$patientid %in% df.mutampN[df.mutampN$n>1,'patientid'],]) #NMD with co-alt

NMD.ref.name <- 'NMD_ctrl'
df.NMD.ref <- df.NMD.ctrl #reference dataset to be used for the comparison, e.g. background control

# NMD.ref.name <- 'NMD_1variant'
# df.NMD.ref <- df.NMD.any[df.NMD.any$patientid %in% df.mutampN[df.mutampN$n<2,'patientid'],]

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
        nmdval.ctrl <- df.NMD.ref[nmdcolname]
      } else {
        nmdval <- df.NMD[df.NMD$cancertype==cancertype,nmdcolname]
        nmdval.ctrl <- df.NMD.ref[df.NMD.ref$cancertype==cancertype,nmdcolname]
      }
      
      if(!all(is.na(nmdval)) && length(nmdval) > 0 && !all(is.na(nmdval.ctrl)) && length(nmdval.ctrl) > 0){
        df1 <- data.frame(val=nmdval, type='NMD')
        colnames(df1) <- c('val','type')
        df2 <- data.frame(val=nmdval.ctrl, type=NMD.ref.name)
        colnames(df2) <- c('val','type')
        df <- rbind(df1, df2)
        
        val <- wilcox.test(df1$val, df2$val)
        diff <- median(df1$val, na.rm=TRUE) - median(df2$val, na.rm=TRUE)
        tmp <- data.frame(cancertype=cancertype, nmdmetric=nmdcolname, diff=diff, pval=val$p.value)
        if(nrow(nmdstat) < 1) { nmdstat <- tmp } else { nmdstat <- rbind(nmdstat, tmp) }
        
        r <- tryCatch({
          val <- wilcox.test(log10(df1$val), log10(df2$val))
          diff <- median(log10(df1$val), na.rm=TRUE) - median(log10(df2$val), na.rm=TRUE)
          tmp <- data.frame(cancertype=cancertype, nmdmetric=nmdcolname, diff=diff, pval=val$p.value)
          if(nrow(nmdstat.log) < 1) { nmdstat.log <- tmp } else { nmdstat.log <- rbind(nmdstat.log, tmp) }
        }, error=catch.errorfunc, finally={})
      }
    }
  }
  
  nmdstat$qval <- p.adjust(nmdstat$pval, method="BH")
  nmdstat.log$qval <- p.adjust(nmdstat.log$pval, method="BH")
  
  #Keep the statistically significant ones and plots those
  nmdstat.log.sig <- nmdstat.log[nmdstat.log$qval < qval_cutoff,]
  
  plotBarsNMD <- function(x){
    dirout <- paste0(dirout, '/sig', qval_cutoff,'_', anlyzname, '/')
    if(!file.exists(dirout)) dir.create(dirout)
    cancertype <- x['cancertype']
    nmdcolname <- x['nmdmetric']
    
    if(cancertype == 'pancancer'){
      nmdval <- df.NMD[nmdcolname]
      nmdval.ctrl <- df.NMD.ref[nmdcolname]
    } else {
      nmdval <- df.NMD[df.NMD$cancertype==cancertype,nmdcolname]
      nmdval.ctrl <- df.NMD.ref[df.NMD.ref$cancertype==cancertype,nmdcolname]
    }
    if(length(nmdval) == 1 && is.na(nmdval)) return()
    if(length(nmdval.ctrl) == 1 && is.na(nmdval.ctrl)) return()
    if(length(nmdval) > 0 && length(nmdval.ctrl) > 0){
      df1 <- data.frame(val=nmdval, type=anlyzname)
      colnames(df1) <- c('val','type')
      df2 <- data.frame(val=nmdval.ctrl, type=NMD.ref.name)
      colnames(df2) <- c('val','type')
      df <- rbind(df1, df2)
      
      r <- tryCatch({
        ggplot(df, aes(as.factor(df$type), val)) + geom_boxplot() + geom_jitter(width=0.2) +
          xlab(cancertype) + ylab(nmdcolname) + scale_y_log10() +
          stat_compare_means(method='wilcox.test',comparisons=list(c(anlyzname,NMD.ref.name)))
        ggsave(file.path(dirout, paste0('nmdlog10_', cancertype,'_',nmdcolname,'.pdf')))
        #note since using scale transform and not coord transform, the statistics for stat_compare is test on the
        #transformed data and not the raw data values
      }, error=catch.errorfunc, finally={})
    }
  }
  
  if(nrow(nmdstat.log.sig)>0) apply(nmdstat.log.sig, 1, plotBarsNMD)
  write.csv(nmdstat.log, file=paste0(dirout, '/nmdstatlog_', anlyzname,'.csv'))
}

#analyze the median difference (NMD - NMD ctrl) for all the indication_metric, colored by significance
derivesign <- function(sig,diff){
  #derive the sign (directionaity) of the delta (NMD - NMD ctrl), marked only for ones that are statistically significant
  #+1 NMD-NMD ctrl > 0; -1 NMD-NMD ctrl < 0
  #0 ambiguous
  #NA not statistically significant
  if(all(sig==FALSE))
    return(NA)
  
  if(all(diff[sig==TRUE]>0))
    return(1)
  else if(all(diff[sig==TRUE]<0))
    return(-1)
  else
    return(0)
}

df.sigsign <- data.frame()
nmdstat.pfx <- 'nmdstatlog'
for(anlyzname in names(toAnalyze)){
  if(nrow(toAnalyze[[anlyzname]]) < 1) next
  nmdstat <- read.csv(file=paste0(dirout, '/', nmdstat.pfx, '_', anlyzname,'.csv'))
  nmdstat$name <- paste0(nmdstat$cancertype, nmdstat$nmdmetric)
  nmdstat$sig <- nmdstat$qval < qval_cutoff
  
  #generate the statistically significant results, marked with directionality of the NMD-NMD_ref
  df <- nmdstat %>% group_by(cancertype) %>% mutate(sigsign=derivesign(sig,diff)) %>% select(cancertype, sigsign) %>% unique() %>% as.data.frame()
  if(nrow(df.sigsign) < 1)
    df.sigsign <- df
  else
    df.sigsign <- df.sigsign %>% left_join(df, by="cancertype")
  
  #box plot of all indications
  ggplot(nmdstat, aes(x=name, y=diff, color=sig)) + geom_point(size=2) + geom_hline(yintercept = 0) +
    labs(x="Indication_NMDmetric", y="delta median (NMD - NMD ctrl)", color=paste0("qval < ",qval_cutoff)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          text = element_text(size=5)) +
    scale_color_manual(values=c("gray","red"))
  ggsave(file.path(dirout, paste0('plot_', nmdstat.pfx, '_',anlyzname,'.pdf')), width=12, height=6)
  
  #box plot of only statistically significant indications
  ggplot(nmdstat[nmdstat$cancertype %in% df[(df$sigsign==1 | df$sigsign==-1) & !is.na(df$sigsign), 'cancertype'],], 
         aes(x=name, y=diff, color=sig)) + geom_point(size=3) + geom_hline(yintercept = 0) +
    labs(x="Indication_NMDmetric", y="delta median (NMD - NMD ctrl)", color=paste0("qval < ",qval_cutoff)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          text = element_text(size=5)) +
    scale_color_manual(values=c("gray","red"))
  ggsave(file.path(dirout, paste0('plot_', nmdstat.pfx, '_',anlyzname,'_sigonly.pdf')), width=12, height=6)
  
}

#heatmap summary
# colnames(df.sigsign) <- c('cancertype', names(toAnalyze))
# df.sigsign[df.sigsign==0] <- 0.2 #stat significant, but ambiguous (>0 and <0 in the median delta)
# df.sigsign[is.na(df.sigsign)] <- -0.2 #not stat significant
# df.tmp <- df.sigsign %>% select(-cancertype)
# row.names(df.tmp) <- df.sigsign$cancertype
# palettecol <- colorRampPalette(c("#4682B4", '#F5F5F5', "grey", "#DC143C"))(n = 4)
# pdf(file.path(dirout, 'heatmap_summary.pdf'), width=12, height=12)
# heatmap.2(as.matrix(df.tmp), 
#           col=palettecol, scale="none",
#           key.xlab=sprintf('Statistically significant, signed\n(+1=sig, NMD-%s>0; -1=sig, NMD-%s<0)',NMD.ref.name,NMD.ref.name), key.title=F,
#           key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,
#           margins = c(8,8),
#           sepwidth=c(0.01, 0.01), sepcolor='white',
#           colsep=1:ncol(df.sigsign), rowsep=1:nrow(df.sigsign))
# dev.off()

#proportion of amp and lof, in each indication
# a <- df.NMD.amp %>% group_by(cancertype) %>% count()
# b <- df.NMD.lof %>% group_by(cancertype) %>% count()
# c <- a %>% left_join(b, by="cancertype")
# colnames(c) <- c('cancertype', 'amp', 'lof')
# c$cancertype <- factor(c$cancertype, levels = c$cancertype[order(c$amp/(c$amp+c$lof), decreasing=TRUE)]) #sort by percent amp
# #c$cancertype <- factor(c$cancertype, levels = c$cancertype[order(c$amp+c$lof, decreasing=TRUE)]) #sort by total
# df1 <- melt(c, id.var='cancertype')
# ggplot(df1, aes(x=cancertype, y=value, fill=variable)) + geom_bar(stat = "identity")
# ggsave(file.path(dirout, paste0('plot_amplof.pdf')), width=12, height=6)

#proportion of amp and del, in each indicationdel
a <- df.NMD.amp %>% group_by(cancertype) %>% count()
b <- df.NMD.del %>% group_by(cancertype) %>% count()
c <- a %>% left_join(b, by="cancertype")
colnames(c) <- c('cancertype', 'amp', 'del')
c$cancertype <- factor(c$cancertype, levels = c$cancertype[order(c$amp/(c$amp+c$del), decreasing=TRUE)]) #sort by percent amp
#c$cancertype <- factor(c$cancertype, levels = c$cancertype[order(c$amp+c$lof, decreasing=TRUE)]) #sort by total
df1 <- melt(c, id.var='cancertype')
ggplot(df1, aes(x=cancertype, y=value, fill=variable)) + geom_bar(stat = "identity")
ggsave(file.path(dirout, paste0('plot_ampdel.pdf')), width=12, height=6)

#save
save(NMD_genes, cna_cutoff, cna_ctrl_cutoff, qval_cutoff, rm.gainlossboth, #settings
     patients.ambi, nmdcolnames, df.mutampN,
     df.NMD.lof, df.NMD.amp, df.NMD.any, df.NMD.ctrl, df.mutcna.n,
     df.NMD.lof.all, df.NMD.amp.all, df.NMD.any.all, df.NMD.mut.all, df.NMD.cna.all,
     file=file.path(dirout, 'results.rda'))


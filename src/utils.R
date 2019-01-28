#utils: helper methods
#------- general helper methods --------
print.seg <- function(txt){
  print("###############################################")
  print(sprintf("## %s", txt))
  print("###############################################")
}

print.seg2 <- function(txt){
  print(sprintf("##---------- %s ----------##", txt))
}

#define some colors
# colorRampPalette(brewer.pal(9, "RdBu"))(100)
# colorRampPalette(c('blue', 'yellow'))(12)

catch.errorfunc <- function(e){
  print(e)
}

calcGeoMean <- function(vals){
  #return(prod(vals)^(1/length(vals)))
  return(exp(mean(log(vals))))
}

fillNAmax <- function(x){
  #fill NAs with the maximum val
  x <- as.numeric(x)
  xmax <- max(x, na.rm=TRUE)
  x[is.na(x)] <- xmax
  return(x)
}

rmNonExtVars <- function(varnames, df){
  #check if the varnames exist in the df colnames, if not remove, return the varnames with only ones exist in df
  if(sum(!(names(varnames) %in% colnames(df))) > 0){
    print("WARNING: some var names do not exist in dataset column names provided...")
    c2exclude <- !(names(varnames) %in% colnames(df))
    print(paste0("Vars not found and to exclude: ", paste(names(varnames)[c2exclude], collapse=',')))
    varnames <- varnames[!c2exclude]
  }
  return(varnames)
}

rmConstVars <- function(df){
  #remove features with no variationss
  varnames <- names(df)
  for(varname in varnames){
    if(nrow(unique(df[varname])) < 2){
      print(sprintf("Removing %s, due to no variation in values for this variable...", varname))
      df <- df[!(colnames(df) %in% c(varname))]
    }
  }
  return(df)
}

getIdxDupPatientIDs <- function(pids){
  #retrieve the index to remove because they are duplicate based on the patient id
  #keep only the lowest tumor type (<10) 
  pids_only <-  sub('(-[[:alnum:]]*){1}$', '', pids) #just patient id
  pids_tumortypevials <-  sub('^([[:alnum:]]*-){3}', '', pids) #just patient ids
  pids_vial <- sub('^([[:alnum:]]){2}', '', pids_tumortypevials) #just vial number
  pids_tumortype <- as.numeric(sub('\\w$', '', pids_tumortypevials)) #just tumor type numbers
  df <- data.frame(id=pids, pid=pids_only, tumortype=pids_tumortype, vial=pids_vial)
  toremovepidx <- c()
  
  #remove duplicates based on patient ids
  if(sum(duplicated(df$pid)) > 0){
    #for ones with duplicate patient ids, just keep the one with the lowest tumor type number, and of just lowest vial
    dupids <- df$pid[duplicated(df$pid)]
    for(dupid in dupids){
      df.tmp <- df[df$pid==dupid,]
      df.tmp <- df.tmp[!duplicated(df.tmp$id),] #remove duplciated ids
      df.tmp.sorted <- df.tmp[order(df.tmp$tumortype, df.tmp$vial, decreasing=F),]
      toremoveid <- df.tmp.sorted[-1,'id']
      toremovepidx <- c(toremovepidx, which(pids %in% toremoveid))
    }
  }
  
  #remove any non-tumor samples (e.g. normal, control are removed)
  toremovepidx <- c(toremovepidx, which(df$tumortype>=10)) 
  
  #remove duplicates based on patientid-tumortype-vial
  toremovepidx <- c(toremovepidx, which(duplicated(pids)))
  
  toremovepidx <- unique(toremovepidx)
}

handleMissingVals <- function(df.merged, dropsmall=TRUE){
  #check for missing values in the dataset, if less than 10 of the samples are missing, drop them (if dropsmall is TRUE)
  missN <- sum(!complete.cases(dplyr::select(df.merged, -cyt.bool)))
  if(missN>0){
    print(sprintf('WARNING: there are %.0f entries with missing data... ',missN))
    if(dropsmall && missN<10){
      print('NOTE: missing data is <10, will use complete cases only')
      df.merged <- df.merged[complete.cases(dplyr::select(df.merged, -cyt.bool)),]
    }
  }
  
  return(df.merged)
}

cleanup.factors <- function(df, colvars=NULL){
  #remove empty levels, give warning if some factor has more than 20 levels
  #if colvars, will get all columns that are factors in given df
  
  if(is.null(colvars)){
    colvars <- colnames(df)[sapply(df,is.factor)]
  }
  
  for(varcat in colvars){
    levcount <- length(levels(df[[varcat]]))
    if(levcount>20){
      print(sprintf("WARNING: variable %s has %.0f levels!", varcat, levcount))
    }
    
    if(levcount > length(unique(df[[varcat]]))){
      print(sprintf("There exists empty levels for variable %s, dropping the empty levels...", varcat))
      df[varcat] <- droplevels(df[[varcat]])
    }
  }
  
  return(df)
}

getdfbool <- function(df.merged){
  #retrieve a reduced df.merged, only consist of cyt.cat of high or low (omit the rest)
  df.merged.bool <- df.merged[!is.na(df.merged$cyt.bool),]
  df.merged.bool$cyt.cat <- droplevels(df.merged.bool$cyt.cat)
  return(df.merged.bool)
}

getdfbool.rm <- function(df.merged){
  #retrieve a slice of dataset, only with high/low CYT
  print('Retrieving dataset slice with just high/low CYT...')
  
  df.merged.bool <- df.merged[!is.na(df.merged$cyt.bool),]
  df.merged.bool <- rmConstVars(df.merged.bool)
  df.merged.bool <- cleanup.factors(df.merged.bool)
  return(df.merged.bool)
}

#------- model helper methods --------
buildModel <- function(mod.fml, mod.data, mod.type, mod.method, params.txt=NULL){
  #Inputs:
  # mod.fml: formula as text
  # mod.data: input df
  # mod.type: regression or classification
  # params.txt: txt of hyperparameter(s) for the model, this can be the optimal parameter based on cross-validation
  mod.data.tmp <<- mod.data #make it global
  
  #linear regression
  if(mod.method == 'lm')
    modtxt <- paste0('mod <- lm(', mod.fml, ', data=mod.data.tmp)')
  
  #logistic regression
  else if(mod.method == 'logit')
    modtxt <- paste0("mod <- glm(", mod.fml, ", data=mod.data.tmp, family=binomial('logit'))") 
  
  #random forest
  else if(mod.method == 'rf'){ 
    if(is.null(params.txt)) params.txt <- 'ntree=1000, mtry=2'
    #modtxt <- paste0("mod <- randomForest(", mod.fml, ", data=mod.data.tmp, ", params.txt,", importance=TRUE)") 
    modtxt <- paste0("mod <- rfPermute(", mod.fml, ", data=mod.data.tmp, ", params.txt,", nrep=50)") 
    
    #logistic regression, regularized - not fully implemented
  } else if(mod.method == 'logit.reg'){ 
    if(is.null(params.txt)) params.txt <- 'alpha=1, nlambda = 100'
    #NOTE: format is more annoying for glmnet, does not accept formula
    #modtxt <- paste0("mod <- glmnet(", mod.fml, ", data=mod.data.tmp, family='binomial', ",params.txt,")") #logistic regressions
  }
  
  eval(parse(text=modtxt))
  
  return(mod)
}

buildModel.cv <- function(mod.fml, mod.data, mod.type, mod.method){
  #Inputs:
  # mod.fml: formula as text
  # mod.data: input df
  # mod.type: regression or classification
  mod.data.tmp <<- mod.data #make it global
  
  mod.metric <- ifelse(mod.type=='reg', 'RMSE', 'ROC')
  params.best <- NULL
  
  #Linear regression
  if(mod.method == 'lm'){
    train_ctrl <- trainControl(method="cv", number=10, savePredictions=TRUE) #10-fold CV
    mod.cv <- train(eval(parse(text=mod.fml)), data=mod.data.tmp, trControl=train_ctrl, method='lm',
                    #preProcess=c("center","scale"),
                    metric=mod.metric)
    
    #Logistic regresion
  } else if(mod.method == 'glm.logistic'){
    train_ctrl <- trainControl(method="cv", number=10, savePredictions=TRUE, #10-fold CV
                               classProbs=TRUE, #estimate class probabilities
                               summaryFunction=twoClassSummary) #evaluate using AUC
    mod.cv <- train(eval(parse(text=mod.fml)), data=mod.data.tmp, trControl=train_ctrl, method='glm', family=binomial(link='logit'),
                    #preProcess=c("center","scale"),
                    metric=mod.metric)
    
    #Random forest
  } else if(mod.method=='rf'){
    #train with cross-validation
    if(mod.type == 'bin')
      train_ctrl <- trainControl(method="cv", number=10, savePredictions=TRUE, #10-fold CV
                                 classProbs=TRUE, #estimate class probabilities
                                 summaryFunction=twoClassSummary) #evaluate using AUC
    else if(mod.type == 'reg')
      train_ctrl <- trainControl(method="cv", number=10, savePredictions=TRUE)
    
    #train with out-of-bag error
    # if(mod.type == 'bin'){
    #   train_ctrl <- trainControl(method="oob", savePredictions=TRUE, #10-fold CV
    #                              classProbs=TRUE) #estimate class probabilities
    #   mod.metric <- 'Accuracy' #for OOB, ROC is not available, use Accuracy as metric
    # } else if(mod.type == 'reg')
    #   train_ctrl <- trainControl(method="oob", number=10, savePredictions=TRUE)
    
    
    mtry.val<- round(sqrt(ncol(mod.data.tmp) - 1))
    tunegrid <- expand.grid(.mtry=expand.grid(mtry=c(round(mtry.val/2), mtry.val, 2*mtry.val) )) #mtry parameter values to try
    mod.cv.l <- NULL
    mod.cv.rdf <- data.frame()
    ntrees.grid <- c(500, 1000) #ntree parameter values to try
    n=1
    for(ntree in ntrees.grid){
      mod.cv <- train(eval(parse(text=mod.fml)), data=mod.data.tmp, trControl=train_ctrl, method='rf', 
                      #preProcess=c("center","scale"),
                      tuneGrid=tunegrid, ntree=ntree,
                      importance=TRUE,
                      metric=mod.metric)
      mod.cv.l[[n]] <-mod.cv
      mod.cv$results$ntree <- ntree
      mod.cv.rdf <- rbind(mod.cv.rdf, mod.cv$results)
      n<-n+1
    }
    
    #choose the best hyperparamter
    if(mod.type == 'bin'){
      ntree.best <- mod.cv.rdf[which.max(mod.cv.rdf$ROC),'ntree']
      mtry.best <- mod.cv.rdf[which.max(mod.cv.rdf$ROC),'mtry']
    } else if(mod.type == 'reg'){
      ntree.best <- mod.cv.rdf[which.min(mod.cv.rdf$RMSE),'ntree']
      mtry.best <- mod.cv.rdf[which.min(mod.cv.rdf$RMSE),'mtry']
    }
    mod.cv <- mod.cv.l[[which(ntree.best %in% ntrees.grid)]] #rewrite mod.cv with the one from the best tuned parameters
    params.best <- sprintf('ntree=%.0f, mtry=%.0f', ntree.best, mtry.best)
    
    
    #Regularized GLM
  } else if(mod.method=='glmnet'){
    if(mod.type == 'bin'){
      train_ctrl <- trainControl(method="cv", number=10, savePredictions=TRUE, #10-fold CV
                                 classProbs=TRUE, #estimate class probabilities
                                 summaryFunction=twoClassSummary) #evaluate using AUC
      mod.cv <- train(eval(parse(text=mod.fml)), data=mod.data.tmp, trControl=train_ctrl, method='glmnet',
                      family='binomial',
                      #preProcess=c("center","scale"),
                      metric="ROC")
    } else if(mod.type == 'reg'){
      train_ctrl <- trainControl(method="cv", number=10, savePredictions=TRUE) #evaluate using AUC
      mod.cv <- train(eval(parse(text=mod.fml)), data=mod.data.tmp, trControl=train_ctrl, method='glmnet',
                      #preProcess=c("center","scale"),
                      metric="RMSE")
    }
  }
  
  return(list(mod.cv, params.best))
}

summarizeModel <- function(mod, mod.cv, yname, mod.data.tval, mod.data.test, mod.name, mod.type, mod.method, dirout=dir.out){
  if(mod.type=='reg')
    summarizeModel.reg(mod, mod.cv, yname, mod.data.tval, mod.data.test, mod.name, mod.method, dirout)
  else if(mod.type=='bin')
    summarizeModel.bin(mod, mod.cv, yname, mod.data.tval, mod.data.test, mod.name, mod.method, dirout)
}

summarizeModel.reg <- function(mod, mod.cv, yname, mod.data.tval, mod.data.test, mod.name, mod.method, dirout=dir.out){
  #summarize the model
  #uses both model (mod) and cross-validated model (mod.cv) object (from train())
  #returns summary statistics
  
  #print summary
  print(summary(mod))
  
  #cross validation
  print('Cross validation...')
  print(mod.cv)
  print(mod.cv$results)
  
  #variable importance
  pdf(file.path(dirout, paste0(mod.name,'_varImp_rel.pdf')), height=7, width=7)
  print(plot(varImp(mod.cv)))
  dev.off()
  
  #method specific plots/analyses
  if(mod.method == 'rf'){
    #errors during training
    pdf(file.path(dirout, paste0(mod.name,'_errors.pdf')))
    plot(mod.cv$finalModel)
    dev.off()
    
    #null distribution
    # varsn <- length(attr(mod$terms,'term.labels'))
    # layout(matrix(1:(ceil(varsn/2)*2), ncol = 2))
    # plotNull(mod.rfp) 
    # layout(matrix(1))
    
    #importance, with signifiance, from rfPermute
    pdf(file.path(dirout, paste0(mod.name,'_varImp_sig.pdf')), height=11, width=7)
    rpimpvals.unscaled <- rp.importance(mod, scale=FALSE)
    plot(rpimpvals.unscaled, alpha=0.05)
    dev.off()
    
    pdf(file.path(dirout, paste0(mod.name,'_varImp_sig_scaled.pdf')), height=11, width=7)
    rpimpvals.scaled <- rp.importance(mod, scale=TRUE)
    plot(rpimpvals.scaled, alpha=0.05)
    dev.off()
    
    #importance, simply from importance(), which is the scaled version
    pdf(file.path(dirout, paste0(mod.name,'_varImp.pdf')), height=7, width=11)
    varImpPlot(mod.cv$finalModel)
    dev.off()
    
    #write tables
    write.csv(rpimpvals.unscaled, file=file.path(dirout, paste0(mod.name,'_RFimp_unscaled.csv'))) #unscaled
    write.csv(rpimpvals.scaled, file=file.path(dirout, paste0(mod.name,'_RFimp_scaled.csv'))) #scaled
  } else {
    pdf(file.path(dirout, paste0(mod.name,'_diagnostics.pdf')), height=7, width=7)
    layout(matrix(c(1,2,3,4),2,2))
    plot(mod)
    dev.off()
  }
  
  #summaries
  summarizeModel.reg.s <- function(a.mod, a.data, a.name, passPredict=TRUE){
    if(passPredict)
      yhat = predict(a.mod, subset(a.data, select=-c(cyt,cyt.bool,cyt.cat)))
    else
      yhat = predict(a.mod)
    
    yobs = a.data[[yname]]
    
    #pred vs obs
    pdf(file.path(dirout, paste0(mod.name,'_predvsobs.pdf')), height=7, width=7)
    plot(yhat, yobs)
    abline(0,1)
    dev.off()
    
    r2 <- round(1 - (sum((yhat-yobs)^2)/sum((yobs-mean(yobs))^2)),2) #calculate R-squared
    rmse <- sqrt(mean((yhat-yobs)^2))
    mae <- mean(abs(yhat-yobs))
    return(r2)
  }
  
  #retrieve metric for training, validation, test
  r2.tval <- summarizeModel.reg.s(mod.cv$finalModel, mod.data.tval, 'train.val_set', passPredict=FALSE)
  if(dbset.trainval=='all' && dbset.test=='all')
    r2.test <- summarizeModel.reg.s(mod.cv$finalModel, mod.data.test, 'test_set', passPredict=FALSE)
  else
    r2.test <- summarizeModel.reg.s(mod.cv$finalModel, mod.data.test, 'test_set')
  
  if(nrow(mod.cv$results) > 1){ #there was parameter tuning, select the optimized one
    #TODO: for now just deal with random forest
    r <- mod.cv$results[mod.cv$results$mtry==mod.cv$bestTune$mtry,]
    r2.cv <- round(r$Rsquared,2)
    r2.cv.sd <- round(r$RsquaredSD,2)
  } else {
    r2.cv <- round(mod.cv$results$Rsquared,2)
    r2.cv.sd <- round(mod.cv$results$RsquaredSD,2)
  }
  
  #return stats
  f<-as.character(mod.cv$terms)
  stats.vals <- list(fml=gsub(' ','',paste0(f[2],f[1],f[3])),
                     r2.tval=r2.tval, #r2 (training)
                     r2.cv=r2.cv, #r2 (validation), averaged
                     r2.cv.sd=r2.cv.sd, #r2 (validation)
                     r2.test=r2.test) #r2 (test)
  
  return(stats.vals)
}

summarizeModel.bin <- function(mod, mod.cv, yname, mod.data.tval, mod.data.test, mod.name, mod.method, dirout=dir.out){
  #summarize the model
  #uses both model (mod) and cross-validated model (mod.cv) object (from train())
  #returns summary statistics
  
  #print summary
  print(summary(mod))
  
  #cross validation
  print('Cross validation...')
  print(mod.cv)
  print(mod.cv$results)
  
  #variable importance
  pdf(file.path(dirout, paste0(mod.name,'_varImp_rel.pdf')), height=13, width=7)
  print(plot(varImp(mod.cv)))
  dev.off()
  
  #summaries
  summarizeModel.bin.s <- function(a.mod, a.data, a.name, passPredict=TRUE){
    if(mod.method == 'rf'){
      #predict
      if(passPredict)
        p.test <- predict(a.mod, subset(a.data,select=-c(cyt,cyt.bool,cyt.cat)), type='response')
      else
        p.test <- predict(a.mod, type='response')
      
      #confusion matrix
      print(caret::confusionMatrix(p.test, a.data[[yname]]))
      
      #ROC curve
      if(passPredict)
        p.test <- predict(a.mod, subset(a.data,select=-c(cyt,cyt.bool,cyt.cat)), type='prob')
      else
        p.test <- predict(a.mod, type='prob')
      pred <- prediction(p.test[,'high'], a.data$cyt.bool)
      perf <- performance(pred, "tpr", "fpr")
      perf.auc <- performance(pred, "auc")
      auc.val <- round(perf.auc@y.values[[1]],2)
      
      pdf(file.path(dirout, paste0(mod.name,'_',a.name,'.pdf')))
      plot(perf, col=1, lwd=2, main=paste0(mod.name,'; AUC: ', auc.val, '; ', a.name))
      abline(0,1,col='gray')
      dev.off()
    } else if(mod.method == 'glm.logistic'){
      #mod.cv
      p.test <- predict(a.mod, subset(a.data,select=-c(cyt,cyt.bool,cyt.cat)), type='response')
      
      #response is probability need to convert to factor
      print(caret::confusionMatrix(as.numeric(p.test>0.5), as.numeric(a.data[[yname]])-1))
      
      #ROC curve
      pred <- prediction(p.test, a.data$cyt.bool)
      perf <- performance(pred, "tpr", "fpr")
      perf.auc <- performance(pred, "auc")
      auc.val <- round(perf.auc@y.values[[1]],2)
      
      pdf(file.path(dirout, paste0(mod.name,'_',a.name,'.pdf')))
      plot(perf, col=1, lwd=2, main=paste0(mod.name,'; AUC: ', auc.val, '; ', a.name))
      abline(0,1,col='gray')
      dev.off()
    }
    
    return(auc.val)
  }
  
  #method specific plots/analyses
  #retrieve metric for training and test
  if(mod.method == 'rf'){
    #errors during training
    pdf(file.path(dirout, paste0(mod.name,'_errors.pdf')))
    plot(mod.cv$finalModel)
    dev.off()
    
    if(nrow(mod$importance) > 1){
      #importance, with signifiance 
      pdf(file.path(dirout, paste0(mod.name,'_varImp_sig.pdf')), height=16, width=11)
      rpimpvals.unscaled <- rp.importance(mod, scale=FALSE)
      plot(rpimpvals.unscaled, alpha = 0.05)
      dev.off()
      
      pdf(file.path(dirout, paste0(mod.name,'_varImp_sig_scaled.pdf')), height=16, width=11)
      rpimpvals.scaled <- rp.importance(mod, scale=TRUE)
      plot(rpimpvals.scaled, alpha = 0.05)
      dev.off()
      
      #write tables
      write.csv(rpimpvals.unscaled, file=file.path(dirout, paste0(mod.name,'_RFimp_unscaled.csv'))) #unscaled
      write.csv(rpimpvals.scaled, file=file.path(dirout, paste0(mod.name,'_RFimp_scaled.csv'))) #scaled
      
      #importance
      pdf(file.path(dirout, paste0(mod.name,'_varImp.pdf')), height=16, width=12)
      varImpPlot(mod.cv$finalModel)
      dev.off()
    }
    
    auc.tval <- summarizeModel.bin.s(mod.cv$finalModel, mod.data.tval, 'train.val_set', passPredict=FALSE)
    if(dbset.trainval=='all' && dbset.test=='all')
      auc.test <- summarizeModel.bin.s(mod.cv$finalModel, mod.data.test, 'test_set', passPredict=FALSE)
    else
      auc.test <- summarizeModel.bin.s(mod.cv$finalModel, mod.data.test, 'test_set')
    
  } else if(mod.method == 'glm.logistic'){
    #for glm.logistic regression, there is error when using mod.cv$finalModel, so will use mod
    #goal is the same, use model built based on mod.data.tval, get train and test error
    auc.tval <- summarizeModel.bin.s(mod, mod.data.tval, 'train.val_set')
    auc.test <- summarizeModel.bin.s(mod, mod.data.test, 'test_set')
  }
  
  #retrieve metric for validation
  if(nrow(mod.cv$results) > 1){ #there was parameter tuning, select the optimized one
    #TODO: for now just deal with random forest
    r <- mod.cv$results[mod.cv$results$mtry==mod.cv$bestTune$mtry,]
    auc.cv <- round(r$ROC,2)
    auc.cv.sd <- round(r$ROCSD,2)
  } else {
    auc.cv <- round(mod.cv$results$ROC,2)
    auc.cv.sd <- round(mod.cv$results$ROCSD,2)
  }
  
  # qplot(perf@x.values[[1]], perf@y.values[[1]]) + geom_abline(color='gray') + geom_line(alpha=0.7) +
  #   xlab(perf@x.name) + ylab(perf@y.name) + ggtitle(paste0('ROC AUC: ', round(perf.auc@y.values[[1]],2))) +
  #   theme(plot.title=element_text(hjust=0.5))
  # ggsave(paste0(dir.out, mod.name,'.pdf'))
  
  #return stats
  f<-as.character(mod.cv$terms)
  stats.vals <- list(fml=gsub(' ','',paste0(f[2],f[1],f[3])),
                     auc.tval=auc.tval, #AUC (training)
                     auc.cv=auc.cv, #AUC (validation), averaged
                     auc.cv.sd=auc.cv.sd, #AUC (validation), standard deviation
                     auc.test=auc.test) #AUC (test)
  
  #method specific stats
  if(mod.method == 'rf'){
    err.obb <- round(tail(mod.cv$finalModel$err.rate,1)[,'OOB'],2)
    stats.vals <- c(stats.vals,
                    oob=err.obb) #
  }
  
  return(stats.vals)
}

assessModel <- function(varnames, yname, mod.data, mod.name, mod.type, mod.method, mod.rd, dirout=dir.out){
  #Builds and analyses (including cross-validation) the given full model, in addition,
  #a reduced model by VIP recursive elimination and AIC stepwise selection is produced and returned
  #Input:
  # varnames: vector of variable names
  # yname: y variable name
  # mod.data: input data frame
  # mod.name: model name
  # mod.type: model type, 'reg' or 'gen'
  # mod.method: model method
  # mod.rd: boolean, if to also generate reduced model
  # dirout: output directory
  #Return: list struct, with stats for full and reduced model
  #Requires: build.modRd
  
  #-- model define/evaluate
  varnames <- rmNonExtVars(varnames, mod.data) #check to make sure the varnames exist in df
  yname <- rmNonExtVars(yname, mod.data) #check to make sure the varnames exist in df
  if(length(varnames)<1){warning('X vars is empty and/or do not exist in df...'); return(NA)}
  if(length(yname)<1){warning('Y var does not exist in df...'); return(NA)}
  mod.fml <- paste0(yname, '~', paste0(varnames, collapse="+")) #formula
  
  #split training/validation/test
  #since using the caret package, the training/validation will be split within train func, just define separate test
  #80%/20% (test)
  tval.idx <- createDataPartition(mod.data[[yname]], p=0.80, list=FALSE)
  if(dbset.trainval=='train.val') mod.data.tval <- mod.data[tval.idx,]
  else if(dbset.trainval=='all') mod.data.tval <- mod.data
  
  if(dbset.test=='test') mod.data.test <- mod.data[-tval.idx,]
  else if(dbset.test=='all') mod.data.test <- mod.data
  
  #build the models
  r <- buildModel.cv(mod.fml, mod.data.tval, mod.type, mod.method) #cross-validation, and hyperparamter tuning, if applicable
  mod.cv <- r[[1]]
  params.best.txt <- r[[2]]
  mod <- buildModel(mod.fml, mod.data.tval, mod.type, mod.method, params.best.txt)
  
  print(sprintf('******** Assessing model | %s |  %s ******** ', mod.name, gsub(' ','',as.character(mod$call)[2]) ))
  
  #stepwise selection based on AIC
  r <- tryCatch({
    print("Stepwise AIC model selection...")
    res.step <- step(mod, direction="both")
    print(res.step)
    print("Final model selected by stepwise (both) selection by AIC:")
    print(summary(res.step)) #summarize final model
  }, error=catch.errorfunc, finally={})
  
  #-- model reduction with VIP recursive elimination, and stepwise AIC selection
  if(mod.rd){
    #VIF for multicollinearity
    #recursive elimination of multicollinear variables
    varnames.rvip <- varnames
    r <- tryCatch({ #need try/catch b/c error will be thrown if the X'X is singular, perfect multicollinearity
      print('Calculating variance inflation factors...')
      
      vifdf <- as.data.frame(car::vif(mod))
      vifdf <- data.frame(var=row.names(vifdf), VIF=vifdf[[1]])
      
      vifdf.all <- vifdf
      mod.tmp <- mod
      while(any(vifdf$VIF>5)){
        varmax.rm <- vifdf[which(vifdf$VIF==max(vifdf$VIF)), 'var']
        varkeep <- vifdf[!(vifdf$var %in% varmax.rm), 'var']
        mod.tmp.fml <- paste0(yname, '~', paste0(varkeep, collapse="+"))
        mod.tmp <- buildModel(mod.tmp.fml, mod.data.tval, mod.type, params.best.txt)
        vifdf <- as.data.frame(car::vif(mod.tmp))
        vifdf <- data.frame(var=row.names(vifdf), VIF=vifdf[[1]])
        vifdf.all <- merge(vifdf.all, vifdf, by='var', all=T)
      }
      colnames(vifdf.all)[-1] <- paste0('VIF_',1:(ncol(vifdf.all)-1))
      write.csv(vifdf.all, file=file.path(dirout, paste0(mod.name, '_vif.csv')))
      
      print('Final model selected by recusive VIF:')
      print(summary(mod.tmp))
      varnames.rvip <- varkeep #non-multicollinear variables
    }, error=function(e){print(e);print('WARNING: error in VIP recursive elimination, keeping all variables...')}, finally={})
    
    #stepwise selection based on AIC, of only variables non-multicollinear, from VIP recursive elimination
    if(length(varnames.rvip) > 1){
      print('Stepwise AIC model selection using non-multicollinear variables, selected based on VIP recursive elimination...')
      mod.rvip.fml <- paste0(yname, '~', paste0(varnames.rvip, collapse="+"))
      mod.rvip <- buildModel(mod.rvip.fml, mod.data.tval, mod.type, mod.method, params.best.txt)
      
      r <- tryCatch({
        mod.rvip.raic <- step(mod.rvip, direction='both')
        varnames.rvip <- attr(mod.rvip.raic$terms,'term.labels')
      }, error=function(e){print(e);print('WARNING: error in AIC selection, keeping all variables...')}, finally={})
      
      if(length(varnames.rvip) > 0){
        mod.rvip.raic.fml <- paste0(yname, '~', paste0(varnames.rvip, collapse="+"))
        r <- buildModel.cv(mod.rvip.raic.fml, mod.data.tval, mod.type, mod.method)
        mod.rvip.raic.cv <- r[[1]]
      }
    }
  }
  
  stats.vals1 <- NA
  stats.vals2 <- NA
  
  if(exists('mod') && exists('mod.cv')){
    print.seg2("Summary of full model")
    stats.vals1 <- summarizeModel(mod, mod.cv, yname, mod.data.tval, mod.data.test, 
                                  paste0(mod.name, '_fl'), mod.type, mod.method, dirout)
    #save model results
    modelobj.save <- list(mod, mod.cv, yname, mod.data.tval, mod.data.test, paste0(mod.name, '_fl'), mod.type, mod.method)
    save(modelobj.save, 
         file=file.path(dirout, paste0(mod.name, '_fl.rda')))
  }
  
  if(exists('mod.rvip.raic') && exists('mod.rvip.raic.cv')){
    print.seg2("Summary of reduced model, by VIP recursive elimination and then by stepwise AIC")
    stats.vals2 <- summarizeModel(mod.rvip.raic, mod.rvip.raic.cv, yname, mod.data.tval, mod.data.test, 
                                  paste0(mod.name, '_rd'), mod.type, mod.method, dirout)
  }
  
  #return  
  return(list(full=stats.vals1, reduced=stats.vals2))
}

genanyzModel <- function(varnames, yname, mod.data, resdf, mod.type, mod.method, mod.rd, pfx, mod.name, dirout=dir.out){
  #generates and analyzes linear regression model, with default dataset named df.merged
  #Inputs:
  # varnames: vector of variable names
  # yname: y variable name
  # mod.data: model data
  # resdf: results table
  # mod.type: model type, 'reg' or 'bin' for regression or binary/classification
  # mod.method
  # pfx: prefix, e.g. cancertype
  # mod.name: model name, the actual model name will be appended with the cancer type
  
  #define model
  res.stats <- NA
  if(length(varnames)>0)
    res.stats <- assessModel(varnames, yname, mod.data, mod.name, mod.type, mod.method, mod.rd, dirout)
  
  #initialize
  suffx <- c('_fml' ,'_vals')
  
  colnames.fl <- paste0(pfx, '_full_', suffx)
  colnames.rd <- paste0(pfx, '_rd_', suffx)
  
  #create tmp results table
  resdf.tmp <- data.frame(t(rep(NA,5)))
  colnames(resdf.tmp) <- c('model_name', colnames.fl, colnames.rd)
  resdf.tmp$model_name <- mod.name
  
  #overwrite mod.name to tag on the cancertype
  mod.name <- paste0(pfx, mod.name)
  
  #fill in for the full model
  if(!is.na(res.stats) && length(res.stats$full) > 1 && !is.na(res.stats$full)){
    vals.list <- res.stats$full[!names(res.stats$full) %in% "fml"]
    vals <- paste0(paste0(names(vals.list),':'), collapse=";", use.names=sprintf("%.2f",vals.list))
    resdf.tmp[1, colnames.fl] <- c(res.stats$full$fml, vals)
  }
  
  #fill in for the reduced model
  if(!is.na(res.stats) && length(res.stats$reduced) > 1 && !is.na(res.stats$reduced)){
    vals.list <- res.stats$reduced[!names(res.stats$reduced) %in% "fml"]
    vals <- paste0(paste0(names(vals.list),':'), collapse=";", use.names=sprintf("%.2f",vals.list))
    resdf.tmp[1, colnames.rd] <- c(res.stats$reduced$fml, vals)
  }
  
  #append to results table
  if(nrow(resdf) < 1) resdf <- resdf.tmp
  else resdf <- rbind(resdf, resdf.tmp)
  
  return(resdf)
}

plot.km <- function(df, varname=NULL, conf.int=FALSE, mark.time=TRUE, legloc='bottomleft', timeunit='Days'){
  #plot kaplan-meier survival estimate
  if(!is.null(varname)){
    eval(parse(text=sprintf(" sf <- survfit(Surv(duration, event==1)~%s, data=df, conf.type='log', type='kaplan-meier')",varname)))
    varname.unique <- unique(df[, varname])
    varname.unique.s <- summary(df[varname])
    plot(sf, xlab=sprintf("Time (%s)", timeunit), ylab='Survival Probability',
         conf.int=conf.int, mark.time=mark.time,
         main='Kaplan-Meier curve', lwd=2, col=1:length(varname.unique))
    legend(legloc,legend=varname.unique.s, col=1:length(varname.unique), bty='n', horiz=FALSE, lty=1, lwd=2)
  } else {
    sf <- survfit(Surv(duration, event==1)~1, data=df, conf.type='log', type='kaplan-meier')
    plot(sf, xlab=sprintf("Time (%s)", timeunit), ylab='Survival Probability',
         conf.int=conf.int, mark.time=mark.time,
         main='Kaplan-Meier curve', lwd=2)
  }
}

#------- dataset helper methods --------
retrieve.dsets.clin <- function(cancertype){
  #retrieves clinical dataset
  df.clin.raw <- read.table(file.path(dir.db, paste0(cancertype, '_clin'), paste0(cancertype, '.clin.merged.picked.txt')), 
                            header=T, stringsAsFactors=FALSE, fill=T, sep='\t', quote='')
  
  #check for mutual exclusivity between the death and follow-up rows
  df.clin.t <- data.frame(t(df.clin.raw),stringsAsFactors=FALSE)
  colnames(df.clin.t) <- df.clin.raw$Hybridization.REF
  df.clin.t <- df.clin.t[-1,]
  d1 <- !is.na(df.clin.t$days_to_death)
  d2 <- !is.na(df.clin.t$days_to_last_followup)
  d.sum <- d1+d2
  if(sum(d.sum>1)){
    print('WARNING: the two columns day to death and days to last follow-up are not mutually exclusive...')
    print('WARNING: something may be wrong, there are entries where both days to death and days to last follow-up are both non-NA...')
  }
  
  if('years_to_birth' %in% colnames(df.clin.t)) df.clin.t$years_to_birth <- as.numeric(df.clin.t$years_to_birth)
  df.clin.t$days_to_death <- as.numeric(df.clin.t$days_to_death)
  df.clin.t$days_to_last_followup <- as.numeric(df.clin.t$days_to_last_followup)
  
  #construct survival data frame
  pid <- gsub('\\.','-',toupper(row.names(df.clin.t)))
  df.clin <- data.frame(patientid=pid,
                        duration=as.numeric(NA),
                        event=as.numeric(as.character(df.clin.t$vital_status)),
                        stringsAsFactors = FALSE)
  df.clin$duration[df.clin.t$vital_status==1] <- df.clin.t$days_to_death[df.clin.t$vital_status==1]
  df.clin$duration[df.clin.t$vital_status==0] <- df.clin.t$days_to_last_followup[df.clin.t$vital_status==0]
  df.clin$duration <- df.clin$duration/365.2
  df.clin <- cbind(df.clin, df.clin.t)
  
  #define must-have fields
  if(!'gender' %in% colnames(df.clin)) df.clin$gender <- NA
  if(!'age' %in% colnames(df.clin)) df.clin$age <- NA
  if(!'tnm_stage' %in% colnames(df.clin)) df.clin$tnm_stage <- NA
  
  #parse parameters
  if('years_to_birth' %in% colnames(df.clin.t)){
    df.clin$age <- as.numeric(df.clin.t$years_to_birth)/365.2
  }
  
  if('pathologic_stage' %in% colnames(df.clin.t)){
    df.clin$tnm_stage <- df.clin$pathologic_stage
    df.clin$tnm_stage <- sub('stage ','',df.clin$tnm_stage)
    df.clin$tnm_stage <- toupper(df.clin$tnm_stage)
    df.clin$tnm_stage <- sub('[A,B,C,D]*$','',df.clin$tnm_stage)
    df.clin$tnm_stage <- factor(df.clin$tnm_stage)
  }
  
  if('gender' %in% colnames(df.clin.t)) df.clin$gender <- factor(df.clin$gender)
  if('race' %in% colnames(df.clin.t)) df.clin$race <- factor(df.clin$race)
  
  #remove patients with missing data for the duration
  df.clin <- df.clin[!is.na(df.clin$duration),]
  
  return(df.clin)
}

append.clinVars <- function(df.merged.clin, df.rnaseq.PDL1){
  #append TMB and PDL1 metrics to the data.frame
  #For TMB, derived from other colunmns in df.merged.clin
  #For PDL1, need as input a data.frame with column patientid and PDL1 values
  
  #append TMB info
  df.merged.clin$TMB <- (df.merged.clin$missense+df.merged.clin$nonstop+df.merged.clin$nonsense+df.merged.clin$frameshift)/38
  df.merged.clin$TMB.cat <- 'low'
  df.merged.clin$TMB.cat[df.merged.clin$TMB > median(df.merged.clin$TMB)] <- 'high'
  df.merged.clin$TMB.cat <- factor(df.merged.clin$TMB.cat, levels=c('low','high'))
  
  #append PDL1 info
  df.merged.clin <- df.merged.clin %>% inner_join(df.rnaseq.PDL1, by='patientid')
  df.merged.clin$PDL1 <- as.numeric(df.merged.clin$PDL1)
  
  df.merged.clin$PDL1.cat0 <- 'negative'
  df.merged.clin$PDL1.cat0[df.merged.clin$PDL1 > 195.664] <- 'positive'
  df.merged.clin$PDL1.cat0 <- factor(df.merged.clin$PDL1.cat0, levels=c('negative','positive'))
  
  q1 <- quantile(df.merged.clin$PDL1)['25%']
  q3 <- quantile(df.merged.clin$PDL1)['75%']
  df.merged.clin$PDL1.cat <- 'med'
  df.merged.clin[df.merged.clin$PDL1 <= q1, 'PDL1.cat'] <- 'low'
  df.merged.clin[df.merged.clin$PDL1 >= q3, 'PDL1.cat'] <- 'high'
  df.merged.clin$PDL1.cat <- factor(df.merged.clin$PDL1.cat, levels=c('low','med','high'))
  
  return(df.merged.clin)
}

retrieveRNASeqGene <- function(cancertypes, geneName){
  #collect rna-seq gene expression for given gene across cancertypes into a single data frame
  
  df.rnaseq.PDL1 <- data.frame()
  for(cancertype in cancertypes){
    load(file.path(dir.out, paste0(cancertype, '.rda')))
    df <- as.data.frame(t(df.rnaseq[df.rnaseq$gene == geneName,]))
    df$pid <- row.names(df)
    colnames(df) <- c('PDL1', 'patientid')
    df <- df[!row.names(df)%in%c('gene.id','gene'),]
    
    if(nrow(df.rnaseq.PDL1) < 1){
      df.rnaseq.PDL1 <- df
    } else {
      df.rnaseq.PDL1 <- rbind(df.rnaseq.PDL1, df)
    }
  }
  
  return(df.rnaseq.PDL1)
}

#------- dataset methods --------
#retrieve data
retrieve.dsets <- function(cancertype, dirout=dir.out, reanalyze=FALSE){
  #Retrieves data
  #Requires: dir.out, dir.db
  #saves the parsed datasets into .rda
  #returns df.merged
  if(!file.exists(dirout)) dir.create(dirout)
  print.seg(sprintf('Retrieving data for cancer type %s', cancertype))
  
  #check if .rda already exists and already been analyzed, if so that reanalyze=F, then skip retrieval
  if(file.exists(file.path(dir.out, paste0(cancertype, '.rda'))) && !reanalyze){
    print('There is already a .rda file for this cancer type, will use this instead. Skipping retrieval method...')
    return(TRUE)
  }
  
  #----- acquire expression RNA-seq data
  #format: each row is gene, each column is patient
  print('Acquiring RNA-seq data...')
  fname <- paste0(cancertype,'_rnaseq/')
  fname.txt <- 'rsem_genes.txt'
  df.rnaseq <- read.table(file.path(dir.db, fname, fname.txt), 
                          header=T, stringsAsFactors=F, fill=T, sep='\t', quote='')
  
  #keep only scaled_estimate, and convert to TPM+0.01
  #TPM = scaled_estimate (tau) * 1e6
  col2keep.n <- c('gene_id', 'scaled_estimate')
  col2keep <- df.rnaseq[1,] %in% col2keep.n
  df.rnaseq <- df.rnaseq[-1,col2keep] #only keep scaled_estimate, then only the extra header (row 1)
  df.rnaseq[-1] <- sapply(df.rnaseq[-1], as.numeric) #convert cols to numeric
  df.rnaseq[-1] <- df.rnaseq[-1]*1e6+0.01 #TPM+0.01
  
  #change header with formatted patient id (TCGA-XX-XXXX)
  pids <- gsub('\\.', '-', colnames(df.rnaseq)) #replace . with -
  pids <- sub('(-[[:alnum:]]*){4}$', '', pids) #remove info after sample number
  colnames(df.rnaseq) <- sub('(-[[:alnum:]]*){1}$', '', pids)
  colnames(df.rnaseq)[1] <- 'gene.id'
  
  #remove duplicates, so only to keep one unique patient
  toremoveidx <- getIdxDupPatientIDs(pids)
  if(length(toremoveidx) > 0){
    print(sprintf('There are %.0f duplicates due to tumor type and/or vial, removing them...', length(toremoveidx)))
    df.rnaseq <- df.rnaseq[, !(1:ncol(df.rnaseq) %in% toremoveidx)]
  }
  
  #extract just the gene name
  df.rnaseq$gene <- sub('\\|[[:alnum:]]*$', '', df.rnaseq$gene.id)
  
  #remove ? genes, these were ones with gene id, but no hugo name
  df.rnaseq <- df.rnaseq[df.rnaseq$gene != '?',]
  dupN <- sum(duplicated(df.rnaseq$gene))
  if(dupN>0){
    print(sprintf('WARNING: there are %.0f duplicated gene names to be removed, double check...', dupN))
    df.rnaseq <- df.rnaseq[!duplicated(df.rnaseq$gene),]
  }
  
  #----- acquire copy number CNA data
  # fname <- paste0(cancertype,'_cna/')
  # fname.txt <- 'all_data_by_genes.txt'
  # df.cna <- read.table(file.path(dir.db, fname, fname.txt),
  #                      header=T, stringsAsFactors=F, fill=T, sep='\t')
  # pids <- sub('[A-Za-z](\\.[[:alnum:]]*){3}$', '', colnames(df.cna)[-(1:3)])
  # pids <- gsub('\\.','-',pids)
  # colnames(df.cna) <- c('gene', 'locusid', 'cytoband', pids)
  print('Acquiring CNA data...')
  fname <- paste0(cancertype,'_cna/')
  fname.txt <- 'all_thresholded.by_genes.txt'
  df.cna <- read.table(file.path(dir.db, fname, fname.txt),
                       header=T, stringsAsFactors=F, fill=T, sep='\t')
  pids <- sub('(\\.[[:alnum:]]*){3}$', '', colnames(df.cna)[-(1:3)]) #to format TCGA-XX-XXXX
  pids <- gsub('\\.','-',pids)
  colnames(df.cna) <- c('gene', 'locusid', 'cytoband', sub('(-[[:alnum:]]*){1}$', '', pids))
  
  #remove duplicates, so only to keep one unique patient
  toremoveidx <- getIdxDupPatientIDs(c('gene', 'locusid', 'cytoband', pids))
  if(length(toremoveidx) > 0){
    print(sprintf('There are %.0f duplicates due to tumor type and/or vial, removing them...', length(toremoveidx)))
    df.cna <- df.cna[, !(1:ncol(df.cna) %in% toremoveidx)] #ignore the first three columns, which are gene, locusid, cytoband
  }
  
  #----- acquire mutation SNV data
  #format: each file is one patient; each row per file is a gene
  print('Acquiring SNV data...')
  fname <- paste0(cancertype,'_mut/')
  flist <- list.files(path=file.path(dir.db,fname), pattern="*.maf.txt")
  
  df.snv <- data.frame() #mutation counts, each row=patient
  df.spgenes <- data.frame() #info on specific genes of interest, each row=patient
  exprWT <- data.frame()
  exprFs <- data.frame()
  exprNs <- data.frame()
  pids <- sub(".maf.txt",'',flist) #patient ids - tumortype
  toremoveidx <- getIdxDupPatientIDs(pids) #only keep one unique patient
  if(length(toremoveidx) >  0){
    print(sprintf('There are %.0f duplicates due to tumor type and/or vial, removing them...', length(toremoveidx)))
    pids <- pids[-toremoveidx]
  }
  for(f in pids){
    #SNV data for this specific patient
    df.snv.p <- read.table(file.path(dir.db, fname, paste0(f, '.maf.txt')), header=T, stringsAsFactors=F, fill=T, sep='\t', quote='')
    
    #derive parameter values from SNV
    patientid <- sub('(-[[:alnum:]]*){1}$', '', f) #patient id (remove the tumortype tag)
    #df.snv.p$vaf <- df.snv.p$t_alt_count / (df.snv.p$t_alt_count + df.snv.p$t_ref_count) #allele variant frequency
    
    #rna-seq data for this specific patient
    df.rnaseq.p <- data.frame(genes=df.rnaseq$gene, expr=NA)
    if(patientid %in% colnames(df.rnaseq))
      df.rnaseq.p$expr <- df.rnaseq[,colnames(df.rnaseq)==patientid]
    
    #CNA data for this specific patient
    df.cna.p.all <- data.frame(genes=df.cna$gene, cna=0)
    if(patientid %in% colnames(df.cna))
      df.cna.p.all$cna <- df.cna[,colnames(df.cna)==patientid]
    df.cna.p <- df.cna.p.all[df.cna.p.all$cna!=0,] #only keep the genes that have CNA
    
    #get mutations
    df.ns <- df.snv.p[df.snv.p$Variant_Classification=='Nonsense_Mutation',]
    df.fs <- df.snv.p[df.snv.p$Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins'),]
    df.mut <- df.snv.p[df.snv.p$Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins',
                                                              'Nonstop_Mutation', 'In_Frame_Del', 'Nonsense_Mutation',
                                                              'Translation_Start_Site', 'Splice_Site'),]
    
    #for PTC-bearing mutations only keep ones with VAF>=20%
    #Not used at the moment
    # df.ns <- df.ns[df.ns$vaf>=0.2,]
    # df.fs <- df.fs[df.fs$vaf>=0.2,]
    # df.mut <- df.mut[df.mut$vaf>=0.2,]
    
    #get WT transcripts (expression)
    exprWT.p <- df.rnaseq.p
    if(nrow(df.mut) > 0)
      exprWT.p[which(exprWT.p$genes %in% df.mut$Hugo_Symbol), 'expr'] <- NA #remove transcripts with mutation
    if(nrow(df.cna.p) > 0)
      exprWT.p[which(exprWT.p$genes %in% df.cna.p$genes), 'expr'] <- NA #remove transcripts with CNA
    exprWT.p <- data.frame(p=exprWT.p[,2])
    colnames(exprWT.p) <- patientid
    
    # exprMut.p <- df.rnaseq.p
    # exprMut.p[which(!(exprMut.p$genes %in% df.mut$Hugo_Symbol)), 'expr'] <- NA
    # exprMut.p <- data.frame(p=exprMut.p[,2])
    # colnames(exprMut.p) <- patientid
    
    #get PTC transcripts (expression) with only fs mutation
    exprFs.p <- df.rnaseq.p
    if(nrow(df.mut) > 0)
      exprFs.p[which(!(exprFs.p$genes %in% df.fs$Hugo_Symbol)), 'expr'] <- NA
    if(nrow(df.cna.p) > 0)
      exprFs.p[which(exprFs.p$genes %in% df.cna.p$genes), 'expr'] <- NA #remove transcripts with CNA
    exprFs.p <- data.frame(p=exprFs.p[,2])
    colnames(exprFs.p) <- patientid
    
    #get PTC transcripts (expression) with only ns mutation
    exprNs.p <- df.rnaseq.p
    if(nrow(df.mut) > 0)
      exprNs.p[which(!(exprNs.p$genes %in% df.ns$Hugo_Symbol)), 'expr'] <- NA
    if(nrow(df.cna.p) > 0)
      exprNs.p[which(exprNs.p$genes %in% df.cna.p$genes), 'expr'] <- NA #remove transcripts with CNA
    exprNs.p <- data.frame(p=exprNs.p[,2])
    colnames(exprNs.p) <- patientid
    
    if(nrow(exprWT) < 1) exprWT <- data.frame(genes=df.rnaseq.p$genes)
    exprWT <- cbind(exprWT, exprWT.p)
    if(nrow(exprFs) < 1) exprFs <- data.frame(genes=df.rnaseq.p$genes)
    exprFs <- cbind(exprFs, exprFs.p)
    if(nrow(exprNs) < 1) exprNs <- data.frame(genes=df.rnaseq.p$genes)
    exprNs <- cbind(exprNs, exprNs.p)
    
    #SNV tally
    df.snv.s <- data.frame(patientid=patientid,
                           silent=sum(df.snv.p$Variant_Classification=='Silent'),
                           missense=sum(df.snv.p$Variant_Classification=='Missense_Mutation'),
                           nonstop=sum(df.snv.p$Variant_Classification=='Nonstop_Mutation'),
                           nonsense=sum(df.snv.p$Variant_Classification=='Nonsense_Mutation'),
                           frameshift=length(grep("Frame_Shift", df.snv.p$Variant_Classification)),
                           SNP=sum(df.snv.p$Variant_Type=='SNP'),
                           DEL=sum(df.snv.p$Variant_Type=='DEL'),
                           INS=sum(df.snv.p$Variant_Type=='INS'),
                           total=nrow(df.snv.p) )
    df.snv <- rbind(df.snv, df.snv.s)
    
    #retrieve info for the genes of interest (expr, cna, mut)
    #spgenes.p <- df.rnaseq.p[which(df.rnaseq.p$genes %in% genesList), c('genes','expr')]
    spgenes.p <- df.cna.p.all[which(df.cna.p.all$genes %in% genesList), c('genes','cna')]
    genes.mut <- c(df.ns$Hugo_Symbol, df.fs$Hugo_Symbol, df.mut$Hugo_Symbol)
    spgenes.p$mut <- sapply(spgenes.p$genes, function(x)(sum(x == genes.mut)))
    vals <- as.vector(as.matrix(t(spgenes.p[,-1])))
    vals.names <-  as.vector(sapply(spgenes.p$genes, function(x)(paste0(x, c('_cna','_mut')))))
    spgenes.p$genes <- droplevels(spgenes.p[['genes']])
    spgenes.ps.t <- data.frame(t(vals))
    spgenes.ps.t <- cbind(patientid=patientid, spgenes.ps.t) 
    colnames(spgenes.ps.t) <- c('patientid', vals.names)
    df.spgenes <- rbind.fill(df.spgenes, spgenes.ps.t)
  }
  
  #get the median WT expr
  exprWT$genes <- as.character(exprWT$genes)
  exprWT.med <- data.frame(genes=exprWT$genes,
                           median=apply(exprWT[,-1], 1, median, na.rm=T),
                           mean=apply(exprWT[,-1], 1, mean, na.rm=T),
                           sd=apply(exprWT[,-1], 1, sd, na.rm=T),
                           n=apply(exprWT[,-1], 1, function(x) sum(!is.na(x))))
  exprWT.med$cv <- exprWT.med$sd/exprWT.med$mean
  exprWT.med$pctile <- (ecdf(exprWT.med$median))(exprWT.med$median) #percentile
  
  exprNs$genes <- as.character(exprNs$genes) #nonsense
  exprFs$genes <- as.character(exprFs$genes) #frameshift
  exprPTC <- exprNs #for PTC, a combined nonsense + frameshift
  exprPTC[!is.na(exprFs)] <- exprFs[!is.na(exprFs)]
  exprPTC[,-1] <- as.data.frame(sapply(exprPTC[,-1], as.numeric))
  
  #filtering
  genes.exclude <- c()
  
  #exclude genes with WT median expression < 5 TPM
  genes.exclude <- c(genes.exclude, exprWT.med[exprWT.med$median < 5, 'genes'])
  
  #exclude genes with large variation in WT expression, CV > 0.5
  genes.exclude <- c(genes.exclude, exprWT.med[exprWT.med$cv > 0.5, 'genes'])
  
  #exclude genes with WT samples/measurements < 10
  genes.exclude <- c(genes.exclude, exprWT.med[exprWT.med$n < 10, 'genes'])
  
  #for the genes to exclude, place NA on the median
  exprWT.med$median[exprWT.med$genes %in% genes.exclude] <- NA
  
  #NMD efficiency, for each individual gene
  #NMD eff = 0:no decay (not efficient); 1:half decayed. Higher val, moree NMD efficient
  exprNs.NMD <- -log2( exprNs[,-1] / exprWT.med$median ) #NMD efficiency -log2(expr nonsense / median expr WT)
  exprFs.NMD <- -log2( exprFs[,-1] / exprWT.med$median ) #NMD efficiency -log2(expr frameshift / median expr WT)
  exprPTC.NMD <- -log2( exprPTC[,-1] / exprWT.med$median ) #NMD efficiency -log2(expr PTC / median expr WT)
  
  #NMD burden, for patient
  getNMDfeat <- function(x){
    wts <- exprWT.med$pctile #get weights, normalize so the wts sum to 1
    wts[is.na(x)] <- NA
    wts[!is.na(wts)] <- wts[!is.na(wts)] / sum(wts, na.rm=T)
    x.wt <- x * wts #weighted by percentile of the gene expression (the weight is normalized)
    r <- data.frame(max=max(x, na.rm=T),
                    min=min(x, na.rm=T),
                    n_decayed=sum(x[!is.na(x)]>0), #number of ones with some decay
                    frac_decayed=sum(x[!is.na(x)]>0)/sum(!is.na(x)), #fraction of mut genes with decay
                    
                    med=median(x, na.rm=T),
                    med.wt=median(x.wt, na.rm=T),
                    mean=mean(x, na.rm=T),
                    mean.wt=mean(x.wt, na.rm=T) )
    
    if(sum(!is.na(x)) < 1) r[,] <- NA #if there are no alterations, NA all
    
    r
  }
  
  nmd.fs <- do.call(rbind, lapply(exprFs.NMD, getNMDfeat)) #frameshift
  nmd.ns <- do.call(rbind, lapply(exprNs.NMD, getNMDfeat)) #nonsense
  nmd.ptc <- do.call(rbind, lapply(exprPTC.NMD, getNMDfeat)) #fs + ns
  
  #append patient id
  nmd.fs$patientid <- row.names(nmd.fs)
  nmd.ns$patientid <- row.names(nmd.ns)
  nmd.ptc$patientid <- row.names(nmd.ptc)
  row.names(nmd.fs) <- NULL
  row.names(nmd.ns) <- NULL
  row.names(nmd.ptc) <- NULL
  
  #for NAs, there are no mutations (based on threshold) for measuring NMD burden, create boolean to mark the NA
  nmd.fs.NA <- data.frame(mut.bool=ifelse(rowSums(is.na(nmd.fs[,1:8])) == 8, 0, 1))
  nmd.ns.NA <- data.frame(mut.bool=ifelse(rowSums(is.na(nmd.ns[,1:8])) == 8, 0, 1))
  nmd.ptc.NA <- data.frame(mut.bool=ifelse(rowSums(is.na(nmd.ptc[,1:8])) == 8, 0, 1))
  
  #for NAs, take on the maximum value, as there is no NMD burden
  nmd.fs <- cbind(as.data.frame(apply(nmd.fs[1:8], 2, fillNAmax)), nmd.fs.NA, data.frame(patientid=nmd.fs$patientid))
  nmd.ns <- cbind(as.data.frame(apply(nmd.ns[1:8], 2, fillNAmax)), nmd.ns.NA, data.frame(patientid=nmd.ns$patientid))
  nmd.ptc <- cbind(as.data.frame(apply(nmd.ptc[1:8], 2, fillNAmax)), nmd.ptc.NA, data.frame(patientid=nmd.ptc$patientid))
  
  #----- cytolytic activity
  #method 1 (Rooney et al): uses geometric mean of GZMA and PRF1
  cyt.genes <- c('GZMA', 'PRF1')
  cyt <- df.rnaseq[df.rnaseq$gene %in% cyt.genes,]
  rownames(cyt) <- cyt$gene
  cyt <- cyt[,grep('TCGA',colnames(cyt))] #only keep TCGA entries, drop gene columns
  cyt <- as.data.frame(t(cyt))
  cyt$cyt <- apply(cyt, 1, calcGeoMean)
  cyt$patientid <- rownames(cyt)
  rownames(cyt) <- NULL
  
  #method 2 (Spranger et al)
  #unclear of the methods
  # tcell.siggenes <- c('CD8A','CCL2','CCL3','CCL4','CXCL9','CXCL10','ICOS','GZMK','IRF1','HLA-DMA',
  #                     'HLA-DMB','HLA-DOA','HLA-DOB')
  # cyt2 <- df.rnaseq[df.rnaseq$gene %in% tcell.siggenes,]
  # rownames(cyt2) <- cyt2$gene
  # cyt2 <- cyt2[,grep('TCGA',colnames(cyt2))] #only keep TCGA entries, drop gene columns
  # cyt2 <- as.data.frame(t(cyt2))
  # cyt2$cyt <- apply(cyt2, 1, calcGeoMean)
  
  #relevant extracted variables - per tumor
  #df.rnaseq, df.cna
  #df.snv, nmd.fs, nmd.ns, nmd.ptc, cyt
  
  #----- microsatellite instability 
  print('Acquiring MSI information...')
  load(file.path(dir.db, 'MSI', 'sampledata_post_review_070716.robj'))
  df.msi <- data.frame(patientid=sampledat$sample_name, msi=sampledat$msi)
  levels(df.msi) <- c('MSS','MSI-H') #change order of levels so reference is MSS
  
  #----- clinical data
  df.clin <- retrieve.dsets.clin(cancertype) #generates duration and event columns
  
  #----- combine the datasets
  exc.idx <- which(colnames(nmd.fs) == 'patientid')
  nmd.ptc.tmp <- nmd.ptc
  nmd.fs.tmp <- nmd.fs
  nmd.ns.tmp <- nmd.ns
  colnames(nmd.ptc.tmp)[-exc.idx] <- paste0('nmdptc_', colnames(nmd.ptc.tmp)[-exc.idx])
  colnames(nmd.fs.tmp)[-exc.idx] <- paste0('nmdfs_', colnames(nmd.fs.tmp)[-exc.idx])
  colnames(nmd.ns.tmp)[-exc.idx] <- paste0('nmdns_', colnames(nmd.ns.tmp)[-exc.idx])
  
  df.merged <- merge(df.snv, cyt, by='patientid')
  df.merged <- merge(df.merged, nmd.ptc.tmp, by='patientid')
  df.merged <- merge(df.merged, nmd.fs.tmp, by='patientid')
  df.merged <- merge(df.merged, nmd.ns.tmp, by='patientid')
  df.merged <- merge(df.merged, df.spgenes, by='patientid')
  df.merged <- merge(df.merged, df.msi, by='patientid', all.x=TRUE) 
  row.names(df.merged) <- df.merged$patientid
  df.merged <- df.merged[-1]
  
  #categorize CYT
  q1 <- quantile(df.merged$cyt)['25%']
  q3 <- quantile(df.merged$cyt)['75%']
  
  #cat: low, med, high
  df.merged$cyt.cat <- 'med'
  df.merged$cyt.cat[df.merged$cyt <= q1] <- 'low'
  df.merged$cyt.cat[df.merged$cyt >= q3] <- 'high'
  df.merged$cyt.cat <- factor(df.merged$cyt.cat, levels=c('low','med','high'))
  
  #bool: boolean of just high or low, med is NA
  df.merged$cyt.bool <- NA
  df.merged$cyt.bool[df.merged$cyt <= q1] <- 0
  df.merged$cyt.bool[df.merged$cyt >= q3] <- 1
  df.merged$cyt.bool <- factor(df.merged$cyt.bool, levels=c(0,1))
  
  #filters
  #has these variables initially, but never really modeled with them, so considered to remove them here...
  df.merged <- df.merged[, !colnames(df.merged) %in% c('SNP','DEL','INS', 'nmdptc_min','nmdns_min','nmdfs_min') ]
  
  #deal with missing data
  levels(df.merged$msi) <- c('MSI-H', 'MSS', 'X')
  df.merged$msi[is.na(df.merged$msi)] <- 'X' #change missing MSI data info marked X's
  df.merged <- handleMissingVals(df.merged) #checking missing data, and if small, drop
  
  #drop invariant variables
  df.merged.bfdrop <- df.merged #save a copy of the data frame before droping invariant variables
  df.merged <- rmConstVars(df.merged)
  
  #-- additional derived dataset from df.merged
  #clinical features merge
  df.merged.clin.bfdrop <- df.merged.bfdrop %>% mutate(patientid=rownames(.)) %>% left_join(df.clin, by="patientid")
  df.merged.clin <- df.merged %>% mutate(patientid=rownames(.)) %>% left_join(df.clin, by="patientid")
  df.merged.clin <- df.merged.clin %>% handleMissingVals %>% rmConstVars %>% cleanup.factors
  
  #separate data.frame with just non-NA CYT.bool values, i.e. just low/high
  df.merged.bool <- getdfbool(df.merged)
  df.merged.bool <- df.merged.bool %>% handleMissingVals %>% rmConstVars %>% cleanup.factors
  
  #save merged dataset (input matrix)
  print('Saving parsed dataset...')
  write.csv(df.merged, file=file.path(dirout, paste0('input_',cancertype,'.csv')))
  
  #save objects
  save(df.rnaseq, df.cna, 
       df.snv, nmd.fs, nmd.ns, nmd.ptc, cyt,
       df.spgenes,
       df.merged.bfdrop, df.merged.clin.bfdrop,
       df.merged, df.merged.bool, df.merged.clin,
       file=file.path(dir.out, paste0(cancertype, '.rda')))
}

#create pan-cancer
createPanCancer <- function(cancertypes, dirout=dir.out, reanalyze=FALSE, appendIndication=FALSE, fname='pancancer'){
  #appendIndication: indicates whether to add a new column called 'cancertype'
  if(!file.exists(dirout)) dir.create(dirout)
  print('Creating pan-cancer dataset...')
  
  #check if .rda already exists and already been analyzed, if so that reanalyze=F, then skip retrieval
  if(file.exists(file.path(dir.out, 'pancancer.rda')) && !reanalyze){
    print('There is already a .rda file for pan-cancer, will use this instead. Skipping retrieval method...')
    return(TRUE)
  }
  
  df.merged.tot <- data.frame()
  df.merged.bool.tot <- data.frame()
  df.merged.clin.tot <- data.frame()
  for(cancertype in cancertypes){
    load(file.path(dir.out, paste0(cancertype, '.rda')))
    
    if(appendIndication)
      df.merged.bfdrop$cancertype <- cancertype
    df.merged.tot <- rbind(df.merged.tot, df.merged.bfdrop)
    df.merged.bool.tot <- rbind(df.merged.bool.tot, getdfbool(df.merged.bfdrop))
    
    if(appendIndication)
      df.merged.clin.bfdrop$cancertype <- cancertype
    df.merged.clin.rd <- df.merged.clin.bfdrop %>% select(colnames(df.merged.tot), 
                                                          patientid, duration, event, age, gender, tnm_stage,
                                                          cancertype) #reduced list of vars
    df.merged.clin.tot <- rbind(df.merged.clin.tot, df.merged.clin.rd)
  }
  
  df.merged <- df.merged.tot
  df.merged.bool <- df.merged.bool.tot
  df.merged.clin <- df.merged.clin.tot
  
  #deal with missing data, invariant variables, and remove empty levels in factors
  df.merged <- df.merged %>% handleMissingVals %>% rmConstVars %>% cleanup.factors
  df.merged.bool <- df.merged.bool %>% handleMissingVals %>% rmConstVars %>% cleanup.factors
  df.merged.clin <- df.merged.clin %>% handleMissingVals %>% rmConstVars %>% cleanup.factors
  
  #save merged dataset (input matrix)
  write.csv(df.merged, file=file.path(dirout, paste0('input_',fname,'.csv')))
  
  save(df.merged, df.merged.bool, df.merged.clin,
       file=file.path(dir.out, paste0(fname, '.rda')))
}

#preprocessing
preprocess.dsets <- function(cancertype, dirout=dir.out, reanalyze=FALSE){
  #preprocess dataset
  #if dbset.preprocess.std, standardize
  #if dbset.preprocess.pca, perform PCA on data
  if(!file.exists(dirout)) dir.create(dirout)
  
  #check if .rda already exists and already been analyzed, if so that reanalyze=F, then skip proprocessing
  if(file.exists(file.path(dir.out, paste0(cancertype, '_prep_base.rda'))) && !reanalyze){
    print('There is already a _prep_base.rda file for this cancer type, will use this instead. Assume other transformation
          if applicable has also been completed. If need to rerun, set reanalyze to TRUE. Skipping preprocessing... ')
    return(TRUE)
  }
  
  #load data, and perform any standard preprocessing
  load(file.path(dir.out, paste0(cancertype, '.rda')))
  
  #clean up
  df.merged <- cleanup.factors(df.merged)
  
  #separate data.frame with just non-NA CYT.bool values, i.e. just low/high
  df.merged.bool <- getdfbool.rm(df.merged)
  
  #filter MSI, if there aren't that many MSI-H and MSS, then exclude MSI
  checkparseMSI <- function(df.merged){
    if('msi' %in% colnames(df.merged)){
      notincMSI <- FALSE
      
      if(sum(df.merged$msi=='MSI-H') < 5){
        print('There are less than 5 cases of MSI-H... will not include MSI as a variable in analysis')
        notincMSI <- TRUE
      }
      
      if(sum(df.merged$msi=='MSS') < 5){
        print('There are less than 5 cases of MSS... will not include MSI as a variable in analysis')
        notincMSI <- TRUE
      }
      
      if(notincMSI) df.merged <- dplyr::select(df.merged, -msi)
      
      if(dbset.preprocess.MSIcomplete){
        df.merged <- df.merged[df.merged$msi %in% c('MSI-H', 'MSS'),]
      }
      
      df.merged <- cleanup.factors(df.merged)
    }
    
    return(df.merged)
  }
  
  df.merged <- checkparseMSI(df.merged)
  df.merged.bool <- checkparseMSI(df.merged.bool)
  df.merged.clin <- checkparseMSI(df.merged.clin)
  
  #save
  save(df.merged, df.merged.bool, df.merged.clin, file=file.path(dir.out, paste0(cancertype, '_prep_base.rda')))
  
  if(dbset.preprocess.std){
    #standardized dataset
    #load data
    load(file.path(dir.out, paste0(cancertype, '_prep_base.rda')))
    
    #define variable names
    varsY <- c('cyt','GZMA','PRF1','cyt.bool','cyt.cat')
    varsX <- colnames(df.merged)
    varsX <- varsX[!(varsX %in% varsY)]
    
    #keep categorical variables out - do not standardize them
    vars2std <- varsX[!varsX %in% c('nmdns_mut.bool', 'nmdfs_mut.bool', 'nmdptc_mut.bool')]
    preproc.params <- caret::preProcess(df.merged[,colnames(df.merged) %in% vars2std], method=c("center", "scale"))
    df.merged[,colnames(df.merged) %in% vars2std] <- predict(preproc.params, df.merged[,colnames(df.merged) %in% vars2std])
    
    # vars2std <- varsX[!varsX %in% c('nmdns_mut.bool', 'nmdfs_mut.bool', 'nmdptc_mut.bool')]
    # preproc.params <- caret::preProcess(df.merged[,colnames(df.merged.bool) %in% vars2std], method=c("center", "scale"))
    # df.merged.bool[,colnames(df.merged.bool) %in% vars2std] <- predict(preproc.params, df.merged.bool[,colnames(df.merged.bool) %in% vars2std])
    
    #clean up
    df.merged <- cleanup.factors(df.merged)
    
    #separate data.frame with just non-NA CYT.bool values, i.e. just low/high
    df.merged.bool <- getdfbool.rm(df.merged)
    
    #save merged dataset (input matrix)
    write.csv(df.merged, file=file.path(dirout, paste0('input_',cancertype,'_prep_std.csv')))
    write.csv(df.merged.bool, file=file.path(dirout, paste0('input_',cancertype,'_prep_std_bool.csv')))
    save(df.merged, df.merged.bool, file=file.path(dir.out, paste0(cancertype, '_prep_std.rda')))
  }
  
  if(dbset.preprocess.pca){
    #load data
    load(file.path(dir.out, paste0(cancertype, '_prep_base.rda')))
    
    if(nrow(df.merged) < 1){
      stop("ERROR: the input dataset is empty, cannot create PCA based on an empty dataset")
    }
    
    #PCA of the data
    varsY <- c('cyt','GZMA','PRF1','cyt.bool','cyt.cat')
    varsX <- colnames(df.merged)
    varsX <- varsX[!(varsX %in% varsY)]
    
    #exclude categorical variables
    varsX <- varsX[!sapply(df.merged[varsX],is.factor)]
    
    #get the X values
    dataX <- df.merged[varsX]
    dataX <- scale(dataX)
    
    df.classes <- rep(0,ncol(dataX))
    df.classes[colnames(dataX) %in%  c('silent','missense','nonstop','nonsense','frameshift','total','SNP','DEL','INS')] <- 1
    df.classes[colnames(dataX) %in%  colnames(dataX)[grep('^nmd' ,colnames(dataX))] ] <- 2
    df.classes[colnames(dataX) %in% c(colnames(dataX)[grep('_mut$' ,colnames(dataX))],
                                      colnames(dataX)[grep('_cna$' ,colnames(dataX))])  ] <- 3
    
    #PCA
    Xpca <- prcomp(dataX)
    PCAscores <- Xpca$x #get scores
    PCAloadings <- Xpca$rotation #get loadings
    df.merged <- cbind(PCAscores, df.merged[varsY]) #create the df.merged based on PCA scores
    
    #clean up
    df.merged <- cleanup.factors(df.merged)
    
    #separate data.frame with just non-NA CYT.bool values, i.e. just low/high
    df.merged.bool <- getdfbool.rm(df.merged)
    
    #save merged dataset (input matrix)
    write.csv(df.merged, file=file.path(dirout, paste0('input_',cancertype,'_prep_pca.csv')))
    write.csv(df.merged.bool, file=file.path(dirout, paste0('input_',cancertype,'_prep_pca_bool.csv')))
    save(df.merged, df.merged.bool, file=file.path(dir.out, paste0(cancertype, '_prep_pca.rda')))
    
    
    #PCA plots
    pdf(file.path(dirout, paste0('PCA_1vs2_lab.pdf')), height=7, width=7)
    plot(PCAscores[,1:2], pch=21, bg=df.merged$cyt.bool,cex=1.5, main="Scores")
    legend('topright',legend=levels(df.merged$cyt.cat), pch=21,col=1:2,pt.cex=1.5)
    #text(PCAscores[, 1:2],  labels=1:nrow(PCAscores))
    dev.off()
    
    pdf(file.path(dirout, paste0('PCA_1vs3_lab.pdf')), height=7, width=7)
    plot(PCAscores[,c(1,3)], pch=21, bg=df.merged$cyt.bool,cex=1.5, main="Scores")
    legend('topright',legend=levels(df.merged$cyt.cat), pch=21,col=1:2,pt.cex=1.5)
    #text(PCAscores[, 1:2],  labels=1:nrow(PCAscores))
    dev.off()
    
    pdf(file.path(dirout, paste0('PCA_1vs4_lab.pdf')), height=7, width=7)
    plot(PCAscores[,c(1,4)], pch=21, bg=df.merged$cyt.bool,cex=1.5, main="Scores")
    legend('topright',legend=levels(df.merged$cyt.cat), pch=21,col=1:2,pt.cex=1.5)
    #text(PCAscores[, 1:2],  labels=1:nrow(PCAscores))
    dev.off()
    
    pdf(file.path(dirout, paste0('PCA_1vs2.pdf')), height=7, width=7)
    plot(PCAloadings[,1:2],  pch=21, bg=df.classes, cex=1.5, main="Loadings" )
    legend('topright',legend=c('mut burden','NMD','mut/cna'), pch=21,pt.cex=1.5, col=1:3)
    #text(PCAloadings[,1:2],  labels=rownames(PCAloadings))
    dev.off()
    
    pdf(file.path(dirout, paste0('PCA_1vs3.pdf')), height=7, width=7)
    plot(PCAloadings[,c(1,3)],  pch=21, bg=df.classes, cex=1.5, main="Loadings" )
    legend('topright',legend=c('mut burden','NMD','mut/cna'), pch=21,pt.cex=1.5, col=1:3)
    #text(PCAloadings[,1:2],  labels=rownames(PCAloadings))
    dev.off()
    
    pdf(file.path(dirout, paste0('PCA_1vs4.pdf')), height=7, width=7)
    plot(PCAloadings[,c(1,4)],  pch=21, bg=df.classes, cex=1.5, main="Loadings" )
    legend('topright',legend=c('mut burden','NMD','mut/cna'), pch=21,pt.cex=1.5, col=1:3)
    #text(PCAloadings[,1:2],  labels=rownames(PCAloadings))
    dev.off()
    
    #alternatively, PCA reconstruction
    # nComps = 1:2
    # Xhat = Xpca$x[,nComps] %*% t(Xpca$rotation[,nComps]) #x is scores (data times loadings), rotation is the loadings
    # mu = colMeans(dataX)
    # Xhat = scale(Xhat, center = -mu, scale = FALSE) #add back the mean
  }
}

#analyze data
analyze.dsets.uni <- function(cancertype, res.corr=NULL, dirout=dir.out, reanalyze=F){
  if(!file.exists(dirout)) dir.create(dirout)
  print.seg(sprintf('Analyzing data for cancer type %s | univariate', cancertype))
  
  #check if this already has been analyzed, if so that reanalyze=F, then skip
  if(file.exists(file.path(dirout, 'summary_correlations.csv')) && !reanalyze){
    print('There is already a summary_correlations.csv file for this cancer type. Univariate analysis should already be done. Skipping... ')
    return(TRUE)
  }
  
  #load data
  if(dbset.preprocess.std) load(file.path(dir.out, paste0(cancertype, '_prep_std.rda')))
  else if(dbset.preprocess.pca) load(file.path(dir.out, paste0(cancertype, '_prep_pca.rda')))
  else load(file.path(dir.out, paste0(cancertype, '_prep_base.rda')))
  
  if(nrow(df.merged) < 1){
    print("WARNING: input data.frame is empty...")
    return(FALSE)
  }
  
  #define variable names
  varsY <- c('cyt','GZMA','PRF1','cyt.bool','cyt.cat')
  varsX <- colnames(df.merged)
  varsX <- varsX[!(varsX %in% varsY)]
  varsX.cat <- c(colnames(df.merged)[grep('_cna$', varsX)], 
                 colnames(df.merged)[grep('_mut$', varsX)],
                 colnames(df.merged)[grep('\\.bool$', varsX)],
                 'msi')
  
  #----- dataset general statistics
  #check for missing data, only present warning if there is (dropsmall is FALSE here)
  handleMissingVals(df.merged, dropsmall=FALSE)
  
  #descriptive statistics
  print(describe(df.merged))
  
  #----- univariate
  print('Performing univariate analysis...')
  ##cytolytic activity plots
  ggplot(df.merged, aes(GZMA,PRF1)) + geom_point(size=3, alpha=0.7)
  ggsave(file.path(dirout, paste0('gzma_prf1_',cancertype,'.pdf')))
  #heatmap.2(as.matrix(cyt2), col=topo.colors(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", labRow=F)
  
  ## correlation/collinearity
  #pairs matrix - very messy, not very useful
  # p<-ggpairs(df.merged)
  # ggsave(file.path(dirout, paste0('pairs_',cancertype,'.pdf')), plot=p)
  
  #spearman correlation
  c <- rcorr(sapply(df.merged, as.numeric), type="spearman")
  
  r <- tryCatch({
    pdf(file.path(dirout, paste0('cormatrix_',cancertype,'.pdf')), width=12, height=12)
    heatmap.2(as.matrix(c$r), 
              col=colorRampPalette(brewer.pal(9, "RdBu"))(100), scale="none", 
              key.xlab='Spearman corr', key.title=F,
              key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5, cexCol=0.5,
              margins = c(8,8))
    dev.off()
  }, error=catch.errorfunc, finally={})
  
  ##continuous CYT, scatter plot with loess
  for(xvar in varsX){
    r <- tryCatch({
      #for variables (e.g. categorical) that does not need x transformation
      if(xvar %in% varsX.cat){
        # ggplot(df.merged, aes(eval(parse(text=xvar)), cyt)) + geom_point(size=3, alpha=0.7) +
        #   xlab(xvar) + ylab('Cytolytic activity') +
        #   geom_smooth()
        # ggsave(file.path(dirout, paste0('cyt_', xvar,'_',cancertype,'.pdf')))
        
        #box-whisker plots
        #specific for MSI
        if(xvar == 'msi'){
          p1 <- ggplot(df.merged, aes(as.factor(eval(parse(text=xvar))), cyt)) + geom_boxplot() + geom_jitter(width=0.2) +
            xlab(xvar) + ylab('Cytolytic activity') + stat_compare_means(method='wilcox.test',comparisons=list(c('MSI-H','MSS')))
          p2 <- ggplot(df.merged, aes(as.factor(eval(parse(text=xvar))), cyt)) + geom_boxplot() + geom_jitter(width=0.2) +
            xlab(xvar) + ylab(bquote(log[10]~'['~.('Cytolytic activity')~']')) + scale_y_log10() + 
            stat_compare_means(method='wilcox.test',comparisons=list(c('MSI-H','MSS')))
          ggarrange(p1, p2, ncol = 2, nrow = 1)
          ggsave(file.path(dirout, paste0('cyt_', xvar,'_',cancertype,'.pdf')))
        } else {
          p1 <- ggplot(df.merged, aes(as.factor(eval(parse(text=xvar))), cyt)) + geom_boxplot() + geom_jitter(width=0.2) +
            xlab(xvar) + ylab('Cytolytic activity')
          p2 <- ggplot(df.merged, aes(as.factor(eval(parse(text=xvar))), cyt)) + geom_boxplot() + geom_jitter(width=0.2) +
            xlab(xvar) + ylab(bquote(log[10]~'['~.('Cytolytic activity')~']')) + scale_y_log10()
          ggarrange(p1, p2, ncol = 2, nrow = 1)
          ggsave(file.path(dirout, paste0('cyt_', xvar,'_',cancertype,'.pdf')))
        }
        
        #for variables (e.g. continuous) that include x transformation
      } else {
        ggplot(df.merged, aes(eval(parse(text=xvar)), cyt)) + geom_point(size=3, alpha=0.7) +
          xlab(xvar) + ylab('Cytolytic activity') + scale_x_continuous(trans='log2') +
          geom_smooth()
        ggsave(file.path(dirout, paste0('cyt_', xvar,'_',cancertype,'.pdf')))
      }
    }, error=catch.errorfunc, finally={})
  }
  
  ##categorized CYT, bar plots and statistics
  res.corr.s <- data.frame(colnames(df.merged), as.data.frame(c$r)$cyt, as.data.frame(c$P)$cyt)
  colnames(res.corr.s) <- c('var', paste0(cancertype,'_rho'), paste0(cancertype,'_rho_pval'))
  med.cat <- aggregate(.~cyt.cat,df.merged.bool,median) #statistics (median) based on low/high CYT
  med.delta <- med.cat[2,-1] - med.cat[1,-1] #difference of median betwee low and high group)
  res.corr.s[paste0(cancertype,'_LH_meddelta')] <- NA
  res.corr.s[paste0(cancertype,'_LH_pval')] <- NA
  
  for(xvar in varsX){
    r <- tryCatch({
      if(pval.ptc.excludeNoMut && grepl("^nmd.*", xvar) && !grepl(".*_mut\\.bool$", xvar)){
        #if the xvar is a NMD metric, and pval.ptc.excludeNoMut is TRUE, then will exclude patients with no corresponding
        #mutation in the plotting and calculation of p-val (based on Mann-Whitney test)
        xvar.nmd <- sub("_.*$",'',xvar) #get which PTC to look at (fs, ns, or ptc)
        df.tmp <- df.merged[df.merged[paste0(xvar.nmd,'_mut.bool')]==1,] #df with just patients with mutations for the PTC type
        df.tmp.bool <- df.merged.bool[df.merged.bool[paste0(xvar.nmd,'_mut.bool')]==1,] #df with just patients with mutations for the PTC type
      } else {
        df.tmp <- df.merged
        df.tmp.bool <- df.merged.bool
      }
      
      #box plots
      #for categorical variables
      if(xvar %in% varsX.cat){
        #chi-squared test, contingency table as 
        st <- chisq.test(table(df.tmp.bool$cyt.cat, df.tmp.bool[[xvar]]))
        
        #plots
        ggplot(df.tmp, aes(cyt.cat, eval(parse(text=xvar))) ) + geom_count() +
          xlab('Cytolytic activity') + ylab(xvar) 
        ggsave(file.path(dirout, paste0('cyt_bool_', xvar,'_',cancertype,'_LH.pdf')))
        
        #for continuous variables
      } else {
        #Mann-Whitney test, grouped by low and high CYT
        st <- wilcox.test(eval(parse(text=xvar))~cyt.bool, data=df.tmp.bool)
        
        #plots
        p1 <- ggplot(df.tmp, aes(cyt.cat, eval(parse(text=xvar))) ) + geom_boxplot() + geom_jitter(width=0.2) +
          xlab('Cytolytic activity') + ylab(xvar) + 
          stat_compare_means(method='wilcox.test',comparisons=list(c('high','low')))
        p2 <- ggplot(df.tmp, aes(cyt.cat, eval(parse(text=xvar))) ) + geom_boxplot() + geom_jitter(width=0.2) +
          xlab('Cytolytic activity') + ylab(bquote(log[10]~'['~.(xvar)~']')) + scale_y_log10() +
          stat_compare_means(method='wilcox.test',comparisons=list(c('high','low')))
        ggarrange(p1, p2, ncol = 2, nrow = 1)
        ggsave(file.path(dirout, paste0('cyt_bool_', xvar,'_',cancertype,'_LH.pdf')))
      }
      
      res.corr.s[res.corr.s$var==xvar, paste0(cancertype,'_LH_meddelta')] <- med.delta[xvar] #CYT low-high (median of low minus median of high)
      res.corr.s[res.corr.s$var==xvar, paste0(cancertype,'_LH_pval')] <- st$p.value #CYT low-high (Mann-Whitney for continuous var, and chi-sq for categorical var)
    }, error=catch.errorfunc, finally={})
  }
  
  #append statistics to table
  if(!is.null(res.corr)){
    if(nrow(res.corr) < 1){
      res.corr <- res.corr.s
      rownames(res.corr) <- colnames(df.merged)
    } else
      res.corr <- cbind(res.corr, res.corr.s[match(res.corr$var, res.corr.s$var), 2:5])
  }
  
  #Hoeffding's D
  print('-- Collinearity - correlation using Hoeffding\'s D --')
  vc <- varclus(data.matrix(df.merged), similarity='hoeffding')
  print(vc, abbrev=FALSE) #print similarity matrix, etc.
  
  pdf(file.path(dirout, paste0('corHoefD_',cancertype,'.pdf')), width=12, height=8)
  plot(vc, abbrev=FALSE, cex=0.7)
  dev.off()
  
  #redundancy analysis
  print('-- Redundancy analysis --')
  print(redun(formula(paste0("~",paste0("",varsX, collapse="+"))), data=df.merged, nk=0))
  
  #write correlation results table for this specific cancer type
  write.csv(res.corr.s, row.names=FALSE, file=file.path(dirout,'summary_correlations.csv'))
  
  return(res.corr.s)
}

#analyze data
analyze.dsets.multi <- function(cancertype, res.model=NULL, dirout=dir.out, reanalyze=FALSE){
  #Inputs
  # cancertype
  # res.model: results data for all cancer types
  if(!file.exists(dirout)) dir.create(dirout)
  print.seg(sprintf('Analyzing data for cancer type %s | multivariate', cancertype))
  
  #check if this already has been analyzed, if so that reanalyze=F, then skip
  if(file.exists(file.path(dirout, 'summary_multivariate.csv')) && !reanalyze){
    print('There is already a summary_multivariate.csv file for this cancer type. Multivariate analysis should already be done. Skipping... ')
    return(TRUE)
  }
  
  #load data
  if(dbset.preprocess.std) load(file.path(dir.out, paste0(cancertype, '_prep_std.rda')))
  else if(dbset.preprocess.pca) load(file.path(dir.out, paste0(cancertype, '_prep_pca.rda')))
  else load(file.path(dir.out, paste0(cancertype, '_prep_base.rda')))
  
  if(nrow(df.merged) < 1){
    print("WARNING: input data.frame is empty...")
    return(FALSE)
  }
  
  #define variable names
  varsY <- c('cyt','GZMA','PRF1','cyt.bool','cyt.cat')
  varsX <- colnames(df.merged)
  varsX <- varsX[!(varsX %in% varsY)] #exclude the Y vars
  
  #----- multivariate
  print('Performing multivariate analysis...')
  #define variables list in models
  if(dbset.preprocess.pca){
    modelvars <- list(#all
      'PC1',
      'PC2',
      'PC3',
      'PC4',
      c('PC1', 'PC2'),
      c('PC1', 'PC3'),
      c('PC2', 'PC3'),
      c('PC1', 'PC2', 'PC3'),
      c('PC1', 'PC2', 'PC3', 'PC4')
    )
    modelslist <- list(
      # list('reg_rf_all', modelvars[[1]], 'reg', 'rf', FALSE),
      # list('reg_rf_mut', modelvars[[2]], 'reg', 'rf', FALSE),
      # list('reg_rf_NMD', modelvars[[3]], 'reg', 'rf', FALSE),
      # list('reg_rf_spgene', modelvars[[4]], 'reg', 'rf', FALSE),
      # list('reg_rf_mut+NMD', modelvars[[5]], 'reg', 'rf', FALSE),
      # list('reg_rf_mut+spgene', modelvars[[6]], 'reg', 'rf', FALSE),
      
      list('bin_rf_pc1', modelvars[[1]], 'bin', 'rf', FALSE),
      list('bin_rf_pc2', modelvars[[2]], 'bin', 'rf', FALSE),
      list('bin_rf_pc3', modelvars[[3]], 'bin', 'rf', FALSE),
      list('bin_rf_pc4', modelvars[[4]], 'bin', 'rf', FALSE),
      list('bin_rf_pc1+pc2', modelvars[[5]], 'bin', 'rf', FALSE),
      list('bin_rf_pc1+pc3', modelvars[[6]], 'bin', 'rf', FALSE),
      list('bin_rf_pc2+pc3', modelvars[[7]], 'bin', 'rf', FALSE),
      list('bin_rf_pc1+pc2+pc3', modelvars[[8]], 'bin', 'rf', FALSE),
      list('bin_rf_pc1+pc2+pc3+pc4', modelvars[[9]], 'bin', 'rf', FALSE)  )
  } else {
    modelvars <- list(#all
      varsX[!varsX %in% c('SNP','DEL','INS', 'nmdptc_min','nmdns_min','nmdfs_min')],
      
      #mutation
      c('silent','missense','nonstop','nonsense','frameshift','total'),
      
      #NMD
      varsX[!varsX %in% c('SNP','DEL','INS', 'nmdptc_min','nmdns_min','nmdfs_min',
                          'silent','missense','nonstop','nonsense','frameshift','total',
                          'msi',
                          colnames(df.merged)[grep('_cna$' ,colnames(df.merged))],
                          colnames(df.merged)[grep('_mut$' ,colnames(df.merged))] )],
      
      #genes mut/cna
      varsX[!varsX %in% c('SNP','DEL','INS', 'nmdptc_min','nmdns_min','nmdfs_min',
                          'silent','missense','nonstop','nonsense','frameshift','total',
                          'msi',
                          colnames(df.merged)[grep('^nmd' ,colnames(df.merged))] )],
      
      #msi
      c('msi'),
      
      
      #mutation+NMD
      varsX[!varsX %in% c('SNP','DEL','INS', 'nmdptc_min','nmdns_min','nmdfs_min',
                          'msi',
                          colnames(df.merged)[grep('_cna$' ,colnames(df.merged))],
                          colnames(df.merged)[grep('_mut$' ,colnames(df.merged))] )],
      
      #mutation+genes mut/cna
      varsX[!varsX %in% c('SNP','DEL','INS', 'nmdptc_min','nmdns_min','nmdfs_min',
                          'msi',
                          colnames(df.merged)[grep('^nmd' ,colnames(df.merged))] )],
      
      #mutation+msi
      varsX[!varsX %in% c('SNP','DEL','INS', 'nmdptc_min','nmdns_min','nmdfs_min',
                          colnames(df.merged)[grep('^nmd' ,colnames(df.merged))],
                          colnames(df.merged)[grep('_cna$' ,colnames(df.merged))],
                          colnames(df.merged)[grep('_mut$' ,colnames(df.merged))] )],
      
      #NMD + genes mut/cna
      varsX[!varsX %in% c('SNP','DEL','INS', 'nmdptc_min','nmdns_min','nmdfs_min',
                          'msi',
                          'silent','missense','nonstop','nonsense','frameshift','total' )], 
      
      #NMD + MSI
      varsX[!varsX %in% c('SNP','DEL','INS', 'nmdptc_min','nmdns_min','nmdfs_min',
                          'silent','missense','nonstop','nonsense','frameshift','total',
                          colnames(df.merged)[grep('_cna$' ,colnames(df.merged))],
                          colnames(df.merged)[grep('_mut$' ,colnames(df.merged))] )],
      
      #genes mut/cna + MSI
      varsX[!varsX %in% c('SNP','DEL','INS', 'nmdptc_min','nmdns_min','nmdfs_min',
                          colnames(df.merged)[grep('^nmd' ,colnames(df.merged))],
                          'silent','missense','nonstop','nonsense','frameshift','total' )]
      
    )
    
    #define models
    #model name, model fml (variables), regression/classification (reg/bin), 
    modelslist <- list(
      list('bin_rf_mut', modelvars[[2]], 'bin', 'rf', FALSE),
      list('bin_rf_NMD', modelvars[[3]], 'bin', 'rf', FALSE),
      list('bin_rf_msi', modelvars[[5]], 'bin', 'rf', FALSE),
      list('bin_rf_mut+NMD', modelvars[[6]], 'bin', 'rf', FALSE),
      list('bin_rf_mut+MSI', modelvars[[8]], 'bin', 'rf', FALSE),
      list('bin_rf_NMD+MSI', modelvars[[10]], 'bin', 'rf', FALSE) )
  }
  
  #build and analyze models
  res.model.s <- data.frame()
  for(moditem in modelslist){
    mod.name <- moditem[[1]]
    varnames <- moditem[[2]]
    mod.type <- moditem[[3]]
    mod.method <- moditem[[4]]
    mod.rd <- moditem[[5]]
    
    if(mod.type=='reg'){
      mod.data <- df.merged
      yname <- 'cyt'
    } else if(mod.type=='bin'){
      mod.data <- df.merged.bool
      yname <- 'cyt.cat'
    }
    
    #if varnames does not exist in data, remove
    varnames <- varnames[varnames %in% colnames(mod.data)]
    
    #check for missing data, only present warning if there is (dropsmall is FALSE here)
    handleMissingVals(mod.data, dropsmall=FALSE)
    
    res.model.s <- genanyzModel(varnames, yname, mod.data, res.model.s, mod.type, mod.method, mod.rd, cancertype, mod.name, dirout)
  }
  
  #append statistics to table (for all cancer types)
  if(!is.null(res.model)){
    if(nrow(res.model) < 1) res.model <- res.model.s
    else res.model <- cbind(res.model, res.model.s[-1])
  }
  
  #write results table for this specific cancer type
  write.csv(res.model.s, row.names=FALSE, file=file.path(dirout,'summary_multivariate.csv'))
  
  return(res.model.s)
}

analyze.dsets.clin <- function(cancertype, dirout=dir.out){
  if(!file.exists(dirout)) dir.create(dirout)
  
  #parameters for library
  outpref.sink <<- TRUE
  outpref.verbose <<- FALSE
  outpref.plot <<- TRUE
  outpref.file <<- TRUE
  param.timeUnit <<- 'Year'
  save.plot <<- TRUE
  outpref.diagnostics <<- FALSE
  
  #load data
  load(file.path(dir.out, paste0(cancertype, '_prep_base.rda')))
  
  if(nrow(df.merged.clin) < 1){
    print("WARNING: input data.frame is empty...")
    return(FALSE)
  }
  
  #define variables to look at
  varnames <- c(colnames(df.merged), 'tnm_stage', 'gender','race','age')
  
  #determine which variables are factor/numeric
  var.cat <- retrieveVarTypeNames(df.merged.clin[, colnames(df.merged.clin) %in% varnames], is.factor)
  var.cont <- retrieveVarTypeNames(df.merged.clin[, colnames(df.merged.clin) %in% varnames], is.numeric)
  
  var.cat <- c(var.cat, var.cont[grepl('\\.bool$',var.cont)]) #take care of .bool labels, these are categorical, not continuous
  var.cont <- var.cont[!grepl('\\.bool$',var.cont)] #and remove them from the continuous list here
  
  if('nmdns_mut.bool' %in% colnames(df.merged.clin)) df.merged.clin$nmdns_mut.bool <- factor(df.merged.clin$nmdns_mut.bool)
  if('nmdfs_mut.bool' %in% colnames(df.merged.clin)) df.merged.clin$nmdfs_mut.bool <- factor(df.merged.clin$nmdfs_mut.bool)
  if('nmdptc_mut.bool' %in% colnames(df.merged.clin)) df.merged.clin$nmdptc_mut.bool <- factor(df.merged.clin$nmdptc_mut.bool)
  
  #univariate categorical KM
  #create continuous variables that will be categorized for Kaplan-Meier
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
  
  outdir <- file.path(dirout, 'uni.km')
  if(!file.exists(outdir)) dir.create(outdir)
  dfsummary <- genStats.uni_km(df.merged.clin, c('1', c(var.cont.km, var.cat)), outdir=outdir)
  
  #univariate Cox regression
  outdir <- file.path(dirout, 'uni.cox')
  if(!file.exists(dirout)) dir.create(dirout)
  outpref.plot <<- FALSE
  dfsummary <- genStats.uni_cox(df.merged.clin, union(union(var.cat, var.cont),var.cont.km), outdir=outdir)
  
  #multivariate Cox regression - controlling for baselines: age, gender, and tnm_stage
  outdir <- file.path(dirout, 'uni.cox.withctrl')
  if(!file.exists(dirout)) dir.create(dirout)
  outpref.plot <<- FALSE
  varnames <- union(union(var.cat, var.cont),var.cont.km)
  varnames <- varnames[!varnames %in% c('age','gender','tnm_stage','race')] #remove the variables to control
  varnames.ctrl <- c('age','gender','tnm_stage')
  varnames.ctrl <- varnames.ctrl[varnames.ctrl %in% colnames(df.merged.clin)]
  dfsummary <- genStats.uni_cox.withctrl(df.merged.clin, varnames, varnames.ctrl, outdir=outdir)
}

analyze.dsets.clin2 <- function(cancertype, dirout=dir.out){
  #survival analyses, this extends analyze.dsets.clin to include TMB and PDL1
  #Requires: cancertypes
  if(!file.exists(dirout)) dir.create(dirout)
  
  #parameters for library
  outpref.sink <<- TRUE
  outpref.verbose <<- FALSE
  outpref.plot <<- TRUE
  outpref.file <<- TRUE
  param.timeUnit <<- 'Year'
  save.plot <<- TRUE
  outpref.diagnostics <<- FALSE
  
  #load data
  load(file.path(dir.out, paste0(cancertype, '_prep_base.rda')))
  if(!exists('df.merged.clin')) return(NULL)
  
  #append TMB and PDL1 information, for PDL1 we first need to the gene expression of it
  if(cancertype == 'pancancer')
    cancertypes.PDL1 <- cancertypes[cancertypes != 'pancancer']
  else
    cancertypes.PDL1 <- cancertypes
  
  df.rnaseq.PDL1 <- retrieveRNASeqGene(cancertypes.PDL1,'CD274')
  df.merged.clin <- append.clinVars(df.merged.clin, df.rnaseq.PDL1)
  
  if(nrow(df.merged.clin) < 1){
    print("WARNING: input data.frame is empty...")
    return(FALSE)
  }
  
  #corrections - group low levels with very few patients
  if(cancertype=='BLCA'){
    #BLCA only has one patient for stage I, group this with stage II together
    df.merged.clin$tnm_stage <- as.character(df.merged.clin$tnm_stage)
    df.merged.clin$tnm_stage[df.merged.clin$tnm_stage=='I'] <- 'I/II'
    df.merged.clin$tnm_stage[df.merged.clin$tnm_stage=='II'] <- 'I/II'
    df.merged.clin$tnm_stage <- factor(df.merged.clin$tnm_stage, levels=c('I/II','III','IV'))
  }
  
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
  
  #univariate categorical KM
  #create continuous variables that will be categorized for Kaplan-Meier
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
  
  outdir <- file.path(dirout, 'uni.km2')
  if(!file.exists(outdir)) dir.create(outdir)
  dfsummary <- genStats.uni_km(df.merged.clin, c('1', c(var.cont.km, var.cat)), outdir=outdir)
  
  #univariate Cox regression
  outdir <- file.path(dirout, 'uni.cox2')
  if(!file.exists(dirout)) dir.create(dirout)
  outpref.plot <<- FALSE
  dfsummary <- genStats.uni_cox(df.merged.clin, union(union(var.cat, var.cont),var.cont.km), outdir=outdir)
  
  #multivariate Cox regression - set up
  outpref.plot <<- FALSE
  varnames <- union(union(var.cat, var.cont),var.cont.km)
  varnames <- varnames[!varnames %in% c('age.cat','gender','tnm_stage','race')] #remove the variables to control
  
  #multivariate Cox regression - controlling for baselines: age, gender, and tnm_stage
  outdir <- file.path(dirout, 'uni.cox.withctrl2')
  if(!file.exists(dirout)) dir.create(dirout)
  varnames.ctrl <- c('age.cat','gender','tnm_stage')
  varnames.ctrl <- varnames.ctrl[varnames.ctrl %in% colnames(df.merged.clin)]
  genStats.uni_cox.withctrl(df.merged.clin, varnames, varnames.ctrl, outdir=outdir, fname='univariate_cox_withctrl.csv')

  #multivariate Cox regression - controlling for baselines: age, gender, and tnm_stage, TMB.cat (>median), PDL1.cat (low/med/hight)
  outdir <- file.path(dirout, 'uni.cox.withctrl2')
  if(!file.exists(dirout)) dir.create(dirout)
  varnames.ctrl <- c('age.cat','gender','tnm_stage','TMB.cat','PDL1.cat')
  varnames.ctrl <- varnames.ctrl[varnames.ctrl %in% colnames(df.merged.clin)]
  genStats.uni_cox.withctrl(df.merged.clin, varnames, varnames.ctrl, outdir=outdir, fname='univariate_cox_withctrl_TMBPDL1cat.csv')

}

summarize.combineResults <- function(cancertypes){
  #summarize results - aggregated
  if(!file.exists(file.path(dir.out, 'summaries'))) dir.create(file.path(dir.out, 'summaries'))
  res.corr <- data.frame() #univariate stat results
  res.model <- data.frame() #multivariate model results
  
  for(cancertype in cancertypes){
    #append res.corr
    f <- file.path(dir.out, paste0('_',cancertype), 'summary_correlations.csv')
    if(file.exists(f)){
      res.corr.s <- read.csv(f)
      if(nrow(res.corr) < 1) res.corr <- res.corr.s
      else res.corr <- cbind(res.corr, res.corr.s[match(res.corr$var, res.corr.s$var), 2:5])
    }
    
    f <- file.path(dir.out, paste0('_',cancertype), 'summary_multivariate.csv')
    if(file.exists(f)){
      #append res.model
      res.model.s <- read.csv(f)
      if(nrow(res.model) < 1) res.model <- res.model.s
      else res.model <- cbind(res.model, res.model.s[-1])
    }
  }
  
  write.csv(res.corr, row.names=FALSE, file=file.path(dir.out, 'summaries', 'univariate.csv')) #univariate statistics
  write.csv(res.model, row.names=FALSE, file=file.path(dir.out, 'summaries', 'multivariate.csv')) #multivariate statistics
  
  return(list(res.corr, res.model))
}

analyze.res.model <- function(cancertype, res.model){
  #parse and analyze the model results table (generated by analyze.dsets.multi)
  #Requires: dir.out
  #SE is calculated from SD for CV, and we know we used a 10-fold CV
  print(sprintf('Summarizing %s..',cancertype))
  if(!file.exists(file.path(dir.out, 'summaries'))) dir.create(file.path(dir.out, 'summaries'))
  vals.auctval <-  sub(';.*','',res.model[[paste0(cancertype,'_full__vals')]])
  vals.auctval <- sub('.*:', '', vals.auctval)
  vals.auccv <-  sub(';auc\\.cv\\.sd:.*','',res.model[[paste0(cancertype,'_full__vals')]])
  vals.auccv <- sub('.*auc\\.cv:', '', vals.auccv)
  vals.auccvsd <-  sub(';auc\\.test:.*','',res.model[[paste0(cancertype,'_full__vals')]])
  vals.auccvsd <- sub('.*auc\\.cv\\.sd:', '', vals.auccvsd)
  vals.oob <-  sub('.*oob:','',res.model[[paste0(cancertype,'_full__vals')]])
  df <- data.frame(mod=sub('.*_','',res.model$model_name),
                   auc.tval=as.numeric(vals.auctval),
                   auc.cv=as.numeric(vals.auccv),
                   auc.cv.sd=as.numeric(vals.auccvsd),
                   auc.cv.se=as.numeric(vals.auccvsd)/sqrt(10), #NOTE: for now, we know we used 10-fold CV, so '10' is hard-coded
                   oob=as.numeric(vals.oob))
  df$mod <- factor(df$mod, df$mod)
  
  df <- df[complete.cases(df),] #only keep non-missing data
  
  if(nrow(df) < 1){
    print('Resulting data.frame is empty, skipping...')
    return(FALSE)
  }
  
  #deal with MSI, if MSI is empty, don't need to include models that are also MSI, redundant
  if(!'msi' %in% df$mod){
    df <- df[!grepl("+MSI$",df$mod),]
  }
  
  #produce the plots without spgene and all
  df <- df[!grepl("spgene",df$mod),]
  df <- df[!grepl("all",df$mod),]
  
  #OOB
  ggplot(data=df, aes(x=mod, y=oob)) + geom_bar(stat='identity') + ylab('Out-of-bag error') + xlab('Models')
  ggsave(file.path(dir.out, 'summaries', paste0('oob_bar_', cancertype, '.pdf')))
  
  #AUROC tval
  ggplot(data=df, aes(x=mod, y=auc.tval)) + geom_bar(stat='identity') + ylab('AUROC (train)') + xlab('Models')
  ggsave(file.path(dir.out, 'summaries', paste0('auroc_tval_bar_', cancertype, '.pdf')))
  
  #AUROC CV +/- SE(CV)
  ylimL <- floor(min(df$auc.cv-df$auc.cv.sd)*10)/10
  ylimH <- ceiling(max(df$auc.cv+df$auc.cv.sd)*10)/10
  ggplot(data=df, aes(x=mod, y=auc.cv)) + 
    geom_bar(stat='identity') + geom_errorbar(aes(ymin=auc.cv-auc.cv.se, ymax=auc.cv+auc.cv.se), width=.1) +
    ylab('AUROC (cross-validation)') + xlab('Models') + 
    coord_cartesian(ylim=c(0.4,0.9))
    #coord_cartesian(ylim=c(ylimL,ylimH))
  ggsave(file.path(dir.out, 'summaries', paste0('auroc_bar_', cancertype, '.pdf')))
}

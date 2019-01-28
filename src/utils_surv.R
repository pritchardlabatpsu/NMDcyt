#utils_surv: survival helper methods
# Statistical analyses of dataset and models
pacman::p_load(plyr, dplyr, broom)
pacman::p_load(survival, survminer, rms, relsurv, epiR, pec, ggplot2, GGally)

#--------- General plot methods ---------
saveFig <- function(figdir, figname, savePlots, height=7, width=7){
  #save figure given directory and filename
  if(savePlots)
    dev.print(pdf, height=height, width=width, filename=file.path(figdir, paste0(figname,'.pdf')))
}

openFigDev <- function(figdir, figname, savePlots, height=7, width=7){
  #opens pdf device
  if(savePlots){
    pdf(file.path(figdir, paste0(figname,'.pdf')), width=width, height=height)
    dev.control("enable") #turn on recording so dev can be copied
  }
}

closeFigDev <- function(savePlots){
  if(savePlots) dev.off()
}

format.varname <- function(varname){
  varname <- sub(' ','_',varname)
  varname <- sub('/','\\.',varname)
  return(varname)
}

#--------- public: general dataset statistics - individual analyses ---------
genStats.uni_km <- function(df, varnames=NULL, subdir='', outdir=NULL){
  #KM estimate, for categorical variables and intercept only model
  if(is.null(outdir)) outdir <- file.path(dir.out, subdir)
  if(!file.exists(outdir)) dir.create(outdir)
  
  if(is.null(varnames)){
    varnames <- retrieveVarTypeNames(df, is.factor)
    varnames <- c('1', varnames)
  }
  
  dfsummary <- data.frame()
  for(varname in varnames){
    if(!varname %in% colnames(df)){
      print(sprintf('%s does not exist in data.frame, skipping...', varname))
      next
    }
    
    #if there is only one level (based on unique) and it is not ~1 full model, can skip - this is same as ~1 full model
    val <- df[[varname]]
    if(varname!=1 && length(unique(val[!is.na(val)]))<2){
      print(sprintf('There are no more than >1 level for this variable %s, skipping...', varname))
      next
    }
    
    smod <- getUnivariate.mod(df, varname, 'km') #create Kaplan-Meier estimatess
    if(outpref.verbose) genStats.survKM.verbose(smod, df, varname) #statistics: verbose
    if(outpref.plot) genStats.survKM.plots(smod, df, varname, subdir, outdir) #statistics: plots
    if(outpref.file){
      res <- genStats.survKM.table(smod, df, subdir, outdir, fname=NULL) #statistics: table
      if(nrow(dfsummary) < 1) dfsummary <- res else dfsummary <- rbind(dfsummary, res)
    }
  }
  if(outpref.file) write.csv(dfsummary, file=file.path(outdir,'univariate_summary_km.csv'), row.names=FALSE)
  
  return(dfsummary)
}

genStats.uni_cox <- function(df, varnames=NULL, subdir='', outdir=NULL){
  #Cox regression, for all variables (categorical+continuous)
  if(is.null(outdir)) outdir <- file.path(dir.out, subdir)
  if(!file.exists(outdir)) dir.create(outdir)
  
  if(is.null(varnames)){
    varnames.cat <- retrieveVarTypeNames(df, is.factor)
    varnames.cont <- retrieveVarTypeNames(df, is.numeric)
    varnames <- c(varnames.cat, varnames.cont)
  }
  
  dfsummary <- data.frame()
  for(varname in varnames){
    if(!varname %in% colnames(df)){
      print(sprintf('%s does not exist in data.frame, skipping...', varname))
      next
    }
    
    #if there is only one level (based on unique), can skip
    val <- df[[varname]]
    if(length(unique(val[!is.na(val)]))<2) next
    
    modobj <- getUnivariate.mod(df, varname, 'cox') #create cox model
    if(outpref.verbose) genStats.survCox.verbose(modobj)
    if(outpref.plot) genStats.survCox.plots(modobj, paste0('uni_',varname), df, model.data.km=df, subdir=subdir, outdir=outdir)
    if(outpref.file){
      res <- genStats.survCox.table(modobj, fname=NULL)
      if(nrow(dfsummary) < 1) dfsummary <- res else dfsummary <- rbind(dfsummary, res)
    }
  }
  if(outpref.file) write.csv(dfsummary, file=file.path(outdir,'univariate_summary_cox.csv'), row.names=FALSE)
  
  return(dfsummary)
}

genStats.uni_cox.withctrl <- function(df, varnames=NULL, varnames.ctrl=c(), subdir='', outdir=NULL, fname='univariate_summary_cox_withctrl.csv'){
  #varnames.ctrl refers to baseline
  if(is.null(outdir)) outdir <- file.path(dir.out, subdir)
  if(!file.exists(outdir)) dir.create(outdir)
  
  if(length(varnames.ctrl) < 1){
    print('Baseline variables are empty, there is nothing to control for, no need to run multivariate Cox model here...')
    return(c())
  }
  
  if(is.null(varnames)){
    varnames.cat <- retrieveVarTypeNames(df, is.factor)
    varnames.cont <- retrieveVarTypeNames(df, is.numeric)
    varnames <- c(varnames.cat, varnames.cont)
  }
  
  #check over the varnames.ctrl (baseline variables)
  ckUniqueInc <- function(varname){
    #return TRUE if varname is good and will be included
    if(!varname %in% colnames(df)){
      FALSE
    } else{
      val <- df[[varname]]
      length(unique(val[!is.na(val)]))>=2
    }
  }
  varnames.ctrl <- varnames.ctrl[sapply(varnames.ctrl, ckUniqueInc)]
  
  if(length(varnames.ctrl) < 1){
    print('Baseline variables are empty after filtering, there is nothing to control for, no need to run multivariate Cox model here...')
    return(c())
  }
  
  dfsummary <- data.frame()
  for(varname in varnames){
    if(!varname %in% colnames(df)){
      print(sprintf('%s does not exist in data.frame, skipping...', varname))
      next
    }
    
    #if there is only one level (based on unique), can skip
    val <- df[[varname]]
    if(length(unique(val[!is.na(val)]))<2) next
    
    modobj <- getMultivariate.mod(df, varname, varnames.ctrl, 'cox') #create cox model
    if(outpref.verbose) genStats.survCox.verbose(modobj)
    if(outpref.plot) genStats.survCox.plots(modobj, paste0('uni_',varname), df, model.data.km=df, subdir=subdir, outdir=outdir)
    if(outpref.file){
      res <- genStats.survCox.table(modobj, fname=NULL)
      if(nrow(dfsummary) < 1) dfsummary <- res else dfsummary <- rbind(dfsummary, res)
    }
  }
  if(outpref.file) write.csv(dfsummary, file=file.path(outdir, fname), row.names=FALSE)
  
  return(dfsummary)
}

#--------- public: general survival univariate statistics - individual analyses ---------
getUnivariate.mod <- function(model.data, varname, modtype){
  #Generates univariate model
  #Requires: param.timeUnit, save.plot, dir.out, outpref.plot
  
  if(modtype=='km'){
    #Kaplan-Meier model
    fmla <- paste0("Surv(duration, event)~", varname)
    print(paste0("Model: ",fmla))
    eval(parse(text=paste0("modobj <- survfit(", fmla, ",data=model.data, conf.type='log', type='kaplan-meier')")))
  } else if(modtype=='cox'){
    #Cox regression
    fmla <- paste0("Surv(duration, event)~", varname)
    print(paste0("Model: ",fmla))
    eval(parse(text=paste0("modobj <- coxph(", fmla, ",data=model.data)")))
  }
  
  return(modobj)
}

getMultivariate.mod <- function(model.data, varname, varnames.ctrl, modtype){
  #Requires: param.timeUnit, save.plot, dir.out, outpref.plot
  #Only for cox regression, control for 
  
  if(modtype=='km'){
    modobj <- NULL
  } else if(modtype=='cox'){
    #Cox regression
    fmla <- paste0("Surv(duration, event)~", paste(c(varname,varnames.ctrl),collapse="+"))
    print(paste0("Model: ",fmla))
    eval(parse(text=paste0("modobj <- coxph(", fmla, ",data=model.data)")))
  }
  
  return(modobj)
}

genStats.survKM.verbose <- function(modobj, model.data, varname){
  #Generate Kaplan-Meier Verbose summary statistics
  #Log-rank test (Mantel-Cox test)
  if(varname!='1' && length(unique(model.data[,varname])) > 1)
    print(survdiff(as.formula(modobj$call$formula), data=model.data)) #default rho=0 -> logrank test, rho=1 -> Gehan-Wilcoxon test
  
  #summarize KM model
  print(modobj)
}

genStats.survKM.plots <- function(modobj, model.data, varname, subdir='', outdir=NULL, conf.int=FALSE){
  #Generate Kaplan-Meier plots
  genPlots.km(modobj, model.data, varname, subdir, outdir, conf.int=conf.int)
}

genStats.survKM.table <- function(modobj, model.data, subdir='', outdir=NULL, fname='stats_km_table'){
  #Generate Kaplan-Meier summary table
  if(is.null(outdir)) outdir <- file.path(dir.out, subdir)
  if(!file.exists(outdir)) dir.create(outdir)
  dfsummary <- getSurv.medprob(modobj, refdata=model.data)
  
  if(!is.null(fname)){
    #if fname is not null, write summary table to csv
    if(is.null(outdir)) outdir <- file.path(dir.out, subdir)
    if(!file.exists(outdir)) dir.create(outdir)
    write.csv(dfsummary, file=file.path(outdir, paste0(fname,'.csv')))
  }
  
  return(dfsummary)
}

#--------- public: general survival multivariate statistics - individual analyses ---------
genStats.survCox.verbose <- function(modobj){
  #Generate Cox model verbose summary statistics
  modph.summary <- summary(modobj)
  print(modph.summary)
  
  if(outpref.diagnostics){
    #proportionality test
    model.phzph <- cox.zph(model.ph)
    print(model.phzph)
  }
}

genStats.survCox.table <- function(modobj, subdir='', outdir=NULL, fname='stats_cox_table'){
  #Generate Cox model summary table
  modph.summary <- summary(modobj)
  nvar <- length(attr(modobj$terms,'term.labels'))
  varnames <- attr(modobj$terms,'term.labels')
  
  #construct table
  res.var <- cbind(as.data.frame(modph.summary$coefficients),
                   as.data.frame(modph.summary$conf.int)[3:4])
  res.var <- cbind(rownames(res.var), res.var)
  rownames(res.var) <- NULL
  colnames(res.var) <- c('var','coef', 'hr', 'se_coef', 'z', 'pval', 'hr_ci95L', 'hr_ci95H')
  res.var[,c('hr','hr_ci95L','hr_ci95H')] <- apply(res.var[,c('hr','hr_ci95L','hr_ci95H')], 2,
                                                   function(x)round(x,1))
  res.var$pval <- formatC(res.var$pval, format = "e", digits = 2)
  res.var$hr_ci95 <- paste0(res.var$hr_ci95L, ' to ', res.var$hr_ci95H)
  res.var <- res.var[,c('var','hr','pval','hr_ci95')]
  
  dfsummary <- res.var #summary table to be returned
  
  if(!is.null(fname)){
    #if fname is not null, write summary table to csv
    if(is.null(outdir)) outdir <- file.path(dir.out, subdir)
    if(!file.exists(outdir)) dir.create(outdir)
    write.csv(dfsummary, file=file.path(outdir, paste0(fname,'.csv')))
  }
  return(dfsummary)
}

#generate Cox model plots
genStats.survCox.plots <- function(modobj, modelname, model.data, newdata=NULL, model.data.km=NULL, subdir='', outdir=NULL){
  #generate plots
  genPlots.cox(modobj, modelname, model.data, newdata, model.data.km, subdir, outdir)
}

#--------- private: general dataset statistics helpers/individual analyses ---------
retrieveVarTypeNames <- function(df, func){
  #Retrieve variable names in given data frame (df) that satisfy @func
  varnames <- sapply(df[,!(colnames(df) %in% c('duration','event', 'tstart', 'tstop'))], func)
  varnames <- names(varnames)[varnames==TRUE]
  
  if(length(varnames) == 1 && varnames == "") varnames <- NULL
  
  return(varnames)
}

#--------- private: general survival statistics ---------
getSurv.med.pa <- function(model, newdata){
  #Retrieve median survival from parametric model
  val <- predict(model, newdata=newdata, type='quantile', p=0.5, se=TRUE) #median survival
  return(data.frame(medsurv=val$fit, ci95L=val$fit-1.96*val$se.fit, ci95U=val$fit+1.96*val$se.fi))
}

getSurv.med <- function(modelobj, newdata=NULL){
  #Retrieve median survival from survfit object
  #modelobj can be either survfit object or obj (e.g. coxph) as input to survfit
  
  sfit <- getSurv.helper(modelobj, newdata)
  if(is.null(sfit)) return(NULL)
  
  # val <- read.table(textConnection(capture.output(sfit)),skip=4,header=TRUE)
  # if('strata' %in% names(sfit)){ strataval <- rownames(val) } else { strataval <- '-' }
  # res <- data.frame(strata=strataval, n=val$n, n_event=val$events, medsurv=val$median, ci95L=val$X0.95LCL, ci95H=val$X0.95UCL)
  
  sfit.tbl <- summary(sfit)$table
  if('strata' %in% names(sfit)){
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
  res[,c('medsurv', 'ci95L', 'ci95H')] <- apply(res[,c('medsurv', 'ci95L', 'ci95H')], 2,
                                                function(x)round(x,1))
  res$ci95 <- paste0(res$ci95L,' to ', res$ci95H)
  return(res)
}

getSurv.prob <- function(modelobj, surv.time=NULL, newdata=NULL,refdata=NULL){
  #Retrieve survival probabilities from survfit object, at the specified time points
  #Returns observed, and if specified, expected and relative survival
  #Inputs
  # modelobj: can be either survfit object or obj (e.g. coxph) as input to survfit
  # surv.time: timepoints
  # newdata: newdata to applied with modelobj
  # survRel: calculate relative survival
  # survAgeStd: calculate age-standardized survival
  # refdata: reference data to get slice data, to use for calculations of expected survival, and age-standardization
  
  sfit <- getSurv.helper(modelobj, newdata)
  if(is.null(sfit)) return(NULL)
  
  #if the timepoints is NULL, define default
  if(is.null(surv.time)) surv.time <- c(0.5, 1, 2, 5) #in years
  
  val <- summary(sfit, time=surv.time, extend=TRUE)
  #note: extend is TRUE, so will print value at all times, even if there are no subjects; without extend=T, in this scenario it would throw an error
  if('strata' %in% names(val)){ strataval <- val['strata'] } else { strataval <- '-' }
  res <- data.frame(strata=strataval, time=val$time,
                    spO=val$surv, spO_se=val$std.err,
                    spO_ci95L=val$lower, spO_ci95H=val$upper)
  colnames(res)[1] <- 'strata'
  
  #rounding
  col2round <- c('spO', 'spO_se', 'spO_ci95L', 'spO_ci95H')
  res[, col2round] <- apply(res[, col2round], 2, function(x)round(x,2))
  res$spO_ci95 <- paste0(res$spO_ci95L,' to ', res$spO_ci95H)
  
  return(res)
  
}

getSurv.helper <- function(modelobj, newdata=NULL){
  #Helper function
  #c the input object, and return a survfit object
  if('coxph' %in% class(modelobj)){
    if(is.null(newdata)) sfit <- survfit(modelobj)
    else sfit <- survfit(modelobj, newdata=newdata)
  } else if('survfit' %in% class(modelobj)){
    sfit <- modelobj
  } else {
    warning('While retrieving median survival, the given input to function is neither a coxph object nor survfit object...')
    return(NULL)
  }
  
  return(sfit)
}

getSurv.medprob <- function(modelobj, surv.time=NULL, newdata=NULL, refdata=NULL){
  #Combined call to getSurv.med() and getSurv.prob()
  #Gets the median survival, and survival probabilities at specified time points
  #Inputs
  # modelobj: can be either survfit object or obj (e.g. coxph) as input to survfit
  # surv.time: timepoints
  # newdata: newdata to applied with modelobj
  # refdata: pass onto getSurv.prob, reference data to be used to get cohort to derive match cohort for expected survival
  sfit <- getSurv.helper(modelobj, newdata)
  if(is.null(sfit)) return(NULL)
  
  #if the timepoints is NULL, define default
  if(is.null(surv.time)) surv.time <- c(0.5, 1, 2, 5) #in years
  
  #summary for median survival
  res.survmed <- getSurv.med(sfit)
  colnames(res.survmed) <- c('strata', 'n', 'n_event', 'medsurv', 'medsurv_ci95L', 'medsurv_ci95H', 'medsurv_ci95')
  res.survmed <- res.survmed[,c('strata', 'n', 'n_event', 'medsurv', 'medsurv_ci95')]
  rownames(res.survmed) <- res.survmed$strata
  res.survmed <- res.survmed %>% select(-strata)
  
  #summary for survival probabilities (at time points)
  res.prob <- getSurv.prob(sfit, surv.time, refdata=refdata)
  res.prob$spO_txt <- sprintf('%0.2f (%s)',
                              res.prob$spO,
                              res.prob$spO_ci95)
  
  #transpose survival prob table
  res.prob.transpose <- function(res.prob, txt.col, col.prefix){
    res.prob.t <- as.data.frame(matrix(NA,
                                       nrow = length(unique(res.prob$strata)),
                                       ncol = length(unique(res.prob$time)))
    )
    rownames(res.prob.t) <- unique(res.prob$strata)
    colnames(res.prob.t) <- unique(res.prob$time)
    for(n in 1:nrow(res.prob)){ #this is not the most elegant way, but table is small
      sel.r <- rownames(res.prob.t)==res.prob[n,'strata']
      sel.c <- colnames(res.prob.t)==res.prob[n,'time']
      res.prob.t[sel.r, sel.c] <- res.prob[n, txt.col]
    }
    
    colnames(res.prob.t) <- sprintf('%s_%syr', col.prefix, colnames(res.prob.t)) #rename colname
    
    return(res.prob.t)
  }
  res.prob.t <- res.prob.transpose(res.prob, 'spO_txt', 'survprobObs')
  
  #combine tables
  dfsummary <- merge(res.survmed, res.prob.t, by='row.names')
  colnames(dfsummary)[1] <- 'strata'
  
  #write summary table to csv (if applicable) and return table
  return(dfsummary)
}

#--------- private: survival plot generations ---------
genPlots.km <- function(sfit, model.data, varname, subdir='', outdir=NULL, conf.int=FALSE){
  #Generate S(t), H(t), h(t), etc plots given survfit
  #Inputs:
  # sfit: survfit object (kmf)
  # model.data: dataset
  # varname: variable name (if value is "1", this is the baseline model with intercept only)
  # subdir: subdirectory
  # outdir: full output directory, if defined, supersedes value in subdir
  #Requires: param.timeUnit, save.plot, dir.out
  
  if(is.null(outdir)) outdir <- file.path(dir.out, subdir)
  if(!file.exists(outdir)) dir.create(outdir)
  
  #initialization
  if(varname == '1'){
    varname <- 'intercept'
    varname.unique <- 1 #placeholder, just need this to be a length of 1
    varname.unique.s <- 1 #placeholder
  } else {
    varname.unique <- unique(model.data[,varname])
    varname.unique.s <- summary(model.data[varname])
  }
  ptitle <- sprintf('Kaplan-Meier Estimate of S(t) | %s', varname)
  
  #basic kaplan-meier curve
  openFigDev(outdir, paste0('uni_cat_',format.varname(varname),'_surv'), save.plot)
  plot(sfit, xlab=sprintf("Time (%s)",param.timeUnit), ylab='Survival Probability',
       conf.int=conf.int, mark.time=TRUE,
       main=ptitle, lwd=2, col=1:length(varname.unique))
  if(varname!='overall')
    legend('topright',legend=varname.unique.s, col=1:length(varname.unique), bty='n', horiz=FALSE, lty=1, lwd=2)
  closeFigDev(save.plot)
  
}

genPlots.cox <- function(model.ph, modelname, model.data, newdata=NULL, model.data.km=NULL, subdir='', outdir=NULL){
  if(is.null(outdir)) outdir <- file.path(dir.out, subdir)
  if(!file.exists(outdir)) dir.create(outdir)
  
  nvar <- length(attr(model.ph$terms,'term.labels'))
  varnames <- attr(model.ph$terms,'term.labels')
  modelname <- format.varname(modelname)
  
  #survfit and survival plots
  r <- tryCatch({
    if(is.null(newdata)) sf <- survfit(model.ph) #survival using the mean of all the covariates
    else sf <- survfit(model.ph, newdata=newdata, id=id) #requires that the id column is also specified
    
    #plot estimated S(t) from Cox regression model using mean covariates
    openFigDev(outdir, paste0(modelname,'_surv'), save.plot)
    plot(sf, xlab=sprintf("Time (%s)",param.timeUnit), ylab='Survival Probability',
         main='S(t) estimate (based on mean val of covariates)', lwd=2)
    closeFigDev(save.plot)
    
    # #plot estimated S(t) from Cox regression model using mean covariates
    #ggsurvplot can't handle some of the survfit models (if more than one curve), comment this out for now
    # openFigDev(outdir, paste0(modelname,'_surv'), save.plot)
    # print(ggsurvplot(sf, data=model.data, conf.int=TRUE, pval=TRUE, risk.table=TRUE,
    #                  main=paste0('S(t) estimate (based on mean val of covariates) | ', modelname),
    #                  xlab=sprintf("Time (%s)",param.timeUnit)))
    # closeFigDev(save.plot)
    
    if(!is.null(model.data.km)){
      #survival plot overlay with KM
      model.kmf <- survfit(Surv(duration,event) ~ 1, data=model.data.km, conf.type='log', type='kaplan-meier')
      
      openFigDev(outdir, paste0(modelname,'_surv_withKM'), save.plot)
      plot(model.kmf, xlab=sprintf("Time (%s)",param.timeUnit), ylab='Survival Probability', main='S(t) (KM and Cox overlay)',
           lwd=1, col='gray')
      lines(sf, lwd=1, col='red')
      legend('topright',
             c(paste0('cox model | ', modelname), 'Kaplan-Meier estimate'),
             col=c('red','gray'), bty='n', lty=1, cex=1)
      closeFigDev(save.plot)
    }
  }, error=catch.errorfunc, finally={})
}

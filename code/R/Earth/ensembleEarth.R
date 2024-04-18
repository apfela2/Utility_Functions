# Do ensemble CV on earth model


ensembleEarth <- function(dat, y, x, pmethod, nfold, ncross, stratify, trace, keepxy = TRUE, degree, glm = list(family = binomial),
                          run.parallel = FALSE, num.boot = 10, outer.loop.cv.nfolds = NA) {
  
  
  num.cors = 1
  #setting up parallelization
  if (is.logical(run.parallel) & run.parallel == TRUE) {
    num.cors = max(1, (detectCores() - 1)) #registerDoMC(cores = min(n.folds, detectCores() -3))
  } else if (is.logical(run.parallel) & run.parallel == FALSE) {
    num.cors = 1
  } else if (is.numeric(run.parallel)) {
    num.cors = run.parallel
  }
  
  #setting up folds for outer cv loop, if supplied
  rep.row.idxs <- NA
  if (is.list(outer.loop.cv.nfolds)) {
    rep.row.idxs <- outer.loop.cv.nfolds
  } else if(!is.na(outer.loop.cv.nfolds)) {
    #create row idxs for repeated cv
    rep.row.idxs <- foldsToIdxs(makeFolds(x = y, ntimes = num.boot, 
                                          nfolds = outer.loop.cv.nfolds, stratified.cv = stratified.cv))
    
    #foldsToIdxs and makeFolds are helper functions futher down in this file
    
    num.boot = num.boot*outer.loop.cv.nfolds
  }
  num.boot <- max(length(rep.row.idxs), num.boot)  
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #the stuff below happens for each bootstrapped iteration
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#  
  model.ens <- mclapply(seq_len(num.boot), function(nb) {
    # this.row.idxs <- NA
    # 
    # if (length(rep.row.idxs) == 1 && is.na(rep.row.idxs)) {  
    #   #we are in scenario where outer loop of cv is not being done
    #   if (is.na(as.numeric(sample.perc))) {
    #     #if resample precentage isn't specified, do bootstrapping
    #     this.row.idxs <- sample(seq_len(nrow(dset)), replace = TRUE)     
    #   } else {
    #     #if resample percentage is specified, do repeated random samples
    #     this.row.idxs <- sample(seq_len(nrow(dset)), size = round(nrow(dset)*sample.perc), replace = FALSE) 
    #   }
    # } else {
    this.row.idxs <- rep.row.idxs[[nb]]
    # }
    
    #
    # if (num.boot == 1) #if not doing any bootstrapping, use all samples
    #   this.row.idxs <- seq_len(nrow(dset))
    
    #we are running repeated glmnet for each bootstrap. setting run.parallel to FALSE because
    #we have already parallelized in the outside layer
    # if (!is.matrix(y) & !is.data.frame(y)) 
    #   y <- as.matrix(y)
    
    models <- earth(y = y[this.row.idxs], x = x[this.row.idxs, , drop = FALSE], pmethod = pmethod, nfold = nfold,
                    ncross = ncross, stratify = stratify, trace = trace, keepxy = keepxy,
                    degree = degree, glm = glm)
    
    
    invisible(gc())
    
    #store the stuff below:
    # the folds,
    # the full model
    # the selected shrinkage value
    # the indices selected in the bootstrap
    return(list(selected.model = models, selected.rows = this.row.idxs))
  }, mc.cores = num.cors) #end parallelization
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #the stuff above happens for each bootstrapped iteration
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  return(list(models = model.ens, dset = dat, response = y, outer.loop.folds = rep.row.idxs))
}

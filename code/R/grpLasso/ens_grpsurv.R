# Depends On: rep_grpsurv





ens_grpsurv <- function(data = mod_dat, preds = PARAMS$preds,  time = PARAMS$time,
                        stratified.cv = PARAMS$stratified.cv, response = PARAMS$response,
                        status = PARAMS$status, group = PARAMS$grp, penalty = PARAMS$penalty,
                        se = PARAMS$se, nlambda = PARAMS$nlambda, lsp = PARAMS$lsp,
                        num.boot = PARAMS$num.boot, outer.loop.cv.nfolds = PARAMS$outer.loop.cv.nfolds,
                        num.repeats = PARAMS$num.repeats, lambda.type = PARAMS$lambda.type,
                        n.folds = PARAMS$nfolds, alpha = PARAMS$alpha, returnY = PARAMS$returnY,
                        numCores = PARAMS$numCores, run.parallel = PARAMS$run.parallel, returnX = PARAMS$returnX) {
  
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
    rep.row.idxs <- foldsToIdxs(makeFolds(x = response, ntimes = num.boot, 
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
    if (!is.matrix(response) & !is.data.frame(response)) 
      response <- as.matrix(response)
    
    models <- rep_grpsurv(data = data[this.row.idxs, ], preds = preds,  time = time,
                          status = status, group = group, penalty = penalty,
                          se = se, nlambda = nlambda, lsp = lsp,
                          num.repeats = num.repeats, lambda.type = lambda.type,
                          n.folds = n.folds, alpha = alpha, returnY = returnY,
                          numCores = numCores, returnX = returnX)
      
    invisible(gc())
    
    #store the stuff below:
    # the folds,
    # the full model
    # the selected shrinkage value
    # the indices selected in the bootstrap
    return(list(enet.model = models$enet.model, mod.fam = models$mod.fam, err = models$err,
                selected.rows = this.row.idxs, selected.mod = models$selected.mod))
  }, mc.cores = num.cors) #end parallelization
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #the stuff above happens for each bootstrapped iteration
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  return(list(models = model.ens, data = data, preds = preds,  time = time, outer.loop.folds = rep.row.idxs))
}


#helper function to make the folds
makeFolds <- function(x, ntimes, nfolds, stratified.cv) {
  if (is.Surv(x)) { #means we have a survival model
    x <- as.matrix(x)[, "status"]
  } else {
    x <- as.matrix(x)
  }
  
  cv.idxs = generateCVRuns(x, ntimes = ntimes, nfold = nfolds, leaveOneOut = FALSE, 
                           stratified = stratified.cv)
  
  #  cv.idxs = vector(mode = "list", length = length(cv.idxs))
  for (i in seq_len(length(cv.idxs))) 
    cv.idxs[[i]] = foldsListToVector(cv.idxs[[i]], length(x))
  
  return(cv.idxs)   
}

#helper function to convert the output of makeFolds to a list of row indeces
foldsToIdxs <- function(folds.list) {
  num.folds <- length(unique(folds.list[[1]]))
  
  if(!(all(sapply(folds.list, function(fl) length(unique(fl))) == num.folds)))
    stop("all runs should have same number of folds (", num.folds, ")\nWe really shouldnt get here!")
  
  
  
  idxs.list <- vector(mode = "list", length = num.folds*length(folds.list))
  names(idxs.list) <- paste0(rep(gsub("  ", "_", names(folds.list)), each = num.folds), 
                             "_Fold_", 
                             rep(seq_len(num.folds), times = length(folds.list)))
  
  for (r.idx in seq_len(length(folds.list))) {
    for (f.idx in sort(unique(folds.list[[r.idx]]))) {
      tn              <- paste0("Run_", r.idx, "_Fold_", f.idx)
      idxs.list[[tn]] <- which(folds.list[[r.idx]] != f.idx)
    }
  }
  
  return(idxs.list)
}


#helper function to turn a list of folds into vector
foldsListToVector <- function(fl, len) {
  folds <- rep(0, length = len)
  for (k in seq_len(length(fl)))
    folds[fl[[k]]] <- k
  
  if (any(folds == 0))
    stop("Splitting samples into folds failed, likely because you are running stratified cross-validation and have too many unique values in your response variable")
  return(folds)
}

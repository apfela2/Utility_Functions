ens_discGlinternet <- function(X, Y, numLevels, n.folds = 10, num.repeats = 30, fold.ids = NULL, lambda.type = "min", lsp = NULL, nlambda = 50, 
                               lambdaMinRatio = 0.01, interactionCandidates = NULL, screenLimit = NULL, tol = 1e-05, maxIter = 5000, 
                               verbose = FALSE, numCores = 1, cens, intervals = min(X[[ Y ]]) : (max(X[[ Y ]])) + 1,
                               cat_X, col_contX, time_effect = T, run.parallel = FALSE, stratified.cv = FALSE, num.boot = 10,
                               sample.perc = .8, outer.loop.cv.nfolds = NA, group.multivar = FALSE) {
  
  
  num.cors = 1
  #setting up parallelization
  if (is.logical(run.parallel) & run.parallel == TRUE) {
    num.cors = max(1, (detectCores() - 1)) #registerDoMC(cores = min(n.folds, detectCores() -3))
    numCores = 1
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
    # if (!is.matrix(response) & !is.data.frame(response)) 
    #   response <- as.matrix(response)
    
    models <- rep_disc_Glinternet(X = X[this.row.idxs, ], Y = Y, numLevels = numLevels,
                                  n.folds = n.folds, num.repeats = num.repeats, fold.ids = NULL, lambda.type = lambda.type,
                                  lsp = lsp, nlambda = nlambda, lambdaMinRatio = lambdaMinRatio, interactionCandidates = interactionCandidates,
                                  screenLimit = NULL, tol = 1e-05, maxIter = maxIter, verbose = verbose, numCores = 1,
                                  cens = cens, intervals = intervals, cat_X = cat_X, col_contX = col_contX, time_effect = time_effect)
    
    
    invisible(gc())
    
    #store the stuff below:
    # the folds,
    # the full model
    # the selected shrinkage value
    # the indices selected in the bootstrap
    
    return(list(enet.model = models$enet.model, err = models$err,
                selected.rows = this.row.idxs, selected.mod = models$selected.mod))
  }, mc.cores = num.cors, mc.preschedule = FALSE) #end parallelization
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #the stuff above happens for each bootstrapped iteration
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  return(list(models = model.ens, X = X, Y = Y, outer.loop.folds = rep.row.idxs, cens = cens,
              cat_X = cat_X, col_cont_X = col_contX, intervals = intervals))
}

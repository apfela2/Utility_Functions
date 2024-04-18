# X is continuous survival data set
# Y is survival outcome
rep_disc_Glinternet <- function(X, Y, numLevels, n.folds = 10, num.repeats = 30, fold.ids = NULL, lambda.type = "min", lsp = NULL, nlambda = 50, 
                                lambdaMinRatio = 0.01, interactionCandidates = NULL, screenLimit = NULL, tol = 1e-05, maxIter = 5000, 
                                verbose = FALSE, numCores = 1, cens, intervals = min(X[[ Y ]]) : (max(X[[ Y ]])) + 1, cat_X, col_contX, time_effect = T) {
  
  
  # (dset.mod, response, mod.fam = "gaussian", lsp = NULL, nlambda = 50, n.folds = 10, 
  #                                num.repeats = 30, lambda.type = "min", numLevels, interactionCandidates = NULL, lambdaMinRatio = 0.01,
  #                                numCores = 1, stratified.cv = FALSE, fold.ids = NULL, verbose = FALSE,
  #                                diag.plots = FALSE) {
  
  #make sure inputted data and selected model family are consistent
  # checked.inputs <- checkPredsResp(x = dset.mod, y = response, 
  #                                  mod.fam = mod.fam)
  # 
  # dset          <- checked.inputs$dset
  # response.var  <- checked.inputs$response
  # mod.fam       <- checked.inputs$mod.fam
  # 
  # rm(checked.inputs)
  # invisible(gc())
  
  # if (mod.fam %in% c("binomial", "multinomial", "cox"))#"cox" | mod.fam == "binomial")
  #   stratified.cv = TRUE
  # 
  # rep.folds <- fold.ids
  # if (is.null(rep.folds) || is.na(rep.folds)) 
  #   rep.folds <- makeFolds(response, ntimes = num.repeats, nfolds = n.folds, stratified.cv = stratified.cv)
  
  
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #full glinternet model goes here
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  # Run full glinternet model to get lsp to be used across all cv samples
  if (is.null(lsp)){
    
    
    # Convert to discrete 
    discrete = disc(dat = X, cens, surv.time = Y, intervals = intervals)
    
    disc_genes = discrete[,col_contX]
    if (time_effect) {
      discrete$timeInt <- as.numeric(as.character(discrete$timeInt))
      disc_genes <- cbind(disc_genes, discrete$timeInt)
      numLevels <- c(numLevels, 1)
    }
    Y_disc = discrete$y
    X_disc = cbind(discrete[cat_X], disc_genes)
    
    full.glinternet <- glinternet(X_disc, Y_disc, numLevels = numLevels, lambda = lsp, nLambda = nlambda, 
                                  lambdaMinRatio = lambdaMinRatio, interactionCandidates = 1, screenLimit, family = "binomial", 
                                  tol = tol, maxIter = maxIter, verbose = verbose, numCores = numCores)
    
    # glinternet(X = dset.mod, Y = response, family = mod.fam, numLevels = numLevels,
    #                             lambdaMinRatio = lambdaMinRatio, lambda = lsp, nLambda = nlambda,
    #                             interactionCandidates = interactionCandidates, verbose = verbose, numCores = numCores)
    
    rm(discrete, disc_genes, X_disc, Y_disc)
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    
    #now get the lambda search path from the full model
    
    lsp <- full.glinternet$lambda
    numLevels <- numLevels[-length(numLevels)]
  }
  
  counter <- 1
  whileLimit <- 1
  models <- list()
  
  while (counter <= num.repeats & whileLimit <= num.repeats*2) {
    cv.res <- try({
      
      #run cv.glinternet
      discrete_glinternet_cv(X = X, Y = Y, numLevels = numLevels, nFolds = n.folds, lambda = lsp, nLambda = nlambda, 
                             lambdaMinRatio = lambdaMinRatio, interactionCandidates = interactionCandidates, screenLimit = NULL, 
                             family = "binomial", tol = tol, maxIter = maxIter, verbose = verbose, numCores = numCores,
                             cens = cens, intervals = intervals, cat_X = cat_X, col_contX = col_contX, time_effect = time_effect)
      
      
    }, silent=TRUE)
    
    if (class(cv.res) != "try-error") {
      
      #get all of the lambdas and the errors
      all.lambdas        <- cv.res$cvErr
      names(all.lambdas) <- cv.res$lambda
      
      
      models[[counter]] <- list(cv.all = all.lambdas)
      print(counter)
      counter <- counter + 1
    }
    
    whileLimit <- whileLimit + 1
    
    if(whileLimit != counter){
      cat(paste0("Counter = ", counter, "\n Total = ", whileLimit, "\n"))
    }
  }
  
  
  #now combine lambdas across all repeats
  err        <- do.call("rbind", lapply(models, function(this) this$cv.all))
  err.median <- apply(err, 2, median)
  err.sd     <- apply(err, 2, sd)
  
  #the valus of lambda are the names or err.median
  lambdas <- as.numeric(names(err.median))
  
  #using glmnet function to get the median error and error at one sd
  lvs <- glmnet::getmin(lambdas, err.median, err.sd)
  
  lambda.min <- lvs$lambda.min
  lambda.1se <- lvs$lambda.1se
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #finding the elbow
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  cvm.dir <- rev(diff(rev(err.median)))
  
  #only considering values of the second derivative between lambda.min and lambda1.se
  cvm.dir[lambdas < lambda.min | lambdas > lambda.1se] <- -Inf
  lambda.steep <- lambdas[which.max(cvm.dir)]
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  #storing lambda, error, and number of non zero coefs 
  #for each of "min", "1se", "elbow"
  cv.stats <- data.frame(min   = c(lambda.min, 
                                   as.numeric(err.median[paste0(lambda.min)])), 
                         # sum(coef.glmnet(full.glinternet, s = lambda.min) != 0)),
                         sparse = c(lambda.1se, 
                                    as.numeric(err.median[paste0(lambda.1se)])),
                         # sum(coef.glmnet(full.glinternet, s = lambda.1se) != 0)),
                         elbow  = c(lambda.steep, 
                                    as.numeric(err.median[paste0(lambda.steep)])))
  # sum(coef.glmnet(full.glinternet, s = lambda.steep) != 0)))
  # rownames(cv.stats) <- c("lambda", "meanCvError", "numNonZero")
  rownames(cv.stats) <- c("lambda", "meanCvError")
  
  
  #selecting the lambda of the specific lambda type
  this.mod <- cv.stats[, lambda.type, drop = FALSE]
  
  rm(models)
  invisible(gc())
  
  #store the stuff below:
  # the folds,
  # the full model
  # the selected shrinkage value
  # return(list(folds = rep.folds, enet.model = full.glinternet, selected.mod = this.mod,
  #             lambdaSearchSpace = lsp))
  return(list(enet.model = full.glinternet, selected.mod = this.mod, err = err,
              lambdaSearchSpace = lsp))
}

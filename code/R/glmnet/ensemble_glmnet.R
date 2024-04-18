# glmnetGridSearch
# ensembleGlmnet
# repeatedGlmnet
# nz.coef.glmnet
# sparse.mat.nz # helper function to get non-zero entities for sparse mat
# getModRefitErrors #calculate rmse, r-squared and adjusted r-squared
# getModRefitErrors_vec #calculate rmse, r-squared and adjusted r-squared
# makeFolds
# foldsToIdxs
# foldsListToVector
# pick.alpha.from.grid
# getEnsCoefs #function to get matrix of coeffecients from the refit models #NOTE: AG, May 8th, 2018 - this function is not yet finished!
# refitModel  #helper function to refit model based upon features selected by glmnet #intended for use inside glmnet functions


#this .R file contains functions (including helper functions)
#to fit a glmnet ensemble, and refit each model within the ensemble
#
#r-squared estimates are calculated as described in:
#http://stats.stackexchange.com/questions/11676/pseudo-r-squared-formula-for-glms

library(caret)
library(glmnet)
library(ggfortify)
library(ggplot2)
library(ComplexHeatmap)
library(pROC)
library(tidyr)
library(tibble)
library(ggfortify)
library(outliers)
library(xtable)
library(mclust)
library(MASS)
library(gridExtra)
library(dplyr)  # attach dplyr environment last so it's first in search path
library(TunePareto)
library(doMC)
library(survival)
library(visreg)
library(colorspace)
library(parallel)
library(foreach)


#first we provide a utlitity (based on Ron Ammar's code)
#for doing a grid search

#dset is the feature/design matrix ready for modelling
#response is a one column matrix with the response (outcome)
#n.folds is the number of folds
#n.grid is the number of point in both the alpha and lambda grid
glmnetGridSearch <- function(dset, response, mod.fam = "gaussian", n.folds = 3, stratified.cv = FALSE,
                             asp = seq(0, 1, length = 100), lsp = 2^seq(13, -13, length = 100), 
                             err.tol = .1, max.alpha = .95, num_cores = 1) {
  
  #if class is data.frame, suggest running makeModelReadyData
  if (class(dset) == "data.frame")
    stop("The input data to this function must be a model ready matrix. Try running makeModelReadyData first")
  
  alphaSearchSpace  <- asp #seq(0, 1, length = n)
  lambdaSearchSpace <- lsp #2^seq(13, -13, length = n)  # trying to use the bounds of glmnet
  NUM_FOLDS         <- n.folds
  
  #make sure inputted data and selected model faily are consistent
  checked.inputs <- checkPredsResp(x = dset, y = response, 
                                   mod.fam = mod.fam)
  
  dset          <- checked.inputs$dset
  response.var  <- checked.inputs$response
  mod.fam       <- checked.inputs$mod.fam
  response.name <- checked.inputs$response.name
  
  rm(checked.inputs)
  invisible(gc())
  
  if (mod.fam %in% c("binomial", "multinomial", "cox"))#"cox" | mod.fam == "binomial")
    stratified.cv = TRUE
  
  testFolds <- makeFolds(response.var, ntimes = 1, nfolds = n.folds, stratified.cv = stratified.cv)[[1]]
  
  #grid_search
  cvfit <- mclapply(alphaSearchSpace, 
                    function(a) cv.glmnet(x = dset, y = response.var,
                                          foldid = testFolds, family = mod.fam,
                                          type.measure = "deviance",
                                          lambda = lambdaSearchSpace, alpha = a),
                    mc.cores = num_cores)
  
  err     <- lapply(cvfit, "[[", "cvm")
  nZero   <- lapply(cvfit, "[[", "nzero")
  lambdas <- lapply(cvfit, "[[", "lambda")
  
  errDF <- data.frame()
  for(i in seq_len(length(err))) 
    errDF <- bind_rows(errDF, data.frame(alpha = rep(alphaSearchSpace[i], length(err[[i]])),
                                         lambda = lambdas[[i]], err = err[[i]],
                                         nzero = nZero[[i]]))
  
  
  best.alpha = pick.alpha.from.grid(errDF, tol = err.tol, max.alpha = max.alpha) #tol controls within what percent of 
  
  p <- ggplot(errDF, aes(x = alpha, y = log(lambda, base = 2), z = err, label = nzero)) + 
    geom_raster(aes(fill = err)) +
    scale_fill_distiller(palette = "Spectral") +
    # Use a formula to make the fewer non-zero coefs larger while making 0 tiny
    geom_text(aes(size = (1/(ifelse(nzero == 0, 100, nzero))))) +
    #geom_contour(aes(col=..level..)) +war
    labs(fill = "Error",   # use newline so not squeezing plot
         title = "glmnet grid search performance annotated with # of non-zero coefficients") +
    guides(size = FALSE)
  
  return(list(best.alpha = best.alpha, plot = p))
}

#response-var: name of response variable
#mod.fam: model family, assumes "gaussain" if not specified
#num.repeats: number of repeated cross-validation runs (i.e. num bootstraps)
#not that any character variables will be converted to factors

#dset is the data.frame of all predictor features, and non of the response features
#response.var is the response variable to be modeled
#sample.perc - percentage of samples used at each repeat, if NA then bootstrapping is done
#              a value of NA will result in bootstrapping
# outer.loop.cv.nfolds - number of folds for outer loop of cv. If left as NA repeated random
#                       samples will be done. Note that if not NA, the number of overall
#                       bootstraps will be num.boot * outer.loop.cv.nfolds
#group.multivar - whether or not to group multivariate response - matters only for multinomial (or multiple continuous responses,
#                  or multiple continuous responses, but this hasen't been used yet)
ensembleGlmnet <- function(dset, response, mod.fam = "gaussian", lsp = NULL, nlambda = 100, n.folds = 5, 
                           num.repeats = 30, glmn.alpha = .95, lambda.type = "sparse", 
                           run.parallel = FALSE, stratified.cv = FALSE, num.boot = 10, sample.perc = NA,
                           outer.loop.cv.nfolds = NA, group.multivar = FALSE, 
                           use.weights = FALSE, use.weights.npred = 1, back.trans = NULL, prop_lambda = 0) {
  
  num.cors = 1
  #setting up parallelization
  if (is.logical(run.parallel) & run.parallel == TRUE) {
    num.cors = max(1, parallel::detectCores() - 1) #registerDoMC(cores = min(n.folds, detectCores() -3))
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
    this.row.idxs <- NA
    
    if (length(rep.row.idxs) == 1 && is.na(rep.row.idxs)) {  
      #we are in scenario where outer loop of cv is not being done
      if (is.na(as.numeric(sample.perc))) {
        #if resample precentage isn't specified, do bootstrapping
        this.row.idxs <- sample(seq_len(nrow(dset)), replace = TRUE)     
      } else {
        #if resample percentage is specified, do repeated random samples
        this.row.idxs <- sample(seq_len(nrow(dset)), size = round(nrow(dset)*sample.perc), replace = FALSE) 
      }
    } else {
      this.row.idxs <- rep.row.idxs[[nb]]
    }
    
    #
    if (num.boot == 1) #if not doing any bootstrapping, use all samples
      this.row.idxs <- seq_len(nrow(dset))
    
    #we are running repeated glmnet for each bootstrap. setting run.parallel to FALSE because
    #we have already parallelized in the outside layer
    if (!is.matrix(response) & !is.data.frame(response)) 
      response <- as.matrix(response)
    
    #getting weights using residual of response regressed onto the use.weight.npred most correlated genes
    this.weights <- NULL
    if (use.weights) {
      #get the weights and residuals
      #w <- getWeightsFromResid(x = dset[this.row.idxs, ], 
      #                       y = response[this.row.idxs, , drop = FALSE], 
      #                       num.preds = use.weights.npred) 
      #dset         <- dset[, -which(colnames(dset) %in% w$used.genes)] #remove the gene used from the dataset
      #this.weights <- w$residuals
      this.weights <- as.numeric(as.matrix(response[this.row.idxs, , drop = FALSE]))
      
    }
    
    models <- repeatedGlmnet(dset[this.row.idxs, ], response = response[this.row.idxs, , drop = FALSE], 
                             mod.fam = mod.fam,
                             n.folds = n.folds, num.repeats = num.repeats, lsp = lsp, nlambda = nlambda,
                             glmn.alpha = glmn.alpha, lambda.type = lambda.type, 
                             run.parallel = ifelse(num.boot == 1, TRUE, FALSE), 
                             stratified.cv = stratified.cv, fold.ids = NULL, diag.plots = FALSE, 
                             group.multivar = group.multivar, weights = this.weights, prop_lambda = prop_lambda)
    
    if (!is.null(back.trans))
      models$back.trans <- back.trans
    
    invisible(gc())
    
    #store the stuff below:
    # the folds,
    # the full model
    # the selected shrinkage value
    # the indices selected in the bootstrap
    return(list(folds = models$folds, enet.model = models$enet.model, mod.fam = models$mod.fam, 
                selected.rows = this.row.idxs, selected.mod = models$selected.mod))
  }, mc.cores = num.cors) #end parallelization
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #the stuff above happens for each bootstrapped iteration
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  return(list(models = model.ens, dset = dset, response = response, outer.loop.folds = rep.row.idxs, back.trans = back.trans))
}


#internal function for running as part of ensemble glmnet, buyer beware
#setting for respeated glmnet, can set both the lambda search space and number of lambdas.
#best to leave lsp as null, and let glmnet decide. lsp will take precedence over nlambda. Using default
#value of nlambda, which is 100
repeatedGlmnet <- function(dset.mod, response, mod.fam = "gaussian", lsp = NULL, nlambda = 100, n.folds = 5, 
                           num.repeats = 30, glmn.alpha = .95, lambda.type = "min",
                           run.parallel = FALSE, stratified.cv = FALSE, fold.ids = NULL,
                           diag.plots = FALSE, group.multivar = TRUE, weights = NULL, prop_lambda = 0) {
  
  #make sure inputted data and selected model family are consistent
  checked.inputs <- checkPredsResp(x = dset.mod, y = response, 
                                   mod.fam = mod.fam)
  
  dset          <- checked.inputs$dset
  response.var  <- checked.inputs$response
  mod.fam       <- checked.inputs$mod.fam
  
  rm(checked.inputs)
  invisible(gc())
  
  err <- "deviance"
  if (mod.fam %in% c("binomial", "multinomial", "cox"))#"cox" | mod.fam == "binomial")
    stratified.cv = TRUE
  
  rep.folds <- fold.ids
  if (is.null(rep.folds) || is.na(rep.folds)) 
    rep.folds <- makeFolds(response, ntimes = num.repeats, nfolds = n.folds, stratified.cv = stratified.cv)
  
  if (is.null(weights))
    weights <- rep(1, nrow(dset.mod))
  
  #running the full model, inputting both lambda search space (lsp) and 
  #number of lambdas. Note that lsp takes precedence over nlambda
  full.glmn <- glmnet(x = dset, y = response.var, family = mod.fam, 
                      alpha = glmn.alpha, lambda = lsp, nlambda = nlambda, 
                      type.multinomial = ifelse(group.multivar, "grouped", "ungrouped"),
                      weights = weights)
  
  if (is.null(lsp))
    lsp <- full.glmn$lambda
  
  #for plotting:
  f.size <- 16 #setting font size
  n.non.zero.coefs <- sapply(lsp, function(s) { return(length(nz.coef.glmnet(full.glmn, s = s)))
  }) #number of non-zero betas at each lambda, adding one
  #for the intercept
  
  #registering parallel backend if we are running in parallel
  if (run.parallel) 
    doMC::registerDoMC(cores = min(n.folds, parallel::detectCores() -3))
  
  models <- lapply(seq_len(num.repeats), function(mi) {
    p <- NULL
    f <- rep.folds[[mi]]
    
    cv.res <- cv.glmnet(x = dset, y = response.var, foldid = f, lambda = lsp, 
                        alpha = glmn.alpha, family = mod.fam, 
                        parallel = run.parallel, type.measure = err, 
                        type.multinomial = ifelse(group.multivar, "grouped", "ungrouped"),
                        weights = weights)
    
    #getting the second dirivative of the mean errors
    #cvm.dir      <- secondDerivative(cv.res$cvm)
    #decided to use first derivative instead, more appropriate
    #using rev because differneces are computed in reverse order
    cvm.dir <- rev(diff(rev(cv.res$cvm)))
    
    #only considering values of the second derivative between lambda.min and lambda1.se
    cvm.dir[cv.res$lambda < cv.res$lambda.min | cv.res$lambda > cv.res$lambda.1se] <- -Inf
    
    lambda.steep <- cv.res$lambda[which.max(cvm.dir)]
    
    cv.stats <- data.frame(min   = c(cv.res$lambda.min, 
                                     cv.res$cvm[which(cv.res$lambda == cv.res$lambda.min)], 
                                     cv.res$nz[which(cv.res$lambda == cv.res$lambda.min)]),
                           sparse = c(cv.res$lambda.1se, 
                                      cv.res$cvm[which(cv.res$lambda == cv.res$lambda.1se)],
                                      cv.res$nz[which(cv.res$lambda == cv.res$lambda.1se)]),
                           elbow  = c(lambda.steep, 
                                      cv.res$cvm[cv.res$lambda == lambda.steep],
                                      cv.res$nz[cv.res$lambda == lambda.steep]))
    rownames(cv.stats) <- c("lambda", "meanCvError", "numNonZero")
    
    all.lambdas        <- cv.res$cvm
    names(all.lambdas) <- cv.res$lambda
    
    #make the plots here, since we have the glmnet.cv object at this point
    #writing my own plotting code because ggortify plot behaves weird when adding 
    #more layers
    if (diag.plots) {
      p[[1]] <- cvCurvePlot(cv.stats, lambdas = cv.res$lambda, nzeros = cv.res$nzero, 
                            err = cv.res$cvm, err.sd = cv.res$cvsd, 
                            err.name = cv.res$name, plot.title = paste0("CV Curve for Repeat ", mi))[[1]]
    }
    
    return(list(cv.stats = cv.stats, cv.all = all.lambdas, err.name = cv.res$name, plots = p))
  })
  
  
  #now combine lambdas across all repeats
  
  #sometimes cv.glmnet drops certain lambdas even if they are in the lambda search path
  #remove lambdas that not kept across at least minimal proportion of repeats
  lambda.count <- sapply(models, function(this) length(this$cv.all))
  
  lambda_matrix        <- do.call("bind_rows", lapply(models, 
                                                      function(this) this$cv.all))
  
  
  prop_miss <- apply(lambda_matrix, 2, function(x) sum(is.infinite(x))/nrow(lambda_matrix))
  
  err <- lambda_matrix[,which(prop_miss <= prop_lambda)]
  err.median <- apply(err, 2, median)
  err.sd     <- apply(err, 2, sd)
  n.non.zero.coefs <- n.non.zero.coefs[which(prop_miss <= prop_lambda)]
  
  #the valus of lambda are the names or err.median
  lambdas <- as.numeric(names(err.median))
  
  #using glmnet function to get the median error and error at one sd
  getmin=function(lambda,cvm,cvsd){
  cvmin=min(cvm,na.rm=TRUE)
  idmin=cvm<=cvmin
  lambda.min=max(lambda[idmin],na.rm=TRUE)
  idmin=match(lambda.min,lambda)
  semin=(cvm+cvsd)[idmin]
  idmin=cvm<=semin
  lambda.1se=max(lambda[idmin],na.rm=TRUE)
  list(lambda.min=lambda.min,lambda.1se=lambda.1se)
}
                     
lvs <- getmin(lambdas, err.median, err.sd)
  
  lambda.min <- lvs$lambda.min
  lambda.1se <- lvs$lambda.1se
  
  
  #get the elbow
  cvm.dir <- rev(diff(rev(err.median)))
  
  #only considering values of the second derivative between lambda.min and lambda1.se
  cvm.dir[lambdas < lambda.min | lambdas > lambda.1se] <- -Inf
  
  lambda.steep <- lambdas[which.max(cvm.dir)]
  
  cv.stats <- data.frame(min   = c(lambda.min, 
                                   as.numeric(err.median[paste0(lambda.min)]), 
                                   length(nz.coef.glmnet(full.glmn, s = lambda.min))),
                         sparse = c(lambda.1se, 
                                    as.numeric(err.median[paste0(lambda.1se)]),
                                    length(nz.coef.glmnet(full.glmn, s = lambda.1se))),
                         elbow  = c(lambda.steep, 
                                    as.numeric(err.median[paste0(lambda.steep)]),
                                    length(nz.coef.glmnet(full.glmn, s = lambda.steep))))
  rownames(cv.stats) <- c("lambda", "meanCvError", "numNonZero")
  
  
  #get errors for each lambda type
  all.lambda.types <- colnames(cv.stats)
  
  errDist <- lapply(all.lambda.types, function(et) {
    return(err[, paste0(cv.stats["lambda", et])])
  })
  names(errDist) <- all.lambda.types
  
  if (!lambda.type %in% all.lambda.types)
    stop("Selected lambda.type must be one of (", paste(all.lambda.types, collapse = ", "), "), but you selected ", lambda.type)
  
  this.mod <- cv.stats
  
  #make diag plots
  p <- NULL
  
  if (diag.plots) {
    #adding summary plots across the repeated runs
    p <- append(p, repeatedCVSummaryPlot(full.glmn, models, f.size = 16, family = mod.fam))
    #consider making a 3D plot
    
    p <- append(p,  
                cvCurvePlot(cv.stats, lambdas, nzeros = n.non.zero.coefs, err = err.median, err.sd = err.sd, 
                            err.name = models[[1]]$err.name, plot.title = paste0("CV Curve for Aggregated Median Errors Across ", 
                                                                                 num.repeats, " Repeats")))
    
    #plots for each repeat
    p.idx = length(p) + 1
    #loop over models and make plots
    for (mi in seq_len(length(models))) {
      p[[p.idx]] <- models[[mi]]$plots[[1]]
      p.idx = p.idx + 1
    }
    
  }
  rm(models)
  invisible(gc())
  
  #store the stuff below:
  # the folds,
  # the full model
  # the selected shrinkage value
  return(list(folds = rep.folds, enet.model = full.glmn, selected.mod = this.mod, mod.fam = mod.fam, 
              lambdaSearchSpace = lsp, plots = p))
}

#helper function to get number of non-zero coefficients
#needed because if family is multionmoial the output of coef.glmnet is a list
nz.coef.glmnet <- function(mod, s, collapse.list = TRUE) {
  cc <- coef.glmnet(mod, s = s)
  if (class(cc) != "list")
    return(rownames(cc)[which(cc != 0)])
  
  #multinomial returns a list, catching that here
  
  if (collapse.list)
    return(unique(unlist(lapply(cc, function(this) rownames(this)[which(this != 0)]))))
  
  return(unique(unlist(sapply(cc, function(this) rownames(this)[which(this != 0)]))))
}

#helper function to get non-zero entities for sparse mat
sparse.mat.nz <- function(x) {
  nz <- x[which(x != 0)]
  names(nz) <- rownames(x)[which(x != 0)]
  return(nz)
}


#calculate rmse, r-squared and adjusted r-squared
getModRefitErrors <- function(mf) {
  rv <- as.character(mf$formula)[2]
  
  y     <- mf$data[, rv]
  y.hat <- predict(mf)
  n     <- length(y)
  p     <- length(coef(mf)) - 1
  
  rmse <- sqrt(sum((y - y.hat)^2)/n)
  
  r.sq     <- 1 - sum((y - y.hat)^2)/sum((y - mean(y))^2)
  adj.r.sq <- 1 - (1 - r.sq)*((n - 1)/(n - p - 1))
  
  return(list(rmse = rmse, r.squared = r.sq, adj.r.squared = adj.r.sq))
}

#calculate rmse, r-squared and adjusted r-squared
getModRefitErrors_vec <- function(y, y.hat, nz) {
  n     <- length(y)
  p     <- nz
  
  rmse <- sqrt(sum((y - y.hat)^2)/n)
  
  r.sq     <- 1 - sum((y - y.hat)^2)/sum((y - mean(y))^2)
  adj.r.sq <- 1 - (1 - r.sq)*((n - 1)/(n - p - 1))
  
  return(list(rmse = rmse, r.squared = r.sq, adj.r.squared = adj.r.sq))
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

#picks alpha that results in smallest error and 
#not the null model
#columns of err.df are: alpha, lambda, misclass (err rate), nzero
pick.alpha.from.grid <- function(err.df, tol = .1, max.alpha = .95) {
  #keep only rows that are within tolerance
  err.df = err.df %>% dplyr:::filter(err <= quantile(err, 1 - tol)) %>%
    dplyr:::filter(alpha <= max.alpha)
  
  #since all of the errors are within our tolerance, we sort
  #by alpha, and pick the largest alpha (sparset model) that 
  #still has some non-zero predictors
  err.df = err.df[order(err.df$alpha, decreasing = TRUE), ]
  
  if (all(err.df$nzero == 0))
    return(1) #this means the l1-only (lasso) model will be run
  
  return(err.df$alpha[which(err.df$nzero != 0)[1]]) #pick the alpha corresponding to the first
  #entry with any terms in the model
}



#function to get matrix of coeffecients from the refit models
#NOTE: AG, May 8th, 2018 - this function is not yet finished!
getEnsCoefs <- function(ens, refit = T, lambda.type = "min") {
  all.coefs <- lapply(ens$models, function(m) {
    tc           <- coef.glmnet(m$enet.model, s = m$selected.mod["lambda", lambda.type])
    tc.nz        <- suppressMessages(tc[tc != 0])
    names(tc.nz) <- rownames(tc)[which(tc != 0)]
    
    intc <- "(Intercept)"
    if (any(names(tc.nz) == intc))
      tc.nz <- tc.nz[-which(names(tc.nz) == intc)]
    
    return(tc.nz)
  })
  
  #cts       <- sapply(all.coefs, length)
  
  unq.coefs  <- unique(unlist(sapply(all.coefs, names)))
  ebetas.mat <- matrix(NA, length(ens$models), length(unq.coefs))
  rownames(ebetas.mat) <- paste0("mod_", seq_len(length(ens$models)))
  colnames(ebetas.mat) <- unq.coefs
  
  betas.mat <- ebetas.mat
  pv.mat    <- ebetas.mat
  
  if (refit == FALSE) {
    betas.mat <- NULL
    pv.mat    <- NULL
  }
  
  #loop over each model, call the index mi, stands for model index
  for (mi in seq_len(length(ens$models))) {
    #get the enet coefs, which we have already calculated
    ebetas.mat[mi, names(all.coefs[[mi]])] <- all.coefs[[mi]]
    
    #if we need to refit the model
    tr <- ens$response[ens$models[[mi]]$selected.rows, , drop = FALSE] 
    tp <- ens$dset[ens$models[[mi]]$selected.rows, names(all.coefs[[mi]]), drop = FALSE]
    
    #refit the model
    refit.mod <- refitModel(response = tr, predictors = tp, family = ens$models[[mi]]$mod.fam)
    
    #store the results
    betas.mat[mi, names(refit.mod$betas)] <- refit.mod$betas
    pv.mat[mi, names(refit.mod$pvs)] <- refit.mod$pvs
  }
  
  #consider adding code to rescale parameters
  #and option to trim based on number of features in model...maybe this happens downstream
  
  return(list(enet.betas  = data.frame(ebetas.mat, stringsAsFactors = FALSE),
              refit.betas = data.frame(betas.mat, stringsAsFactors = FALSE),
              refit.pvs   = data.frame(pv.mat, stringsAsFactors = FALSE)))
}

#helper function to refit model based upon features selected by glmnet
#intended for use inside glmnet functions
refitModel <- function(response, predictors, family) {
  if (family == "cox")
    stop("Have not implemented refitting for coxph!")
  
  
  this.dset <- data.frame(cbind(response, predictors), stringsAsFactors = FALSE)
  
  refit.model    <- glm(as.formula(paste0(names(response), " ~ .")), data = this.dset, family = family)
  refit.mod.coef <- summary(refit.model)$coefficients
  
  intc <- "(Intercept)"
  if (any(rownames(refit.mod.coef) == intc))
    refit.mod.coef <- refit.mod.coef[-which(rownames(refit.mod.coef) == intc), ]
  
  #now we get the betas and the pvalues
  #lines below may change based on model type
  
  #get the pvalue column
  pvc <- grep("Pr", colnames(refit.mod.coef))
  if (length(pvc) > 1)
    stop("Can't have more than one p-value column. This is a silly error")
  
  ret <- list(betas = refit.mod.coef[, "Estimate", drop = T],
              pvs   = refit.mod.coef[, pvc, drop = T])
  return(ret)
}

## Abraham Apfel
## July 18, 2018


# Create function to do repeated CV for glinternet based off of Alex's repeated glmnet function


# Needs the following scripts
# source("submodules/P00465_predictive_modeling_framework/ensemble_glmnet.R")
# source("submodules/P00465_predictive_modeling_framework/ensemble_glmnet_plotting.R")
# source("submodules/P00465_predictive_modeling_framework/data_cleaning.R")

#template function for repeatedGlinternet
repeatedGlinternet <- function(dset.mod, response, mod.fam = "gaussian", lsp = NULL, nlambda = 50, n.folds = 10, 
                               num.repeats = 30, lambda.type = "min", numLevels, interactionCandidates = NULL, lambdaMinRatio = 0.01,
                               numCores = 1, stratified.cv = FALSE, fold.ids = NULL, verbose = FALSE,
                               diag.plots = FALSE) {
  
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
  
  if (mod.fam %in% c("binomial", "multinomial", "cox"))#"cox" | mod.fam == "binomial")
    stratified.cv = TRUE
  
  # rep.folds <- fold.ids
  # if (is.null(rep.folds) || is.na(rep.folds)) 
  #   rep.folds <- makeFolds(response, ntimes = num.repeats, nfolds = n.folds, stratified.cv = stratified.cv)
  
  
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #full glinternet model goes here
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  # numLevels <- rep(NA, ncol(dset.mod))
  # for (i in 1:ncol(dset.mod)) {
  #   if (class(dset.mod[,i]) == "factor") {
  #     levels(dset.mod[,i]) <- 1:length(levels(dset.mod[,i]))
  #     numLevels[i] = length(levels(dset.mod[,i]))
  #   } else if (class(dset.mod[,i]) == "numeric") {
  #     numLevels[i] = 1
  #   } else {
  #     stop("Error.  Columns must be numeric or factor")
  #   }
  # }
  
  
  full.glinternet <- glinternet(X = dset.mod, Y = response, family = mod.fam, numLevels = numLevels,
                                lambdaMinRatio = lambdaMinRatio, lambda = lsp, nLambda = nlambda,
                                interactionCandidates = interactionCandidates, verbose = verbose, numCores = numCores)
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  #now get the lambda search path from the full model
  if (is.null(lsp))
    lsp <- full.glinternet$lambda
  
  
  counter <- 1
  whileLimit <- 1
  models <- list()
  
  while (counter <= num.repeats & whileLimit <= num.repeats*2) {
    cv.res <- try({
      
      #run cv.glinternet
      glinternet.cv(X = dset.mod, Y = response, family = mod.fam, numLevels = numLevels,
                    lambdaMinRatio = lambdaMinRatio, lambda = lsp, nLambda = nlambda, nFolds = n.folds,
                    interactionCandidates = interactionCandidates, verbose = verbose, numCores = numCores)
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
  return(list(enet.model = full.glinternet, selected.mod = this.mod, err = err, mod.fam = mod.fam,
              lambdaSearchSpace = lsp))
}

############################################################################################
# survAssess <- function(pred.mat, timepoint, marker, newTime, newEvent, trainTime, trainEvent)
#   
#   
#   full_temp <- contToDisc(pred.mat$dset, timeColumn = pred.mat$time,  intervalLimits = pred.mat$intervals)
# 
# tpr <- list()
# fpr <- list()
# ROCcurves <- list()
# 
# for (repeats in 1:ncol(pred.mat$pred.mat)){
#   for (fold in 1: (length(pred.mat$idxs)/ncol(pred.mat$pred.mat))){
#     full_train <- full_temp[pred.mat$idxs[[ repeats*fold]],]
#     full_test <- full_temp[-pred.mat$idxs[[ repeats*fold]],]
#     
#     quant <- c(5, 13, 21)
#     
#     tpr[[(5*(repeats-1)+fold)]] <- list()
#     fpr[[(5*(repeats-1)+fold)]] <- list()
#     ROCcurves[[(5*(repeats-1)+fold)]] <- list()
#     
#     for (time in 1:length(quant)) {
#       
#       tpr[[(5*(repeats-1)+fold)]][[time]] <- tprUnoShort(timepoint = quant[[time]], marker = pred.mat$pred.mat[-pred.mat$idxs[[ repeats*fold]], repeats], newTime = full_test$timeDisc, newEvent = full_test$PFSINV.EVENT,
#                                                trainTime = full_train$timeDisc, trainEvent = full_train$PFSINV.EVENT)
#       fpr[[(5*(repeats-1)+fold)]][[time]] <- fprUnoShort(timepoint = quant[[time]], marker = pred.mat$pred.mat[-pred.mat$idxs[[ repeats*fold]], repeats], newTime = full_test$timeDisc)
#       ROCcurves[[(5*(repeats-1)+fold)]][[time]] <- aucUno(tpr[[(5*(repeats-1)+fold)]][[time]], fpr[[(5*(repeats-1)+fold)]][[time]])
#     }
#   }
#   # quant_temp <- quantile(as.numeric(as.character(full_temp$timeDisc)), probs = seq(.2, .8, .2))
#   # 
#   # for (time in 1:length(quant_temp)) {
#   #   fpr[[repeats]][time] <- fprUnoShort(timepoint = quant_temp[[time]], marker = pred.mat$pred.mat[, repeats], newTime = full_temp$timeDisc)
#   # }
# }
#

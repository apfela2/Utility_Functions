## Abraham Apfel
## September 20, 2019


# Function to run inner loop of repeated cvgrpsurv
# data <- dat_scale
# preds <- c(genes, INT_genes, mod_DISETI, mod_trt, mod_clin)
# time <- "OSDAYS"
# status <- "OSEVENT"
# penalty = "grLasso"
# lsp <- NULL
# nlambda <- 100
# n.folds = 5
# num.repeats = 30
# alpha = .95
# # gamma = 0.8
# lambda.type = "min"
# # standardize = F
# run.parallel = FALSE
# stratified.cv = FALSE
# fold.ids = NULL
# diag.plots = FALSE
# weights = NULL
# prop_lambda = 0
# numCores = 30
# returnX = T
# returnY = T
# se = "bootstrap"
# 
# 
# group <- c(rep(1:length(genes), times = 2), rep(length(genes) + 1, times = length(mod_DISETI)),
#          seq(length(genes) + 2, length(genes) + 1 + length(mod_trt) + length(mod_clin)))
# 
# 

rep_grpsurv <- function(data, preds, time, status, group, penalty = "grLasso", lsp = NULL, nlambda = 50,
                        n.folds = 5, fold.ids, se = 'quick', num.repeats = 30, alpha = .95,
                        lambda.type = "min", numCores = 1, tau=1/3,
                        group.multiplier, warn=TRUE, returnX=FALSE, returnY=FALSE, trace=FALSE) {
  
  # #make sure inputted data and selected model family are consistent
  # checked.inputs <- checkPredsResp(x = dset.mod, y = response, 
  #                                  mod.fam = mod.fam)
  # 
  # dset          <- checked.inputs$dset
  # response.var  <- checked.inputs$response
  # mod.fam       <- checked.inputs$mod.fam
  # 
  # rm(checked.inputs)
  # invisible(gc())
  
  # Configure data for SGL function
  
  # X <- dat[, X]
  # X <- as.matrix(X)
  
  # Check to make sure there is no missing
  if(any(is.na(data[, c(preds, time, status)])))
    stop("Remove missing data before running model")
  
  X <- as.matrix(data[, preds])
  y <- as.matrix(data.frame(data[,time],
                            data[,status]))
  
  
  # err <- "Negative Log Likelihood"
  # if (type %in% c("binomial", "multinomial", "cox"))#"cox" | mod.fam == "binomial")
  #   stratified.cv = TRUE
  # 
  # rep.folds <- fold.ids
  # if (is.null(rep.folds) || is.na(rep.folds)) 
  #   rep.folds <- makeFolds(mod_dat$status, ntimes = num.repeats, nfolds = n.folds, stratified.cv = stratified.cv)
  # 
  
  #running the full model, inputting both lambda search space (lsp) and 
  #number of lambdas. Note that lsp takes precedence over nlambda
  full.grpsurv <- grpsurv(X = X, y = y, group = group, penalty = penalty,
                          alpha = alpha, nlambda = nlambda, returnX = returnX)
  
  
  if (is.null(lsp))
    lsp <- full.grpsurv$lambda
  
  #for plotting:
  f.size <- 16 #setting font size
  n.non.zero.coefs <- NULL
  for(s in 1:ncol(full.grpsurv$beta)) {
    n.non.zero.coefs[s] <- length(full.grpsurv$beta[,s][which(full.grpsurv$beta[,s] != 0)])
  } #number of non-zero betas at each lambda, adding one
  #for the intercept
  
  
  #registering parallel backend if we are running in parallel
  # if (run.parallel) 
  #   doMC::registerDoMC(cores = min(n.folds, parallel::detectCores() -3))
  
  #   models <- lapply(seq_len(num.repeats), function(mi) {
  #     p <- NULL
  #     f <- rep.folds[[mi]]
  #     
  #     cv.res <- cvSGL(data = mod_dat, index = index, lambdas = lsp, nfold = n.folds,
  #                         alpha = alpha, gamma = gamma, step = step, reset = reset,
  #                     type = type, nlam = nlambda, standardize = standardize)
  #     
  #     #getting the second dirivative of the mean errors
  #     #cvm.dir      <- secondDerivative(cv.res$cvm)
  #     #decided to use first derivative instead, more appropriate
  #     #using rev because differneces are computed in reverse order
  #     cvm.dir <- rev(diff(rev(cv.res$cvm)))
  #     
  #     #only considering values of the second derivative between lambda.min and lambda1.se
  #     cvm.dir[cv.res$lambda < cv.res$lambda.min | cv.res$lambda > cv.res$lambda.1se] <- -Inf
  #     
  #     lambda.steep <- cv.res$lambda[which.max(cvm.dir)]
  #     
  #     cv.stats <- data.frame(min   = c(cv.res$lambda.min, 
  #                                      cv.res$cvm[which(cv.res$lambda == cv.res$lambda.min)], 
  #                                      cv.res$nz[which(cv.res$lambda == cv.res$lambda.min)]),
  #                            sparse = c(cv.res$lambda.1se, 
  #                                       cv.res$cvm[which(cv.res$lambda == cv.res$lambda.1se)],
  #                                       cv.res$nz[which(cv.res$lambda == cv.res$lambda.1se)]),
  #                            elbow  = c(lambda.steep, 
  #                                       cv.res$cvm[cv.res$lambda == lambda.steep],
  #                                       cv.res$nz[cv.res$lambda == lambda.steep]))
  #     rownames(cv.stats) <- c("lambda", "meanCvError", "numNonZero")
  #     
  #     all.lambdas        <- cv.res$cvm
  #     names(all.lambdas) <- cv.res$lambda
  #     
  #     #make the plots here, since we have the glmnet.cv object at this point
  #     #writing my own plotting code because ggortify plot behaves weird when adding 
  #     #more layers
  #     if (diag.plots) {
  #       p[[1]] <- cvCurvePlot(cv.stats, lambdas = cv.res$lambda, nzeros = cv.res$nzero, 
  #                             err = cv.res$cvm, err.sd = cv.res$cvsd, 
  #                             err.name = cv.res$name, plot.title = paste0("CV Curve for Repeat ", mi))[[1]]
  #     }
  #     
  #     return(list(cv.stats = cv.stats, cv.all = all.lambdas, err.name = cv.res$name, plots = p))
  #   })
  #   
  #   
  #   #now combine lambdas across all repeats
  #   
  #   #sometimes cv.glmnet drops certain lambdas even if they are in the lambda search path
  #   #remove lambdas that not kept across at least minimal proportion of repeats
  #   lambda.count <- sapply(models, function(this) length(this$cv.all))
  #   
  #   lambda_matrix        <- do.call("bind_rows", lapply(models, 
  #                                                       function(this) this$cv.all))
  #   
  #   
  #   prop_miss <- apply(lambda_matrix, 2, function(x) sum(is.infinite(x))/nrow(lambda_matrix))
  #   
  #   err <- lambda_matrix[,which(prop_miss <= prop_lambda)]
  #   err.median <- apply(err, 2, median)
  #   err.sd     <- apply(err, 2, sd)
  #   n.non.zero.coefs <- n.non.zero.coefs[which(prop_miss <= prop_lambda)]
  #   
  #   #the valus of lambda are the names or err.median
  #   lambdas <- as.numeric(names(err.median))
  #   
  #   #using glmnet function to get the median error and error at one sd
  #   lvs <- glmnet::getmin(lambdas, err.median, err.sd)
  #   
  #   lambda.min <- lvs$lambda.min
  #   lambda.1se <- lvs$lambda.1se
  #   
  #   
  #   #get the elbow
  #   cvm.dir <- rev(diff(rev(err.median)))
  #   
  #   #only considering values of the second derivative between lambda.min and lambda1.se
  #   cvm.dir[lambdas < lambda.min | lambdas > lambda.1se] <- -Inf
  #   
  #   lambda.steep <- lambdas[which.max(cvm.dir)]
  #   
  #   cv.stats <- data.frame(min   = c(lambda.min, 
  #                                    as.numeric(err.median[paste0(lambda.min)]), 
  #                                    length(nz.coef.glmnet(full.glmn, s = lambda.min))),
  #                          sparse = c(lambda.1se, 
  #                                     as.numeric(err.median[paste0(lambda.1se)]),
  #                                     length(nz.coef.glmnet(full.glmn, s = lambda.1se))),
  #                          elbow  = c(lambda.steep, 
  #                                     as.numeric(err.median[paste0(lambda.steep)]),
  #                                     length(nz.coef.glmnet(full.glmn, s = lambda.steep))))
  #   rownames(cv.stats) <- c("lambda", "meanCvError", "numNonZero")
  #   
  #   
  #   #get errors for each lambda type
  #   all.lambda.types <- colnames(cv.stats)
  #   
  #   errDist <- lapply(all.lambda.types, function(et) {
  #     return(err[, paste0(cv.stats["lambda", et])])
  #   })
  #   names(errDist) <- all.lambda.types
  #   
  #   if (!lambda.type %in% all.lambda.types)
  #     stop("Selected lambda.type must be one of (", paste(all.lambda.types, collapse = ", "), "), but you selected ", lambda.type)
  #   
  #   this.mod <- cv.stats
  #   
  #   #make diag plots
  #   p <- NULL
  #   
  #   if (diag.plots) {
  #     #adding summary plots across the repeated runs
  #     p <- append(p, repeatedCVSummaryPlot(full.glmn, models, f.size = 16, family = mod.fam))
  #     #consider making a 3D plot
  #     
  #     p <- append(p,  
  #                 cvCurvePlot(cv.stats, lambdas, nzeros = n.non.zero.coefs, err = err.median, err.sd = err.sd, 
  #                             err.name = models[[1]]$err.name, plot.title = paste0("CV Curve for Aggregated Median Errors Across ", 
  #                                                                                  num.repeats, " Repeats")))
  #     
  #     #plots for each repeat
  #     p.idx = length(p) + 1
  #     #loop over models and make plots
  #     for (mi in seq_len(length(models))) {
  #       p[[p.idx]] <- models[[mi]]$plots[[1]]
  #       p.idx = p.idx + 1
  #     }
  #     
  #   }
  #   rm(models)
  #   invisible(gc())
  #   
  #   #store the stuff below:
  #   # the folds,
  #   # the full model
  #   # the selected shrinkage value
  #   return(list(folds = rep.folds, enet.model = full.glmn, selected.mod = this.mod, mod.fam = mod.fam, 
  #               lambdaSearchSpace = lsp, plots = p))
  # }
  
  
  models <- mclapply(seq_len(num.repeats), function(x) {
    
    
    # for (i in 1:num.repeats){
    #run cv.glinternet
    cv.res <- cv.grpsurv(X = X, y = y, group = group, penalty = penalty, nfolds = n.folds, returnY = returnY,
                         alpha = alpha, nlambda = nlambda, returnX = returnX, se = se)
    
    
    #get all of the lambdas and the errors
    all.lambdas        <- cv.res$cve
    names(all.lambdas) <- cv.res$lambda
    
    return(list(cv.all = all.lambdas))
  }, mc.cores = numCores)
  
  
  
  
  #now combine lambdas across all repeats
  err        <- do.call("rbind", lapply(models, function(this) this$cv.all))
  err.median <- apply(err, 2, median)
  err.sd     <- apply(err, 2, sd)
  
  #the valus of lambda are the names or err.median
  lambdas <- as.numeric(names(err.median))
  
  #using glmnet function to get the median error and error at one sd
  lambda.min <- lambdas[which.min(err.median)]
  lambda.1se <- lambdas[which.max(as.numeric(as.character(names(err.median[err.median <= (err.median[which.min(err.median)] + err.sd[which.min(err.median)])]))))]
  
  # lvs <- glmnet::getmin(lambdas, err.median, err.sd)
  # 
  # lambda.min <- lvs$lambda.min
  # lambda.1se <- lvs$lambda.1se
  
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
  return(list(enet.model = full.grpsurv, selected.mod = this.mod, err = err, mod.fam = "cox",
              lambdaSearchSpace = lsp))
}

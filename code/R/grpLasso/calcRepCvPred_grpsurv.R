calcRepCvPred_grpsurv <- function(ens, lambda_type = "min", run.parallel = TRUE, refit = FALSE) {
  num.cores <- ifelse(run.parallel, parallel::detectCores() - 1, 1)
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #   some sanity checks
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  
  #1) the object ens.cont$outer.loop.folds contains the indices used for training in each fold, making sure
  # these are the same as the selected.rows field in each model object
  diff.idxs <- sapply(seq_len(length(ens$models)), function(idx) sum(ens$models[[idx]]$selected.rows != ens$outer.loop.folds[[idx]]))
  if (sum(diff.idxs) != 0)
    stop("Some indices don't match up....something is very wrong")
  #again, all good so far, no obvious bugs
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  #we will loop over each repeated run of 5 fold cv, and for each one we calculate performance on
  #test set in two ways
  #   1) constrained coefficients from glmnet
  
  
  #for book keeping create a data frame which maps each model which cv repeat and fold it comes from
  mod.map <- data.frame(cv.rep = as.numeric(sapply(strsplit(names(ens$outer.loop.folds), "_"), function(s) s[[2]])),
                        fold = as.numeric(sapply(strsplit(names(ens$outer.loop.folds), "_"), function(s) s[[4]])))
  rownames(mod.map) <- names(ens$outer.loop.folds)
  
  preds <- mclapply(seq_len(length(ens$models)), function(idx){ 
    
    
    tm <- ens$models[[idx]] #get the model at this idx
    
    #get the indices used for training
    train.rows <- tm$selected.rows
    test.rows  <- seq_len(nrow(ens$data))[-train.rows]
    
    
    this.s <- tm$selected.mod["lambda", lambda_type] #get value of lambda based on selected lambda type
    lambda_ind <- which(round(tm$enet.model[["lambda"]], digits = 6) == round(tm$selected.mod["lambda", lambda_type], digits = 6))
    
    if (!refit) {
      this.pred <- predict(object = tm$enet.model, X = as.matrix(ens$data[test.rows, ens$preds]), type = "link", lambda = tm$enet.model[["lambda"]][lambda_ind])
      
      return(list(pred = this.pred, test.rows = test.rows))
    } else {
      #get selected features at chosen lambda
      this.feat <- names(getNZCoefs(coef.glmnet(tm$enet.model, s = this.s), remove_intercept = T))
      
      #get the  data
      this.md <- cbind(ens$response[train.rows, , drop = F], ens$data[train.rows, this.feat])
      
      #build lm model
      this.mod <- glm(as.formula(paste0(names(ens$response), " ~ .")), data = data.frame(this.md, stringsAsFactors = FALSE))
      
      #make prediction
      this.pred <- predict(this.mod, newdata = data.frame(ens$data[test.rows, ]))
      
      
      return(list(pred = this.pred, test.rows = test.rows))
    }
  }, mc.cores = num.cores)
  
  #setting up matrices to store results
  perf.mat <- matrix(NA, nrow(ens$data), max(mod.map$cv.rep))
  
  for (idx in seq_len(length(ens$models))) {
    #cat(".")
    cv.rep    <- mod.map$cv.rep[idx] #which cv.rep are we on
    test.rows <- preds[[idx]]$test.rows
    
    #make sure none oe the test rows for this repeat of cv have results,
    #if they do we've double counted somehwere and have an error
    if (any(!is.na(perf.mat[test.rows, cv.rep])))
      stop("We've already calculated performance for at least one of the test indices...there is a bug somewhere")
    
    
    perf.mat[test.rows, cv.rep] <- preds[[idx]]$pred
    
    rm(cv.rep, test.rows) 
  }
  
  coefs <- mclapply(ens$models, function(m) {
    lambda_ind <- which(round(m$enet.model[["lambda"]], digits = 6) == round(m$selected.mod["lambda", lambda_type], digits = 6))
    return(coef(m$enet.model, lambdaIndex = lambda_ind))
  }, mc.cores = num.cores)
  
  return(list(pred.mat = perf.mat, coefs = coefs, idxs = ens$outer.loop.folds))
}

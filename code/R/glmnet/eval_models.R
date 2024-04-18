# getNZCoefs
# getAllCoefs
# recal
# predict.repglmnet
# evalPred #helper function to calculate prediction error for continuous outcome
# compBinPerf #function to compute performance (ppv, npv, opv) between two binary variables
# evalBinPred #helper function to calculate prediction error for binary outcomes
# evalModel #evaluate an enet model on the full data hyperparameters fit using repeated-cv
# calcRepCvPred_enet #evaluate the full process, creating out of sample predictions for the repeated folds
# evalEnsModel #calculate the performance metrics given a matrix of out of sample predictions
# calcRepCvPred_lm #eval a pre-specifiied linear model with repeated-cv #dset - dataset with all of the features to be used in model building #response - response #outer.loop.folds - name list which contains training indices for each model, for each cv
# evalLmRepCv
# calcBootPred_enet #evaluate the full process, creating out of sample predictions for the repeated folds
# getWeightsFromResid #helper functions to get abs of residuals to use as weights
# makeThresholdTables #helper function to make table of performance metrics and threhsold where npv, ppv cross eachother, and where npv, pvv cross the selected cutoff values, first two cutoffs are npv, the second two are ppv agg.dat is an aggregated data frame of the results from repeated cv
# getComparableCoefs #helper function to get scaled coefs #inputs: 1) features - names of the selected coefficients to re-scale 2) dset - dataset of predictors (features) 3) response - the response variable 4) mod.fam - model family (e.g. gaussian, bionomial, etc.) 5) scale.type - one of "IQR" or "SD", which scales the (refit) betas by either IQR or SDs. #outputs: 1) scaled.betas - vector of selected features refit to full dataset and scaled as specified #note that these rescaled parameters are invalid for inference or further model testing. These betas can be used solely to compare the relative importance of features





#helper function to get non-zero coefficints from elasticnet
getNZCoefs <- function(x, remove_intercept = TRUE, intercept_string = "(Intercept)") {
  xx <- x[which(x[, 1] != 0), ]
  
  if (remove_intercept & any(names(xx) == intercept_string))
    xx <- xx[-which(names(xx) == intercept_string)]
  
  return(xx)
}

#helper function to get all coefficients
getAllCoefs <- function(xx, remove_intercept = TRUE, intercept_string = "(Intercept)" ) {
  if (remove_intercept & any(rownames(xx) == intercept_string))
    xx <- xx[-which(rownames(xx) == intercept_string), , drop = FALSE]
  
  return(rownames(xx))
}

#helper function to recalibrate a model
#this.y - resposne variable
#this.y.hat - predicted/fitted values of y
recal <- function(this.y, this.y_hat, method = "lm", test.dat = NULL) {
  if (method == "lm") {
    this.lm <- lm(this.y ~ this.y_hat)
    
    #return the fitted values if no test data supplied
    if (is.null(test.dat))
      return(this.lm$fitted.values)
    
    #make predictions on test data
    return(suppressWarnings(predict(this.lm, newdata = data.frame(test.dat))))
  }
  else if (method == "rlm") {
    this.rlm <- rlm(this.y ~ this.y_hat)
    
    #return the fitted values if no test data supplied
    if (is.null(test.dat))
      return(this.rlm$fitted.values)
    
    #make predictions on test data
    return(suppressWarnings(predict(this.rlm, newdata = data.frame(test.dat))))
  } else if (method == "errors.in.x") {
    this.mcr <- mcr::mcreg(x = this.y_hat, y = this.y, method.reg = "PaBa", method.ci = "analytical")
    
    if (is.null(test.dat))
      return(mcr::MCResult.getFitted(this.mcr)[, "y_hat"])
    
    return(mcr::MCResultAnalytical.calcResponse(this.mcr, test.dat)[, "Y"])
  }
}

#helper function to make predictions from repeated glmnet results
predict.repglmnet <- function(em, this_dat, resp_var, lambda_type, recalibrate) {
  all_coefs  <- getAllCoefs(xx = coef.glmnet(em$enet.model, s = em$selected.mod["lambda", lambda_type]),
                            remove_intercept = TRUE)
  
  #get the predictors and their coefficietns from the model
  #nz_coefs <-  getNZCoefs(coef.glmnet(em$enet.model, s = em$selected.mod["lambda", lambda_type]), 
  #                        remove_intercept = TRUE) 
  
  y_hat <- as.numeric(predict(em$enet.model, s = em$selected.mod["lambda", lambda_type], 
                              type = "response", newx = as.matrix(this_dat[, all_coefs])))
  
  if (!recalibrate)
    return(y_hat)
  
  return(recal(this.y = this_dat[, resp_var], this.y_hat = tmb_hat)) 
}


#helper function to calculate prediction error for continuous outcome
evalPred <- function(y, y_hat) {
  if (length(unique(y)) > 2 & length(unique(y)) <= sqrt(length(y)))
    stop("Can only handle continouous outcomes, looks like you have multinomial or ordinal")
  
  if (is.numeric(y) & length(unique(y)) > 2)
    return(c("R-squared" = cor(y, y_hat, use = "complete.obs")^2))
}

#function to compute performance (ppv, npv, opv) between 
#two binary variables
compBinPerf <- function(actual, pred) {
  tp = sum(pred == 1 & actual == 1) #true positive
  fp = sum(pred == 1 & actual == 0) #false positive
  tn = sum(pred == 0 & actual == 0) #true negative
  fn = sum(pred == 0 & actual == 1) #false negative
  
  ppv <- tp/(tp + fp)
  npv <- tn/(tn + fn)
  opv <- (tp + tn)/(tp + fp + tn + fn)
  return(c(tp = tp, fp = fp, tn = tn, fn = fn, ppv = ppv, npv = npv, opv = opv))
}

#helper function to calculate prediction error for binary outcomes
evalBinPred <- function(y, y_hat) {
  if (length(unique(y)) > 2)
    stop("Can only handle binary outcomes")
  
  #calculate AUROC, ROC, opv, ppv, npv
  # roc_obj <- pROC::roc(y, y_hat)
  # 
  # auc(roc_obj)
  # ## Area under the curve: 0.825
  # roc_df <- data.frame(TPR    = rev(roc_obj$sensitivities), 
  #                      FPR    = rev(1 - roc_obj$specificities), 
  #                      labels = roc_obj$response, 
  #                      scores = roc_obj$predictor)
  
  idxs = seq_len(length(y_hat))
  if (any(is.na(y_hat)))
    idxs = which(!is.na(y_hat))
  
  y     <- y[idxs]
  y_hat <- y_hat[idxs]
  
  rocr.pred <- ROCR::prediction(predictions = y_hat, labels = y)
  
  pred.mat <- dplyr::inner_join(ggplot2::fortify(ROCR::performance(rocr.pred, "tpr", "fpr")),
                                ggplot2::fortify(ROCR::performance(rocr.pred, "ppv", "npv")),
                                by = "Cutoff") %>%
    dplyr::inner_join(ggplot2::fortify(ROCR::performance(rocr.pred, "acc")), by = "Cutoff") %>%
    dplyr::select(-matches("Repetition.Number")) %>%
    dplyr::mutate(Specificity = 1 - False.positive.rate)
  
  #adding prevelance
  this.n              <- if (is.numeric(y)) length(y) else nrow(y)
  pred.mat$Prevalence <- sapply(as.vector(pred.mat$Cutoff), function(x) (sum(rocr.pred@predictions[[1]] >= x))/this.n)
  
  
  #adding tp, fp, tn, fn
  my.pred <- data.frame(do.call("rbind", lapply(pred.mat$Cutoff,
                                                function(co) {
                                                  return(c("Cutoff" = co, compBinPerf(actual = y, pred = ifelse(y_hat >= co, 1, 0))))
                                                })), stringsAsFactors = FALSE)
  
  
  
  #check that npv, ppv, and opv are idencitical
  num.not.equal.preds <- sum(my.pred$npv != pred.mat$Negative.predictive.value, na.rm = T) + 
    sum(my.pred$ppv != pred.mat$Positive.predictive.value, na.rm = T) + 
    sum(my.pred$opv != pred.mat$Accuracy, na.rm = T)
  if ( num.not.equal.preds > 0)
    warning("Accuracy metrics between Alex's code and ROCR don't match up!!!")
  
  if (num.not.equal.preds == 0)
    pred.mat <- pred.mat %>%
    dplyr::inner_join(dplyr::select(my.pred, one_of("tp", "fp", "tn", "fn", "Cutoff")), by = "Cutoff")
  
  ret.list <- list()
  #ret.list$pred.obj <- rocr.pred
  ret.list$pred.mat <- pred.mat
  ret.list$auc      <- ROCR::performance(rocr.pred, "auc")@y.values[[1]]
  
  return(ret.list)
}

#evaluate an enet model on the full data
#hyperparameters fit using repeated-cv
evalModel <- function(em, this_dat, resp_var = "MutLoad", mod_name = "model", 
                      lambda_type = "min", recalibrate = FALSE, ml_co = 200, font_size) {
  
  #setting up lists for storing results 
  plots       <- list()
  #models      <- list()
  performance <- list()
  coefs       <- list()
  
  #adding plots
  plots[1:4] <- em$plots[1:4]
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #evaluate (in-sample) the constrained linear regression model built to predict TMB
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  #get all predictors in the dataset
  all_coefs  <- getAllCoefs(xx = coef.glmnet(em$enet.model, s = em$selected.mod["lambda", lambda_type]),
                            remove_intercept = TRUE)
  
  #get the predictors and their coefficietns from the model
  nz_coefs <-  getNZCoefs(coef.glmnet(em$enet.model, s = em$selected.mod["lambda", lambda_type]), 
                          remove_intercept = TRUE) 
  
  coefs[["all"]]                <- all_coefs
  coefs[["coefs_nonZero_orig"]] <- nz_coefs
  
  
  tmb_hat <- as.numeric(predict(em$enet.model, s = em$selected.mod["lambda", lambda_type], 
                                type = "response", newx = as.matrix(this_dat[, all_coefs])))
  
  
  tmb_hat_re <- recal(this_dat[, resp_var], tmb_hat)  
  
  hat_dat <- data.frame(tmb           = this_dat[, resp_var],
                        #tmb.co        = ifelse(this_dat[, resp_var] >= log10(ml_co), 1, 0),
                        tmb.hat       = tmb_hat,
                        #tmb.hat.co    = ifelse(tmb_hat >= log10(ml_co), 1, 0),
                        tmb.hat.re    = tmb_hat_re) #,
  #tmb.hat.re.co = ifelse(tmb_hat_re >= log10(ml_co), 1, 0))
  
  
  performance[["cont"]] <- cor(hat_dat$tmb, hat_dat$tmb.co)
  performance[["bin"]]  <- compBinPerf(actual = hat_dat$tmb.co, 
                                       pred = hat_dat$tmb.hat.co)
  
  
  plots[[length(plots) + 1]] <- ggplot(hat_dat) + 
    geom_point(aes(y = tmb.hat, x = tmb.hat.re)) +
    geom_abline() +
    ylab("Mutational Load (log10 Missense Count") +
    xlab("Predicted Mutational Load (constrained model)") +
    ggtitle(paste0("Actual vs Predicted TMB R^2 = ", round(cor(hat_dat$tmb, hat_dat$tmb.hat)^2, 3))) +
    geom_vline(aes(xintercept = log10(ml_co)), size = 1.25, linetype = "dashed", col = "darkred") + 
    geom_hline(aes(yintercept = log10(ml_co)), size = 1.25, linetype = "dashed", col = "darkred") + 
    theme_bw(font_size) +
    xlim(.9, 3.25) + 
    ylim(.9, 3.25) + 
    coord_fixed() +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  #plots without vlines or hlines
  plots[[length(plots) + 1]] <- ggplot(hat_dat) + 
    geom_point(aes(y = tmb, x = tmb.hat)) +
    geom_abline() +
    ylab("Mutational Load (log10 Missense Count") +
    xlab("Predicted Mutational Load (constrained model)") +
    ggtitle(paste0("Actual vs Predicted TMB R^2 = ", round(cor(hat_dat$tmb, hat_dat$tmb.hat)^2, 3))) +
    #geom_vline(aes(xintercept = log10(ml_co)), size = 1.25, linetype = "dashed", col = "darkred") + 
    #geom_hline(aes(yintercept = log10(ml_co)), size = 1.25, linetype = "dashed", col = "darkred") + 
    theme_bw(font_size) +
    xlim(.9, 3.25) + 
    ylim(.9, 3.25) + 
    coord_fixed() +
    theme(plot.title = element_text(hjust = 0.5))
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  td <- cbind(this_dat, hat_dat)
  
  #make os and pfs models
  
  pfs_term <- "Surv(PFSIRC1, 1 - PFSIRC1.CNSR)"
  
  #creating the terms
  terms <- paste0(colnames(hat_dat), "*trt")
  
  #creating the vector specifying which terms to put in each models
  m <- list()
  m[[1]] <- 1
  m[[2]] <- 2
  m[[3]] <- 3
  m[[4]] <- 4
  m[[5]] <- c(1, 3)
  m[[6]] <- c(2, 4)
  names(m) <- sapply(m, function(this) paste0(terms[this], collapse = "_"))
  
  
  models[["pfs"]] <- lapply(m, function(term_idxs) {
    return(coxph(as.formula(paste0(pfs_term, " ~ ", paste0(terms[term_idxs], collapse = " + "))), data = td))
  })
  
  
  os_term <- "Surv(OS, 1 - OS.CNSR)"
  models[["os"]] <- lapply(m, function(term_idxs) {
    return(coxph(as.formula(paste0(os_term, " ~ ", paste0(terms[term_idxs], collapse = " + "))), data = td))
  })
  
  
  ret <- list(models = models, plots = plots, performance = performance, coefs = coefs)
  return(ret)
}


#helper function to build model of PFS and OS (seperately) where the covariates
#are the trt arm, as an interaction, and the tmb and tmb^hat
# clinRespModels <- function(this.dset, tmb, tmb.hat, tte, tte.cens, trt.arm) {
#   
# }

#evaluate the full process, creating out of sample predictions for the repeated folds
calcRepCvPred_enet <- function(ens, lambda_type = "min", run.parallel = TRUE, refit.mod = FALSE, recal.mod = FALSE,
                               recal.method = "lm", back.trans = NULL) {
  num.cores <- ifelse(run.parallel, parallel::detectCores() - 2, 1)
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #   some sanity checks
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  #1) are all of models fit using the same model family, which should be gaussian. There is no good
  #   reason why this shouldnt be the case.
  if (length(unique(table(sapply(ens$models, function(m) m$mod.fam)))) != 1)
    stop("more than one model type was used...something is very wrong!")
  #ok, we're good
  
  #2) the object ens.cont$outer.loop.folds contains the indices used for training in each fold, making sure
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
    cat(idx)
    tm <- ens$models[[idx]] #get the model at this idx
    
    #get the indices used for training
    train.rows <- tm$selected.rows
    test.rows  <- seq_len(nrow(ens$dset))[-train.rows]
    
    
    this.s <- tm$selected.mod["lambda", lambda_type] #get value of lambda based on selected lambda type
    
    #in case some feature are missing
    #this.dset = ens$dset
    this.coef = getAllCoefs(coef.glmnet(object = tm$enet.model, s = this.s))
    
    # missing.feat = colnames(this.dset)[!colnames(this.dset) %in% this.coef]
    # if (length(missing.feat) > 0) {
    #   zeroes <- matrix(0, nrow(this.dset), length(missing.feat))
    #   colnames(zeroes) <- missing.feat
    #   this.dset <- cbind(this.dset, zeroes)
    #   ens$dset <- this.dset
    # }
    # 
    ens$dset <- ens$dset[, this.coef]
    if (!refit.mod) {
      this.pred <- predict(object = tm$enet.model, s = this.s, newx = ens$dset[test.rows, ])
      
      #use back transform if response was transformed                                                            
      if (!is.null(back.trans))
        this.pred <- back.trans(this.pred)
      
      if (!recal.mod) 
        return(list(pred = this.pred, test.rows = test.rows))
      
      #if we do recalibrate, then fit take the predicted values,
      #refit, and using the refit (secondary model), make a prediction
      iss.pred <- predict(object = tm$enet.model, s = this.s, newx = ens$dset[train.rows, ])
      y.train  <- ens$response[train.rows, ]
      
      return(list(pred = recal(this.y = y.train, this.y_hat = iss.pred, method = recal.method, test.dat = data.frame(this.y_hat = this.pred)),
                  test.rows = test.rows))
      
    } else {
      #get selected features at chosen lambda
      this.feat <- names(getNZCoefs(coef.glmnet(tm$enet.model, s = this.s), remove_intercept = T))
      
      #get the  data
      this.md <- cbind(ens$response[train.rows, , drop = F], ens$dset[train.rows, this.feat])
      
      #build lm model
      this.mod <- glm(as.formula(paste0(names(ens$response), " ~ .")), data = data.frame(this.md, stringsAsFactors = FALSE))
      
      #make prediction
      this.pred <- predict(this.mod, newdata = data.frame(ens$dset[test.rows, ]))
      
      #use back transform if response was transformed                                                            
      if (!is.null(back.trans))
        this.pred <- back.trans(this.pred)
      
      
      return(list(pred = this.pred, test.rows = test.rows))
    }
  }, mc.cores = num.cores)
  
  #setting up matrices to store results
  perf.mat <- matrix(NA, nrow(ens$dset), max(mod.map$cv.rep))
  
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
  
  coefs <- mclapply(ens$models, function(m) return(getNZCoefs(x = coef.glmnet(m$enet.model, s = m$selected.mod["lambda", lambda_type]),
                                                              remove_intercept = TRUE)), mc.cores = num.cores)
  
  return(list(pred.mat = perf.mat, coefs = coefs, idxs = ens$outer.loop.folds))
}


#calculate the performance metrics given a matrix of out of sample predictions
evalEnsModel <-  function(ens, resp_var = "MutLoad", mod_name, lambda_type = "min", recalibrate = FALSE, ml_co = 200, font_size = 20) {
  #setting up lists for storing results 
  plots       <- list()
  models      <- list()
  performance <- list()
  coefs       <- list()
  
  
  cv.rep.preds <- calcRepCvPred_enet(ens, resp_var = "MutLoad", mod_name, lambda_type = "min", recalibrate = FALSE)
  perf.mat     <- cv.rep.preds$perf.mat
  
  label.fs = 7 #label font size
  
  perf.res <- data.frame(tmbHatCor = apply(perf.mat, 2, function(x) cor(x, ens$response[, resp_var])))
  
  med.perf <- median(perf.res$tmbHatCor)
  
  plots[[length(plots) + 1]] <- ggplot(perf.res) + 
    geom_density(aes(x = tmbHatCor, y = ..scaled..), fill = "darkgrey", col = "darkgrey") +
    theme_bw(font_size) +
    #geom_vline(aes(xintercept = med.perf), size = 1.25, linetype = "dashed") +
    #geom_label(label = paste0("median =\n", round(med.perf, 2)), 
    #           x = med.perf, y = .75, size = label.fs, fontface = "bold") +
    labs(title = "Distribution of Out-of-Sample Correlation",
         subtitle = paste0("Across ", PARAMS$num_repeats, " Repeats of ", PARAMS$outer_cv_folds, "-fold CV")) +
    xlab("Correlation") +
    ylab("Density") + 
    theme(plot.title = element_text(hjust = .5, face = "bold"),
          plot.subtitle = element_text(hjust = .5, size = rel(1.1)))
  
  
  #perf.refit.res <- apply(perf.refit.mat, 2, function(x) cor(x, ens$response[, resp_var]))
  
  #calculte npv, pvv, etc.
  perf.co.res <- t(apply(perf.mat, 2, function(x){
    x.co <- ifelse(x >= log10(ml_co), 1, 0)
    y <- ifelse(ens$response[, resp_var] >= log10(ml_co), 1, 0)
    
    return(compBinPerf(pred = x.co, actual = y))  
  }))
  
  performance[["contin"]] <- perf.res
  performance[["bin"]]    <- perf.co.res
  
  plot_dat <- data.frame(perf.co.res) %>%
    dplyr::select(ppv, npv, opv) %>%
    gather(key = "metric", value = "value")
  
  
  plots[[length(plots) + 1]] <- ggplot(plot_dat, aes(x = value)) + 
    geom_density(aes(fill = metric, col = metric), alpha = .6) + 
    #geom_histogram(aes(fill = metric, col = metric), alpha = .6)
    theme_bw(font_size) +
    labs(title = "Distribution of Out-of-Sample Binary Classifcation Accuracy",
         subtitle = paste0("Across ", PARAMS$num_repeats, " Repeats of ", PARAMS$outer_cv_folds, "-fold CV")) +
    xlab("Accuracy") +
    ylab("Density") + 
    theme(plot.title = element_text(hjust = .5, face = "bold"),
          plot.subtitle = element_text(hjust = .5, size = rel(1.1)))
  
  p.idx  <- which.min(abs(median(perf.res$tmbHatCor) - perf.res$tmbHatCor))
  
  models[["medModel"]] <- list()
  
  models[["medModel"]][["tmb.hat"]]    <- perf.mat[, p.idx]
  models[["medModel"]][["tmb.hat.co"]] <- ifelse(models[["medModel"]][["tmb.hat"]] >= log10(ml_co), 1, 0)
  models[["medModel"]][["tmb.hat.cor"]]     <- perf.res$tmbHatCor[p.idx]
  models[["medModel"]][["tmb.hat.co.conc"]] <- perf.co.res[p.idx, ]
  
  oss_plot_dat <- data.frame(tmb = ens$response[, resp_var],
                             med.tmb.hat = perf.mat[, p.idx])
  
  #plot of actual vs predicted for median
  plots[[length(plots) + 1]] <- ggplot(oss_plot_dat) + 
    geom_point(aes(y = tmb, x = med.tmb.hat)) +
    geom_abline() +
    ylab("Mutational Load (log10 Missense Count") +
    xlab("Predicted Mutational Load (constrained model)") +
    labs(title = paste0("Actual vs Predicted TMB R^2 = ", round(cor(oss_plot_dat$tmb, oss_plot_dat$med.tmb.hat)^2, 3)),
         subtitle = "Median Correlated Out-of-Sample Prediction") +
    geom_vline(aes(xintercept = log10(ml_co)), size = 1.25, linetype = "dashed", col = "darkred") + 
    geom_hline(aes(yintercept = log10(ml_co)), size = 1.25, linetype = "dashed", col = "darkred") + 
    theme_bw(font_size) +
    xlim(.9, 3.25) + 
    ylim(.9, 3.25) + 
    coord_fixed() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = rel(1.1)))
  
  #get the coefs
  coefs <- lapply(ens$models, function(m) return(getNZCoefs(x = coef.glmnet(m$enet.model, s = m$selected.mod["lambda", lambda_type]),
                                                            remove_intercept = TRUE)))
  
  ret <- list(models = models, plots = plots, performance = performance, coefs = coefs)
  return(ret)
}



#eval a pre-specifiied linear model with repeated-cv
#dset - dataset with all of the features to be used in model building
#response - response 
#outer.loop.folds - name list which contains training indices for each model, for each cv
calcRepCvPred_lm <- function(dset, response, outer.loop.folds, run.parallel = T) {
  num.cores <- ifelse(run.parallel, parallel::detectCores() - 2, 1)
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #we will loop over each repeated run of 5 fold cv, and for each one we calculate performance on
  #test set in two ways
  #   1) constrained coefficients from glmnet
  
  #for book keeping create a data frame which maps each model which cv repeat and fold it comes from
  mod.map <- data.frame(cv.rep = as.numeric(sapply(strsplit(names(outer.loop.folds), "_"), function(s) s[[2]])),
                        fold = as.numeric(sapply(strsplit(names(outer.loop.folds), "_"), function(s) s[[4]])))
  rownames(mod.map) <- names(outer.loop.folds)
  
  #parallelizing the code to make predictions from models
  preds <- parallel::mclapply(seq_len(nrow(mod.map)), function(idx){
    print(".")
    train.rows <- outer.loop.folds[[idx]]
    test.rows  <- seq_len(nrow(dset))[-train.rows]
    
    #make data frame of training data
    m.dat <- data.frame(cbind(dset[train.rows, , drop = FALSE], response[train.rows, , drop = FALSE]))
    
    #fit the model
    m <- lm(as.formula(paste0(colnames(response), "~ .")), data = m.dat)
    
    #get the predictions 
    m.pred <- predict(object = m, newdata = data.frame(dset[test.rows, ]))
    
    return(list(oss.pred = m.pred, test.rows = test.rows, coefs = summary(m)$coefficients))
  }, mc.cores = num.cores)
  
  #now put the predictions and coefs in their proper storage in the matrix
  perf.mat  <- matrix(NA, nrow(dset), max(mod.map$cv.rep))
  coefs     <- vector(mode = "list", length = length(outer.loop.folds))
  
  for (idx in seq_len(nrow(mod.map))) {
    cat(".")
    cv.rep <- mod.map$cv.rep[idx] #which cv.rep are we on
    
    test.rows  <- preds[[idx]]$test.rows
    
    #make sure none oe the test rows for this repeat of cv have results,
    #if they do we've double counted somehwere and have an error
    if (any(!is.na(perf.mat[test.rows, cv.rep])))
      stop("We've already calculated performance for at least one of the test indices...there is a bug somewhere")
    
    #evaluating on test data
    perf.mat[test.rows, cv.rep] <- preds[[idx]]$oss.pred
    
    coefs[[idx]] <- preds[[idx]]$coefs
    
    rm(cv.rep, test.rows) 
    #invisible(gc())
  }
  
  return(list(pred.mat = perf.mat, coefs = coefs, idxs = outer.loop.folds))
}


evalLmRepCv <-  function(dset, response, mod_name, outer.loop.folds, ml_co = 200, font_size = 20) {
  #setting up lists for storing results 
  plots       <- list()
  models      <- list()
  performance <- list()
  #coefs       <- list()
  
  
  
  cv.rep.preds <- calcRepCvPred_lm(dset = dset, response = response, outer.loop.folds = outer.loop.folds)
  perf.mat     <- cv.rep.preds$pred.mat
  
  label.fs = 7 #label font size
  
  perf.res <- data.frame(tmbHatCor = apply(perf.mat, 2, function(x) cor(x, response[, 1])))
  
  med.perf <- median(perf.res$tmbHatCor)
  
  plots[[length(plots) + 1]] <- ggplot(perf.res) + 
    geom_density(aes(x = tmbHatCor, y = ..scaled..), fill = "darkgrey", col = "darkgrey") +
    theme_bw(font_size) +
    #geom_vline(aes(xintercept = med.perf), size = 1.25, linetype = "dashed") +
    #geom_label(label = paste0("median =\n", round(med.perf, 2)), 
    #           x = med.perf, y = .75, size = label.fs, fontface = "bold") +
    labs(title = "Distribution of Out-of-Sample Correlation",
         subtitle = paste0("Across ", PARAMS$num_repeats, " Repeats of ", PARAMS$outer_cv_folds, "-fold CV")) +
    xlab("Correlation") +
    ylab("Density") + 
    theme(plot.title = element_text(hjust = .5, face = "bold"),
          plot.subtitle = element_text(hjust = .5, size = rel(1.1)))
  
  
  #perf.refit.res <- apply(perf.refit.mat, 2, function(x) cor(x, ens$response[, resp_var]))
  
  #calculte npv, pvv, etc.
  perf.co.res <- t(apply(perf.mat, 2, function(x){
    x.co <- ifelse(x >= log10(ml_co), 1, 0)
    y <- ifelse(ens$response[, resp_var] >= log10(ml_co), 1, 0)
    
    return(compBinPerf(pred = x.co, actual = y))  
  }))
  
  performance[["contin"]] <- perf.res
  performance[["bin"]]    <- perf.co.res
  
  plot_dat <- data.frame(perf.co.res) %>%
    dplyr::select(ppv, npv, opv) %>%
    gather(key = "metric", value = "value")
  
  
  plots[[length(plots) + 1]] <- ggplot(plot_dat, aes(x = value)) + 
    geom_density(aes(fill = metric, col = metric), alpha = .6) + 
    #geom_histogram(aes(fill = metric, col = metric), alpha = .6)
    theme_bw(font_size) +
    labs(title = "Distribution of Out-of-Sample Binary Classifcation Accuracy",
         subtitle = paste0("Across ", PARAMS$num_repeats, " Repeats of ", PARAMS$outer_cv_folds, "-fold CV")) +
    xlab("Accuracy") +
    ylab("Density") + 
    theme(plot.title = element_text(hjust = .5, face = "bold"),
          plot.subtitle = element_text(hjust = .5, size = rel(1.1)))
  
  
  
  
  # perf.co.refit.res <- t(apply(perf.refit.mat, 2, function(x){
  #   x.co <- ifelse(x >= log10(ml_10Mut_co), 1, 0)
  #   return(compBinPerf(pred = x.co, actual = dat.all$MutLoad_ten))  
  # }))
  # apply(perf.co.refit.res, 2, mean)
  # 
  
  p.idx  <- which.min(abs(median(perf.res$tmbHatCor) - perf.res$tmbHatCor))
  
  models[["medModel"]] <- list()
  
  models[["medModel"]][["tmb.hat"]]    <- perf.mat[, p.idx]
  models[["medModel"]][["tmb.hat.co"]] <- ifelse(models[["medModel"]][["tmb.hat"]] >= log10(ml_co), 1, 0)
  models[["medModel"]][["tmb.hat.cor"]]     <- perf.res$tmbHatCor[p.idx]
  models[["medModel"]][["tmb.hat.co.conc"]] <- perf.co.res[p.idx, ]
  
  oss_plot_dat <- data.frame(tmb = ens$response[, resp_var],
                             med.tmb.hat = perf.mat[, p.idx])
  
  #plot of actual vs predicted for median
  plots[[length(plots) + 1]] <- ggplot(oss_plot_dat) + 
    geom_point(aes(y = tmb, x = med.tmb.hat)) +
    geom_abline() +
    ylab("Mutational Load (log10 Missense Count") +
    xlab("Predicted Mutational Load (constrained model)") +
    labs(title = paste0("Actual vs Predicted TMB R^2 = ", round(cor(oss_plot_dat$tmb, oss_plot_dat$med.tmb.hat)^2, 3)),
         subtitle = "Median Correlated Out-of-Sample Prediction") +
    geom_vline(aes(xintercept = log10(ml_co)), size = 1.25, linetype = "dashed", col = "darkred") + 
    geom_hline(aes(yintercept = log10(ml_co)), size = 1.25, linetype = "dashed", col = "darkred") + 
    theme_bw(font_size) +
    xlim(.9, 3.25) + 
    ylim(.9, 3.25) + 
    coord_fixed() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = rel(1.1)))
  
  #get the coefs
  coefs <- lapply(ens$models, function(m) return(getNZCoefs(x = coef.glmnet(m$enet.model, s = m$selected.mod["lambda", lambda_type]),
                                                            remove_intercept = TRUE)))
  
  ret <- list(models = models, plots = plots, performance = performance, coefs = coefs)
  return(ret)
}


#evaluate the full process, creating out of sample predictions for the repeated folds
calcBootPred_enet <- function(ens, lambda_type = "min", run.parallel = TRUE, refit.mod = FALSE, recal.mod = FALSE, recal.method = "lm") {
  num.cores <- ifelse(run.parallel, parallel::detectCores() - 2, 1)
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #   some sanity checks
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  #1) are all of models fit using the same model family, which should be gaussian. There is no good
  #   reason why this shouldnt be the case.
  if (length(unique(table(sapply(ens$models, function(m) m$mod.fam)))) != 1)
    stop("more than one model type was used...something is very wrong!")
  #ok, we're good
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  #we will loop over each repeated run of 5 fold cv, and for each one we calculate performance on
  #test set in two ways
  #   1) constrained coefficients from glmnet
  
  
  preds <- mclapply(seq_len(length(ens$models)), function(idx){ 
    tm <- ens$models[[idx]] #get the model at this idx
    
    #get the indices used for training
    train.rows <- unique(tm$selected.rows)
    test.rows  <- seq_len(nrow(ens$dset))[-train.rows]
    
    
    this.s <- tm$selected.mod["lambda", lambda_type] #get value of lambda based on selected lambda type
    
    if (!refit.mod) {
      this.pred <- predict.glmnet(object = tm$enet.model, s = this.s, newx = ens$dset[test.rows, ])
      
      if (!recal.mod) 
        return(list(pred = this.pred, test.rows = test.rows))
      
      #if we do recalibrate, then fit take the predicted values,
      #refit, and using the refit (secondary model), make a prediction
      iss.pred <- predict.glmnet(object = tm$enet.model, s = this.s, newx = ens$dset[train.rows, ])
      y.train  <- ens$response[train.rows, ]
      
      return(list(pred = recal(this.y = y.train, this.y_hat = iss.pred, method = recal.method, test.dat = this.pred),
                  test.rows = test.rows))
      
    } else {
      #get selected features at chosen lambda
      this.feat <- names(getNZCoefs(coef.glmnet(tm$enet.model, s = this.s), remove_intercept = T))
      
      #get the  data
      this.md <- cbind(ens$response[train.rows, , drop = F], ens$dset[train.rows, this.feat])
      
      #build lm model
      this.mod <- glm(as.formula(paste0(names(ens$response), " ~ .")), data = data.frame(this.md, stringsAsFactors = FALSE))
      
      #make prediction
      this.pred <- predict(this.mod, newdata = data.frame(ens$dset[test.rows, ]))
      
      
      return(list(pred = this.pred, test.rows = test.rows))
    }
  }, mc.cores = num.cores)
  
  #setting up matrices to store results
  perf.mat <- matrix(NA, nrow(ens$dset), length(preds))
  
  for (idx in seq_len(length(preds))) {
    #cat(".")
    test.rows <- preds[[idx]]$test.rows
    
    #make sure none oe the test rows for this repeat of cv have results,
    #if they do we've double counted somehwere and have an error
    if (any(!is.na(perf.mat[test.rows, idx])))
      stop("We've already calculated performance for at least one of the test indices...there is a bug somewhere")
    
    
    perf.mat[test.rows, idx] <- preds[[idx]]$pred
    
    rm(idx, test.rows) 
  }
  
  coefs <- mclapply(ens$models, function(m) return(getNZCoefs(x = coef.glmnet(m$enet.model, s = m$selected.mod["lambda", lambda_type]),
                                                              remove_intercept = TRUE)), mc.cores = num.cores)
  
  return(list(pred.mat = perf.mat, coefs = coefs, idxs = ens$outer.loop.folds))
}


#helper functions to get abs of residuals to use as weights
getWeightsFromResid <- function(x, y, num.preds = 1, mod.fam = "gaussian") {
  if (mod.fam != "gaussian")
    stop("only gaussian family implemented")
  
  checked.inputs <- checkPredsResp(x = x, y = y, mod.fam = mod.fam)
  
  this.x <- checked.inputs$dset
  this.y <- checked.inputs$response
  
  cors      <- apply(this.x, 2, function(xx) cor(xx, this.y))
  top.genes <- names(sort(abs(cors), decreasing = T))
  
  use.genes <- top.genes[1:num.preds]
  
  this.dset <- cbind(this.x[, use.genes, drop = FALSE], data.frame(y = this.y))
  
  m <- glm(y ~ ., data = this.dset)
  r <- m$residuals
  return(list(residuals = abs(r), used.genes = use.genes))
}


#helper function to make table of performance metrics and
#threhsold where npv, ppv cross eachother, and where npv, pvv cross the
#selected cutoff values, first two cutoffs are npv, the second two are ppv
#agg.dat is an aggregated data frame of the results from repeated cv
#this.title
makeThresholdTables <- function(agg.dat, cutoffs = c(.8, .9, .8, .9), this.title = "Table", this.subt = "Cross-validated performance at select cutoffs") {
  p <- list()
  
  cos <- c()
  a <- agg.dat %>%
    dplyr::filter(Metric == "Negative predictive value") %>%
    dplyr::select(Metric, Performance_median, Cutoff_median)
  
  b <- agg.dat %>%
    dplyr::filter(Metric == "Positive predictive value") %>%
    dplyr::select(Metric, Performance_median, Cutoff_median)
  
  differences <- abs((a$Performance_median - b$Performance_median))
  co.idx <- which(differences == min(differences, na.rm = T))
  cos <- c(cos, b$Cutoff_median[co.idx])
  rm(co.idx, differences)
  
  diffs <- abs(a$Performance_median - cutoffs[1])
  co.idx <- which(diffs == min(diffs, na.rm = T))[1]
  cos <- c(cos, a$Cutoff_median[co.idx])
  rm(co.idx, diffs)
  
  diffs <- abs(a$Performance_median - cutoffs[2])
  co.idx <- which(diffs == min(diffs, na.rm = T))[1]
  cos <- c(cos, a$Cutoff_median[co.idx])
  
  diffs <- abs(b$Performance_median - cutoffs[3])
  co.idx <- which(diffs == min(diffs, na.rm = T))[1]
  cos <- c(cos, b$Cutoff_median[co.idx])
  
  diffs <- abs(b$Performance_median - cutoffs[4])
  co.idx <- which(diffs == min(diffs, na.rm = T))[1]
  cos <- c(cos, b$Cutoff_median[co.idx])
  rm(co.idx, diffs)
  
  x <- agg.dat[agg.dat$Cutoff_median %in% cos, ]
  xx <- x %>%
    dplyr::select(Cutoff_median, Metric, Performance_median) %>%
    spread(key = Metric, value = Performance_median) %>%
    remove_rownames() %>%
    dplyr::select(-one_of("tp", "fp", "tn", "fn")) %>%
    mutate(Cutoff_median = format(round(Cutoff_median, 2), 1)) %>%
    filter(!duplicated(Cutoff_median)) %>%
    column_to_rownames("Cutoff_median")
  
  
  
  xxx <- gather(xx) 
  
  xxx$cutoff <- rep(rownames(xx), ncol(xx))
  xxx$value <- as.character(format(round(xxx$value, 2), 1))
  xxx <- xxx[xxx$key != "False positive rate", ]
  
  #add the cutoffs to the table
  tmp       <- dplyr::filter(xxx, key == xxx$key[1]) 
  tmp$key   <- rep("Threshold", nrow(tmp))
  tmp$value <- tmp$cutoff
  
  xxx <- rbind(xxx, tmp)
  xxx <- xxx %>%
    dplyr::mutate(key = factor(key))
  
  #xxx$value <- format(xxx$value, digits = 3)
  
  xxx$key <- factor(xxx$key, levels(xxx$key)[rev(c(7, 4, 6, 5, 3, 2, 1))])
  
  p[[1]] <-  ggplot(xxx) + 
    geom_text(aes(y = key, x = cutoff, label = value), size = 7.5) +
    theme_bw(22) +
    labs(title = this.title,
         subtitle = this.subt) + 
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = .5),
          plot.subtitle = element_text(hjust = .5)) 
  return(p)
}  


#hlper function to get model annotation based on file name (dangerous, I now!)
getModAnno <- function(x, nc) {
  mod.name <- gsub("\\.rds", "", strsplit(x, "/")[[1]][[length(strsplit(x, "/")[[1]])]])
  mod.anno <- strsplit(mod.name, "_")[[1]]
  
  if (mod.anno[2] == "lmMod")
    mod.anno <- rev(mod.anno)
  
  #do we have an elasticnet model
  is.en <- length(grep("alpha",  mod.name)) != 0
  
  #set model type accordingly
  mod.type    <- ifelse(is.en, "ElasticNet", "pre-specified linear model")
  
  #det selected alpha if we have an elasticnet model
  sel.alpha   <- ifelse(is.en, 
                        gsub("alpha0", "\\.", mod.anno[which(sapply(lapply(mod.anno, function(z) grep("alpha", z)), function(zz) return(length(zz) != 0)))]), #automagically finding which entry in mod.anno has alpha 
                        NA)
  
  
  transform.str <- "log10"
  if (length(grep("sqrt", x)) != 0)
    transform.str <- "sqrt"
  
  feature.set = "all genes"
  if (length(grep("cinTwenty", x)) != 0)  {
    feature.set <- "CIN 20 gene signature"
  } else if (length(grep("cinTen", x)) != 0) {
    feature.set <- "CIN 10 gene signature"
  } else {
    feature.set <- paste0(feature.set, " (", nc," seleted)")
  }
  
  weights.str = ""
  if (length(grep("_weights_", x)) != 0)
    weights.str = " (weighted)"
  
  return(list(mod.type = mod.type, feature.set = feature.set, transform = transform.str, sel.alpha = sel.alpha, weights = weights.str))
}


#helper function to get scaled coefs
#inputs:
#1) features - names of the selected coefficients to re-scale
#2) dset - dataset of predictors (features)
#3) response - the response variable
#4) mod.fam - model family (e.g. gaussian, bionomial, etc.)
#5) scale.type - one of "IQR" or "SD", which scales the (refit) betas by either IQR or SDs.
#
#outputs:
#1) scaled.betas - vector of selected features refit to full dataset and scaled as specified
#note that these rescaled parameters are invalid for inference or further model testing. These
#betas can be used solely to compare the relative importance of features
getComparableCoefs <- function(features, dset, response, mod.fam = "gaussian", scale.type = "IQR") {
  scale.fun <- ifelse(scale.type == "IQR", function(x) return(quantile(x, .75) - quantile(x, .25)), function(x) sd(x))
  scales    <- apply(dset[, features], 2, scale.fun)
  
  md   <- data.frame(dset[, features], stringsAsFactors = FALSE, check.names = FALSE)
  md$y <- as.vector(response)
  
  res       <- glm(y ~ ., md, family = mod.fam)
  new.coefs <- getNZCoefs(as.matrix(res$coefficients))
  #need to remove some weird characters that happen sometimes
  names(new.coefs) <- gsub("\\`", "", names(new.coefs))
  
  if (!all(names(new.coefs) %in% names(scales)))
    stop("At least one of the features to be scaled does not have a scaling factor. This should never happen")
  
  ret <- new.coefs*scales[names(new.coefs)]
  
  #sort in order of decreasing absolute value
  ret <- ret[order(abs(ret), decreasing = T)]
  return(ret)
}
# Depends On: rep_grpreg


# data = mod_dat
#                x = c("TRT2", clin, INT_cat, Mstage, INT_Mstage, Region, INT_Region, bm, INT_bm, Ind, INT_Ind)
#                y = "PRIMRES_na"
#                group = grp
#                penalty = "grLasso"
#                lambda_seq = NULL
#                nreps_inner    = 2 #inner loop of cv to fit lambda
#                nfold_inner         = 5  #inner loop of cv to fit lambda
#                nreps_outer    = 3 #inner loop of cv to fit lambda
#                nfold_outer         = 4  #inner loop of cv to fit lambda
#                # nlambda        = 100, #number of lambdas to evaluate
#                # lambda.type    = "min", #which lambda to select from based on inner loop cv to apply to full model within each outer loop (not full-full model if there's outer loop present)
#                # num_outer_rep  = 100, #number of repeats for outer loop for performance estimation
#                # outer_cv_nfolds = 5,  #number of folds for outer loop of cross-validation
#                font_size   = 18
#                # run.parallel    = T,
#                ncores = 1             # response = surv.resp,
#                # cont_X = str_subset(colnames(dat), "_Q"),
#                alpha = 0.95
#                verbose = T
#                family = "binomial"
#                penalty = "grLasso"



ens_grpreg <- function(data = PARAMS$dat, x = PARAMS$x, y = PARAMS$y, group = PARAMS$group, family = "binomial", penalty = "grLasso",
                             alpha = 0.95, nfold_inner = 5, nfold_outer = 10, nreps_inner = PARAMS$nreps_inner,
                             nreps_outer = PARAMS$nreps_outer, lambda_seq = NULL, ncores = 1) {
  
  # function(data, x, y, group = NULL, family = "binomial", penalty = "grLasso",
  #                            alpha = 0.95, nfold_inner = 5, nfold_outer = 10, nreps_inner = 10,
  #                            nreps_outer = 10, lambda_seq = NULL, ncores = detectCores() - 1) {


    # Extract the predictor and response variables
  # X <- data[, x, drop = FALSE]
  # y <- data[, y, drop = FALSE]
  
  # Set up the parallel cluster
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  
  # Define the function for parallel computation
  preds <- foreach(
    s = 1:nreps_outer,
    .inorder = F,
    .packages = c("grpreg", "parallel"),
    .combine = cbind,
    .verbose = T,
    .export = c("rep_grpreg")) %dopar% {
    
    pred_matrix <- matrix(NA, nrow = nrow(data), ncol = 1)  
    # Split the data into training and test sets for the outer loop
    # for (j in fun) {
    #   cat("#### Starting REPEAT ", j)
      
    folds_outer <- sample(rep(1:nfold_outer, length.out = nrow(data)))
    cv_errors_outer <- rep(0, nfold_outer)
    
    
    for (i in 1:nfold_outer) {
      train <- data[folds_outer != i, ]
      test <- data[folds_outer == i, ]
 
      # Perform the inner loop of cross-validation on the training set
      models <- rep_grpreg(data = train, x = x, y = y, group = group, family = family, penalty = penalty,
                           alpha = alpha, nfold = nfold_inner, nreps = nreps_inner, lambda_seq = lambda_seq,
                           ncores_rep = 1)
     
      full_model <- models$full_model
  
      # Predict the response variable for each repeat of the test set
      # for (j in 1:nreps_outer) {
        pred_matrix[folds_outer == i, 1] <- predict(full_model, X = as.matrix(test[, x, drop = FALSE]),
                                                    type = "link", lambda = models$best_lambda)
        # cat("Finished Fold ", i)
    }
    # cat("#### Finished REPEAT ", fun)
      # }
    
    # Compute the mean test error and prediction matrix across the outer loop
    # mean_cv_error <- mean(cv_errors_outer)
    # mean_pred_matrix <- apply(pred_matrix, 1, mean)
    
    # Return the mean test error and prediction matrix
    return(pred_matrix)
    
    
  }
  
  
  # Run the function in parallel for each repeat
    
    # results <- mclapply(1:nreps_outer, outer_rep_grpreg_parallel, mc.cores = ncores)
  
  # Stop the parallel cluster
  
  stopCluster(cl)
  
  # preds <- do.call(cbind, results)
  
  # Extract the mean test errors and prediction matrices from the results
  
  # mean_cv_errors <- sapply(results, function(result) result$mean_cv_error)
  # mean_pred_matrices <- lapply(results, function(result) result$mean_pred_matrix)
  # 
  # Return the mean test errors and prediction matrices
  
  return(preds)
  }

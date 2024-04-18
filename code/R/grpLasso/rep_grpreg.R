# grpreg
# parallel


library(grpreg)
library(parallel)

# Define repeated_grpreg function
rep_grpreg <- function(data, x, y, group = NULL, family = "binomial", penalty = "grLasso",
                       alpha = 0.95, nfold = 5, nreps = 10, lambda_seq = NULL,
                       ncores_rep = detectCores() - 1) {

  
  # Subset data to x and y variables
  x_data <- data[, x, drop = F]
  y_data <- data[, y, drop = F]
  
  # Check for missing data in input dataframe
  if (any(is.na(x_data)) || any(is.na(y_data))) {
    stop("Input dataframe contains missing data")
  }
  
  # Choose lambda sequence if not provided
  if (is.null(lambda_seq)) {
    fit_full <- grpreg(x_data, y_data, family = family, group = group,
                       penalty = penalty, alpha = alpha)
    lambda_seq <- fit_full$lambda
  }

  # Initialize objects
  models <- vector("list", nreps)
  cv_errors <- vector("list", nreps)
  
  # Parallelize fitting and cross-validation
  if(ncores_rep > 1){
  cl <- makeCluster(ncores_rep)
  }
  # clusterExport(cl, c("x_data", "y_data", "group", "family", "nfold", "lambda_seq"))
  # registerDoParallel
  
  # Check the number of cores being used by the parallel backend
  # nc <- getOption("cl.cores")
  # cat("Using", nc, "cores\n")
  
  results <- mclapply(1:nreps, function(i) {
    # fit <- grpreg(x_data, y_data, family = family, group = group, penalty = penalty,
                  # alpha = alpha, lambda = lambda_seq)
    cv <- cv.grpreg(x_data, y_data, family = family, group = group, penalty = penalty,
                  alpha = alpha, nfold = nfold, lambda = lambda_seq)
    list(cv_errors = cv$cve)
  }, mc.cores = ncores_rep)
  
  if(ncores_rep > 1){
  stopCluster(cl)
  }
  
  # Store models and cross-validation errors
  for (i in 1:nreps) {
    # models[[i]] <- results[[i]]$model
    cv_errors[[i]] <- results[[i]]$cv_errors
  }
  
  # Combine cross-validation errors across repeats
  cv_errors_df <-do.call(rbind, cv_errors)
  median_cv_errors <- apply(cv_errors_df, 2, median)
  
  # Choose best lambda and final model
  
  best_lambda <- lambda_seq[which.min(median_cv_errors)]
  # final_model <- models[[which.min(mean_cv_errors)]]
  
  # Return list of results
  cat("Finished")
  
  list(lambda_seq = lambda_seq, cv_errors = cv_errors_df, best_lambda = best_lambda, full_model = fit_full, train = data)
  
}

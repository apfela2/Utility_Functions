# Adapt built-in glinternet.cv function to work with discrete data
discrete_glinternet_cv <- function (X, Y, numLevels, nFolds = 10, lambda = NULL, nLambda = 50, 
                                    lambdaMinRatio = 0.01, interactionCandidates = NULL, screenLimit = NULL, 
                                    family = c("gaussian", "binomial"), tol = 1e-05, maxIter = 5000, 
                                    verbose = FALSE, numCores = 1, cens, cat_X, col_contX,
                                    intervals = min(X[[ Y ]]) : (max(X[[ Y ]])), time_effect = T) 
{
  
  thisCall = match.call()
  family = match.arg(family)
  n = length(X[[Y]])
  pCat = sum(numLevels > 1)
  pCont = length(numLevels) - pCat
  # stopifnot(n == nrow(X), pCat + pCont == ncol(X), family ==
  #             "gaussian" || family == "binomial")
  folds = sample(rep(1:nFolds, ceiling(n/nFolds)), n, replace = FALSE)
  X$folds = folds
  discrete = disc(dat = X, cens = cens, surv.time = Y, intervals = intervals)
  folds = discrete$folds
  disc_genes = discrete[,col_contX]
  if (time_effect) {
    discrete$timeInt <- as.numeric(as.character(discrete$timeInt))
    disc_genes <- cbind(disc_genes, discrete$timeInt)
    numLevels <- c(numLevels, 1)
  }
  
  Y_disc = discrete$y
  X_disc = cbind(discrete[[cat_X]], disc_genes)
  fullfitted = glinternet(X = X_disc, Y = Y_disc, numLevels = numLevels, lambda = lambda, nLambda = nLambda, 
                          lambdaMinRatio = lambdaMinRatio, interactionCandidates = 1, screenLimit = screenLimit, family = family, 
                          tol = tol, maxIter = maxIter, verbose = verbose, numCores = numCores)
  if (verbose) {
    cat("\n Done fit on all data\n")
  }
  lambda = fullfitted$lambda
  nlambda = length(lambda)
  compute_loss = function(y, yhat, family) {
    if (family == "gaussian") {
      return(sum((y - yhat)^2)/(2 * length(y)))
    }
    yhat = sapply(yhat, function(x) min(max(1e-15, x), 1 - 
                                          1e-15))
    -(t(y) %*% log(yhat) + t(1 - y) %*% log(1 - yhat))/length(y)
  }
  loss = matrix(0, nFolds, nlambda)
  X_disc = as.matrix(X_disc)
  
  for (fold in 1:nFolds) {
    testIndex = (folds == fold)
    trainIndex = !testIndex
    fitted = glinternet(X = X_disc[trainIndex, , drop = FALSE], Y = Y_disc[trainIndex], 
                        numLevels = numLevels, lambda = lambda, nLambda = nlambda, lambdaMinRatio = lambdaMinRatio,
                        interactionCandidates = interactionCandidates, 
                        screenLimit = screenLimit, numToFind = NULL, family = family, tol = tol, maxIter = maxIter, 
                        verbose = verbose, numCores = numCores)
    YtestHat = predict(fitted, X_disc[testIndex, , drop = FALSE], 
                       "response")
    loss[fold, ] = apply(YtestHat, 2, function(yhat) compute_loss(Y_disc[testIndex], 
                                                                  yhat, family))
    if (verbose) {
      cat("\n Done fold", fold, "\n")
    }
  }
  cv = apply(loss, 2, mean)
  cvStd = apply(loss, 2, sd)
  bestIndex1Std = which(cv <= min(cv) + cvStd[which.min(cv)])
  bestIndex = which.min(cv)
  if (length(bestIndex1Std) == nlambda) {
    lambdaHat1Std = lambdaHat = lambda[bestIndex]
  }
  else {
    lambdaHat1Std = lambda[bestIndex1Std[1]]
    lambdaHat = lambda[bestIndex]
  }
  output = list(call = thisCall, glinternetFit = fullfitted, 
                fitted = fullfitted$fitted[, bestIndex], activeSet = fullfitted$activeSet[bestIndex], 
                betahat = fullfitted$betahat[bestIndex], lambda = lambda, 
                lambdaHat = lambdaHat, lambdaHat1Std = lambdaHat1Std, 
                cvErr = cv, cvErrStd = cvStd, family = family, numLevels = numLevels, 
                nFolds = nFolds)
  class(output) = "glinternet.cv"
  return(output)
}
# checkPredsResp #helper function used in only modelling functions that checks whether input dset, input response, and mod.fam are of the correct types # the intended use-case of this function is within a glmnet wrapper function
# makeModelReadyData
# makeModMatrix
# dataFrameToModMatrix
# filterBadVars #helper function to filter variables that either have nas, zero variance, or too few members in the minor class #mc.co is the cutoff for the smallest minor class size, will default to max(3, .01*nrow(x))
# mapFactorNames #helper function to map factor names between the data.frame and the model matrix
# check.features #helper function to check each feature to make sure there are sufficient observations for each class for factors/characters the minor class class has to have more than min.cat occurences for each class, and for numeric features there must be more unique values than min.num, a general cutoff is sqrt(N) #function assumes that x in N x p matrix where N is number of observations, and p is number of features
# get.high.cor #get pairs of variables that have high correlation, above cutoff cor.tol #the Nx3 table of highly correlated variables is returned, as will as the correlation matrix (optional)
# cluster.genes #helper function that clusters genes that are highly correlated above a given cutoff, and for each cluster (connected component) calculates VIF, selecting the gene to keep as the one with the highest VIF #returns a vector of genes that can be removed. These are genes that are in a connected component where at least one pair of correlations is geater than cor.thresh, and these genes do not have highest VIF
# num.empty #function that counts number of empty as either NA or ""



#helper function used in only modelling functions that
#checks whether input dset, input response, and mod.fam are
#of the correct types
#the intended use-case of this function is within a glmnet wrapper function
checkPredsResp <- function(x, y, mod.fam) {
  #checks for outcome type
  if (mod.fam == "cox" & !is.Surv(y))
    stop("selected model family is ", mod.fam, " , however the response is not a Surv object. Please create and supply a Surv object")
  
  if (is.Surv(y) & mod.fam != "cox") {
    mod.fam = "cox"
    warning("selected model family is ", mod.fam, ", however the response is a Surv object. The mod.fam will be changed to \"cox\" ")
  }
  
  if (class(y) %in% c("data.frame", "matrix") && ncol(y) > 1)
    stop("your response matrix has more than one column, and family is not cox...we don't support response matrices just yet")
  
  response.var <- NA
  if (class(y) %in% c("data.frame", "matrix")) {
    if (!is.null(colnames(y)))
      response.name <- colnames(y)[1]
    
    response.var = as.vector(y[, 1])
  } else if (!class(y) %in% c("data.frame", "matrix", "Surv")) {
    response.var <- matrix(y)
    colnames(response.var) <- "outcome"
  } else if (class(y) == "Surv") {
    response.var <- y
  } else {
    stop("Function checkPredAndResp encountered model type it doesn't know how to deal with. Shouldnt happen!")
  }
  
  #if (!is.data.frame(x))
  #  x <- data.frame(x)
  
  if (!is.matrix(x))
    stop("the input data set must be a matrix, try using makeModelReadyData() first")
  
  return(list(dset = x, response = response.var, mod.fam = mod.fam))
}



#high level function to make model ready data frame
makeModelReadyData <- function(dset, verbose = TRUE) {
  filt.res <- filterBadVars(x = dset, verbose = verbose)
  #removed  <- filt.res$removed
  #dset     <- filt.res$cleanDset
  
  converted.dset <- dataFrameToModMatrix(filt.res$cleanDset)
  #dset.mod.full  <- converted.dset$mod.mat
  #fac.map        <- converted.dset$fac.map
  
  return(list(dset = converted.dset$mod.mat,
              factor.map = converted.dset$fac.map,
              removed.vars = filt.res$removed))
}

#helper function for creating model matrix
makeModMatrix <- function(x) {
  classes <- sapply(x, class)
  
  if (any(classes == "character"))
    stop("Convert your characters to factors before running this function")
  
  #no work to do if there are no factors
  if (any(classes == "factor") == FALSE)
    return(as.matrix(x))
  
  num.mat <- as.matrix(x[, classes != "factor", drop = FALSE])
  #fac.mat <- model.matrix( ~ . -1, x[, classes == "factor", drop = FALSE])
  fac.mat <- model.matrix( ~ ., x[, classes == "factor", drop = FALSE])
  
  int.idx <- which(colnames(fac.mat) == "(Intercept)")
  if (length(int.idx) == 1)
    fac.mat <- fac.mat[, -int.idx, drop = F]
  
  ret           <- cbind(fac.mat, num.mat)
  colnames(ret) <- make.names(colnames(ret))
  
  invisible(gc())
  return(ret)
}

#takes in a data.frame
#converts to a modelMatrix and keeps 
#track of which factor levels were set to which values
dataFrameToModMatrix <- function(x) {
  if (class(x) != "data.frame")
    stop("this function requires a data frame; we should never get here")
  
  #get clases
  classes <- sapply(x, class) 
  #convert characters to factors
  if (any(classes == "character"))
    x    <- dplyr:::mutate_if(x, is.character, factor) 
  classes <- sapply(x, class) 
  
  #glmnet requires a data.frame, creating the model.frame here
  dset.mod.full           <- makeModMatrix(x)
  colnames(dset.mod.full) <- make.names(colnames(dset.mod.full))
  
  fac.map  <- mapFactorNames(dat = colnames(x), mod.dat = colnames(dset.mod.full), classes = classes) 
  #the factor names get all messed up, creating a mapping between the model.matrix names and 
  
  return(list(mod.mat = dset.mod.full, fac.map = fac.map))
}

#helper function to filter variables that either have nas, zero variance, or too few members in the minor class
#mc.co is the cutoff for the smallest minor class size, will default to max(3, .01*nrow(x))
filterBadVars <- function(x, verbose = TRUE, mc.co = NULL) {
  rns <- rownames(x)
  
  #remove columns with nas
  nas.ct <- sapply(x, num.empty) #summarise_each(x, funs(sum(is.na(.))))
  nas.rm <- colnames(x)[which(nas.ct > 0)]
  x      <- dplyr:::select(x, -one_of(nas.rm))
  if (length(nas.rm) > 0 & verbose)
    warning("Removed the variables below because they have missing (NA or empty) values:\n", paste(nas.rm, collapse = ", "))
  
  #remove columns with near zero variance
  nzv  <- colnames(x)[caret::nearZeroVar(x)]
  x    <- dplyr:::select(x, -one_of(nzv))
  if (length(nzv) > 0 & verbose) 
    warning("Removed the variables below, because they have variance near zero:\n", paste(nzv, collapse = ", "))
  
  #remove categorical columns whose minor class has too few vals
  ct.co <- mc.co
  if (is.null(ct.co))
    ct.co <- max(3, .01*nrow(x))
  
  classes <- sapply(x, class) #as.character(summarise_each(x, funs(class(.))))  
  
  #for the factors and character vars count the number of occurences in the minor class
  #minClass.ct <- summarise_each(x[, classes %in% c("factor", "character")], funs(min(table(.))))
  minClass.ct <- sapply(x[, classes %in% c("factor", "character"), drop = FALSE], function(this) min(table(this)))
  
  too.few <- names(minClass.ct)[which(minClass.ct < ct.co)]
  x       <- dplyr:::select(x, -one_of(too.few))
  
  if (length(too.few) > 0 & verbose)  
    warning("Removed the variables below, because the minor class has fewer than ", ct.co," observations:\n", paste(too.few, collapse = ","))
  
  rownames(x) <- rns
  
  removed = list(nas = nas.rm, zeroVariance = nzv, smallMinorClass = too.few)
  return(list(cleanDset = x, removed = removed))
}


#helper function to map factor names between the data.frame and the model matrix
mapFactorNames <- function(dat, mod.dat, classes) {
  if (sum(classes == "factor") == 0)
    return("")
  
  factors = dat[classes == "factor"]
  
  mod.fac = do.call("c", lapply(factors, function(f){ res <- grep(f, mod.dat, value = T)
  ff <- rep(f, length(res))
  names(ff) <- res
  return(ff)
  }))
  return(mod.fac)
}

#helper function to check each feature to make sure there are
#  sufficient observations for each class
#for factors/characters the minor class class has to have
#  more than min.cat occurences for each class, and for numeric
#  features there must be more unique values than min.num, a general
#  cutoff is sqrt(N)
#function assumes that x in N x p matrix where N is number of
#   observations, and p is number of features
check.features <- function(x, min.cat = 3, min.num = NULL) {
  if (is.null(min.num))
    min.num = round(sqrt(nrow(x)))
  
  classes = sapply(x, class)
  cat     = which(classes %in% c("factor", "character"))
  num     = seq_len(length(classes))[-cat]
  
  cat.occ = lapply(x[, cat], table)
  num.occ = lapply(x[, num], len.unq)
  
  rem = c(sapply(cat.occ, function(this) any(this < min.cat)),
          sapply(num.occ, function(this) any(this < min.num)))
  
  return(names(which(rem == TRUE)))
}

#get pairs of variables that have high correlation, above cutoff cor.tol
#the Nx3 table of highly correlated variables is returned, as will as the
#correlation matrix (optional)
get.high.cor <- function(x, cor.tol = .85, return.cor = FALSE) {
  classes <- NA
  
  if (is.matrix(x)) {
    classes <- apply(x, 2, class)
  } else {
    classes <- sapply(x, class)
  }
  xc <- coop::pcor(x[, classes %in% c("integer", "numeric"), drop = FALSE])
  
  #xc = cor(x)
  
  diag(xc) = 0
  xc[lower.tri(xc)] = 0
  
  high.cor     = data.frame(which(xc >= cor.tol, arr.ind = T), row.names = NULL, stringsAsFactors = FALSE)
  
  high.cor$rn = rownames(xc)[high.cor$row]
  high.cor$cn = colnames(xc)[high.cor$col]
  high.cor$cor = xc[which(xc >= cor.tol)]
  
  
  high.cor = dplyr:::select(high.cor, -one_of("row", "col"))
  
  
  return(list(high.cor.tbl = high.cor, cor.mat = ifelse(return.cor, xc, NA)))
}

#helper function that clusters genes that are highly correlated above a given cutoff, and 
#for each cluster (connected component) calculates VIF, selecting the gene to keep as the one with the highest VIF
#returns a vector of genes that can be removed. These are genes that are in a connected component where at least one
#pair of correlations is geater than cor.thresh, and these genes do not have highest VIF
cluster.genes <- function(dat, cor.thresh = .85) {
  cm <- coop::pcor(dat)
  
  cm[which(lower.tri(cm))] <- 0
  diag(cm) <- 0
  cm[cm < cor.thresh] <- 0
  
  rc <- which(cm > cor.thresh, arr.ind = T)
  rownames(rc) <- NULL
  rc[, "row"] <- rownames(cm)[as.numeric(rc[, "row"])]
  rc[, "col"] <- colnames(cm)[as.numeric(rc[, "col"])]
  
  g <- igraph::graph.edgelist(rc, directed = FALSE)
  
  g.clust <- igraph::clusters(g)
  
  cluster.size <- g.clust$csize
  cluster.membership <- g.clust$membership
  names(cluster.membership) <- igraph::V(g)$name
  
  to.keep = c()
  for(cluster in 1:length(cluster.size)){
    cluster.genes  = names(which(cluster.membership == cluster))
    cluster.mat    = dat[, colnames(dat) %in% cluster.genes]
    cluster.cormat = cor(cluster.mat)
    
    #if determinant is near zero, means we can't invert so take gene closest to median as the
    #representative for that cluster
    if (det(cluster.cormat) <= .Machine$double.eps) {
      #if we can't invert the cor matrix and caclulate variance inflation factor,
      #we'll just take the column that has the highest total correlation
      cor.sums <- apply(cluster.cormat, 2, function(this) sum(abs(this)))
      to.keep = append(to.keep, names(sort(cor.sums, decreasing = T))[1])
    } else {
      cluster.vif = diag(solve(cluster.cormat))
      to.keep = append(to.keep, names(sort(cluster.vif, decreasing = T))[1])
    }
  }
  
  to.remove <- unique(as.vector(rc))
  to.remove <- to.remove[-which(to.remove %in% to.keep)]
  
  return(to.remove)
}

#function that counts number of empty as either NA or ""
num.empty <- function(x) {return(sum(is.na(x) | x == ""))}



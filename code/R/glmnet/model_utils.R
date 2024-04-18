# sparseClust # sparse clusteirng, based on BIC
# plot.sparseClust #seperate plotting function to plot spare clustering results #sc.res is object returned from sparse clustering
# eucDist #helper function to calculate euclidean distance
# logPlus
# eigthRoot
# normality.test #helper function to test for normality
# rsquared #computes r-squared for a model fit
# non.zero
# len.unq
# secondDerivative


#  sparse clusteirng, based on BIC
sparseClust = function(x = t(sig.scores), method = 'ward.D2', maxK = NULL, 
                       makePlots = TRUE, agg.func = median, less.sparse = TRUE) {
  #do the clustering calculating bics
  
  ##cluster rows of x
  clust = hclust(dist(x), method = method)
  ##if no K is specified, search up to sqrt(K), which is a high discretization rule of thumb
  ##for the # of bins (usually lower than this)
  if (is.null(maxK))  
    maxK = round(sqrt(nrow(x))) 
  
  bics = vector('numeric', maxK)
  for (k in 1:maxK) {
    pred = cutree(clust, k)
    
    if (k == 1) {
      bics[k] = mean(apply(x, 2, function(this) BIC(lm(this ~ 1))))
    } else {
      bics[k] = mean(apply(x, 2, function(this) BIC(lm(this ~ as.factor(pred)))))
    }
    
  }
  rm(k)
  
  
  bestK = which.min(bics)#which(bics == min(bics))
  ##Adjust to be one higher than the optimum in order to stay conservative 
  ##in reporting representative networks when the basin of attraction is broad and shallow
  
  if (less.sparse == TRUE)
    bestK = bestK + 1
  
  classification = cutree(clust, bestK)
  
  
  
  #find euclidian distances between each cluster and the mediods
  #return object that has the representative element for each cluster
  centers = vector(mode = "list", length = bestK)
  names(centers) = 1:bestK
  
  for (k in 1:bestK) {
    tx = x[classification == k, , drop = FALSE]
    
    if (nrow(tx) == 1) {
      centers[[k]] = rownames(tx)
      next;
    }
    
    clust.center = apply(tx, 2, agg.func)
    
    dist.center = apply(tx, 1, function(this) eucDist(this, clust.center))
    centers[k]  = names(which.min(dist.center))
  }
  
  #making plot of K vs BIC
  bicdf = data.frame(k = 1:maxK, BIC = bics)
  
  # if (makePlots == TRUE) {
  #   print(ggplot(data = bicdf, aes(x = k, y = BIC)) + 
  #           geom_point() +
  #           xlab("# of clusters") + 
  #           ylab("Mean of variable BIC") +
  #           ggtitle("Number of Clusters vs. Mean Variable BIC") + 
  #           theme_bw() + 
  #           theme(plot.title = element_text(size = 14)))
  #   
  #   
  #   
  #   #plot the heatmap using ComplexHeatmap
  #   #use colorspace to pick colors
  #   cols = rainbow_hcl(n = bestK)
  #   names(cols) = names(centers)
  #   
  #   ha.df           = data.frame(cluster = classification)
  #   rownames(ha.df) = names(classification)
  #   
  #   csl = list(cluster = cols)
  #   
  #   ha = HeatmapAnnotation(df = ha.df, col = csl, which = "row")
  #   ht = Heatmap(x, show_row_names = ifelse(nrow(x) < 50, TRUE, FALSE), show_column_names = FALSE, name = "Value")
  #   
  #   draw(ha + ht)
  # }
  # 
  return(list(classification = classification, 
              K = bestK,
              clust = clust,
              centers = centers,
              bicdf = bicdf))
}

#seperate plotting function to plot spare clustering results
#sc.res is object returned from sparse clustering
plot.sparseClust <- function(x, sc.res) {
  print(ggplot(data = sc.res$bicdf, aes(x = k, y = BIC)) + 
          geom_point() +
          xlab("# of clusters") + 
          ylab("Mean of variable BIC") +
          ggtitle("Number of Clusters vs. Mean Variable BIC") + 
          theme_bw() + 
          theme(plot.title = element_text(size = 14)))
  
  if (ncol(x) >= 5000) {
    sds <- unlist(parallel::mclapply(x, sd))
    x <- x[, order(sds, decreasing = T)[1:1000]]
  }
  
  #plot the heatmap using ComplexHeatmap
  #use colorspace to pick colors
  cols = rainbow_hcl(n = sc.res$K)
  names(cols) = names(sc.res$centers)
  
  ha.df           = data.frame(cluster = sc.res$classification)
  rownames(ha.df) = names(sc.res$classification)
  
  csl = list(cluster = cols)
  
  ha = HeatmapAnnotation(df = ha.df, col = csl, which = "row")
  ht = Heatmap(x, show_row_names = ifelse(nrow(x) < 50, TRUE, FALSE), show_column_names = FALSE,
               name = "Value", use_raster = TRUE)
  
  draw(ha + ht)
  
}

#helper function to calculate euclidean distance
eucDist = function(x, y) {
  return(sqrt(sum((x - y)^2)))
}


logPlus   <- function(x) {return(log(x + 1))}
eigthRoot <- function(x) {return(x^(1/8))}

#helper function to test for normality
normality.test <- function(dat, trans.func = logPlus, sig.cut = .05) {
  classes       = as.chracter(summarise_each(dat, funs(class(.))))
  sw.res        = vector(mode = "list", length = sum(classes != "character"))
  names(sw.res) = names(which(classes != "character"))
  
  #func = ifelse(trans.func == "log", c(logPlus), c(eigthRoot))[[1]]
  
  for (nm in names(sw.res)) {
    this.x = dset[, nm]
    
    res     = shapiro.test(this.x)
    res.log = shapiro.test(trans.func(this.x))
    this.sw = rbind(data.frame(W = res$statistic, p = res$p.value), 
                    data.frame(W = res.log$statistic, p = res.log$p.value))
    rownames(this.sw) = c("raw", "trans")
    
    sw.res[[nm]] = this.sw
    rm(this.x, res, res.log, this.sw)
  }
  
  do.trans = sapply(sw.res, function(this) return(this["raw", "p"] <= sig.cut &
                                                    this["trans", "p"] > sig.cut))
  
  ret = list("shapiro.wilk.test" = sw.res, "to.transform" = names(which(do.trans == TRUE)))
}



#computes r-squared for a model fit
rsquared <- function(mod.fit, adjusted =  T) {
  yhat = fitted(mod.fit)
  y = mod.fit$model[, all.vars(lm.res$call)[1]]
  rsq = sum((yhat - mean(y))^2)/sum((y -mean(y))^2)
  
  if (adjusted == FALSE)
    return(rsq)
  
  n = nrow(mod.fit$model)
  p = length(mod.fit$coefficients) - 1 #don't count the intercept
  
  return(rsq - (1 - rsq)*(p/(n - 1 - p)))
}  



non.zero <- function(x) {
  idxs = which(x != 0)
  
  y = x[idxs]
  
  if (is.null(dim(x))) {
    names(y) = names(x)[idxs]
  } else {
    names(y) = rownames(x)[idxs]
  }
  
  return(y)
}

len.unq <- function(x) {
  return(length(unique(x)))
}


secondDerivative <- function(x) {
  #https://stackoverflow.com/questions/4471993/compute-the-elbow-for-a-curve-automatically-and-mathematically
  ret = vector(mode = "numeric", length = length(x))
  ret[c(1, length(ret))] = 0
  
  for (i in 2:length(ret))
    ret[i] = x[i+1] + x[i-1] - 2*x[i]
  
  return(ret)
}
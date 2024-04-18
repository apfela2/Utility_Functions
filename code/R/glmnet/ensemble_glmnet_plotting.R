# plotEnsembleStats
# plotEnsembleFeatureDistribution #type - one of "p-value", "beta", "beta-disc"
# plotSingleModel
# plotSingleModel_enet
# cvCurvePlot
# repeatedCVSummaryPlot


#plotting models
plotEnsembleStats <- function(ens, fname = "ensembleStats.pdf") {
  models    <- ens$models
  ens.stats <- data.frame(cv.err     = sapply(models, function(this) this$cv.stats["meanCvError", "min"]),
                          log.lambda = sapply(models, function(this) log(this$cv.stats["lambda", "min"])),
                          adj.rsq    = sapply(models, function(this) this$mod.summ$adj.r.squared),
                          rmse       = sapply(models, function(this) this$mod.summ$rmse),
                          num.beta   = sapply(models, function(this) length(this$refit.mod$coefficients) - 1))
  
  
  pdf(file = fname, height = 8.5, width = 11)  
  print(ggplot(ens.stats, aes(y = cv.err)) + 
          geom_point(aes(x = num.beta), position = "identity", size = 4.5, alpha = .3, color = "darkblue") +
          ggtitle("Number of Parameters vs Cross-Validiation Error Acrros the Ensemble") +
          ylab(expression(paste("Cross Validation Error (at best ", lambda, ")"))) +
          xlab(expression(paste("Number of Non-zero Parameters (at best ", lambda, ")"))) + 
          theme_bw(14))
  
  print(ggplot(ens.stats, aes(y = cv.err)) + 
          geom_point(aes(x = log.lambda), position = "identity", size = 4.5, alpha = .3, color = "darkblue") +
          ggtitle(expression(paste("log(", lambda, ") vs Cross Validation Error Across the Ensemble"))) +
          xlim(min(0, min(ens.stats$log.lambda)), max(ens.stats$log.lambda)) +
          ylab(expression(paste("Cross Validation Error (ar best ", lambda, ")"))) +
          xlab(expression(paste("log(", lambda, ") Picked at Each Repeat"))) +
          theme_bw(14))
  
  print(ggplot(ens.stats, aes(x = log.lambda)) + 
          geom_histogram(fill = "darkblue", bins = sqrt(nrow(ens.stats))) +
          xlab(expression(paste("log(", lambda, ") (selected at each repeated CV run)"))) +
          ggtitle(expression(paste("Distribution of Optimized log(", lambda, ") across the ensemble"))) +
          theme_bw(14))
  
  lines.df <- data.frame(value = c(mean(ens.stats$cv.err), median(ens.stats$cv.err)), 
                         Statistic = c("mean", "median"))
  
  print(ggplot(ens.stats, aes(x = cv.err)) + 
          geom_histogram(fill = "darkblue", bins = sqrt(nrow(ens.stats))) +
          xlab(expression(paste("Cross-validation Error (selected at each repeated CV run)"))) +
          ggtitle(expression(paste("Distribution of Cross-Validation Error across the ensemble"))) +
          geom_vline(data = lines.df, aes(xintercept = value, color = Statistic), size = 2, alpha = 1) +
          theme_bw(14))
  
  
  dev.off()
}

#type - one of "p-value", "beta", "beta-disc"
plotEnsembleFeatureDistribution <- function(ens, type = "p-value", freq.co = 0, fname = "ensemblePlots.pdf", 
                                            rem.null.mods = FALSE, show.rd = TRUE, feat.font.size = 8) {
  models    <- ens$models
  #ens.stats <- data.frame(cv.err     = sapply(models, function(this) this$cv.stats["meanCvError", "min"]),
  #                        log.lambda = sapply(models, function(this) log(this$cv.stats["lambda", "min"])),
  #                        adj.rsq    = sapply(models, function(this) this$mod.summ$adj.r.squared),
  #                        rmse       = sapply(models, function(this) this$mod.summ$rmse),
  #                        num.beta   = sapply(models, function(this) length(this$refit.mod$coefficients) - 1))
  
  all.coefs <- unique(unlist(lapply(models, function(this) rownames(this$mod.sum$coefficients)[-1])))
  
  ens.mat <- matrix(0, length(models), length(all.coefs))
  colnames(ens.mat) <- all.coefs
  rownames(ens.mat) <- paste0("model_", 1:nrow(ens.mat))
  
  beta.mat <- ens.mat
  
  for (i in 1:length(models)) {
    this.coefs <- models[[i]]$mod.sum$coefficients[-1, , drop = FALSE]
    if (nrow(this.coefs) == 0)
      next
    
    ens.mat[i, rownames(this.coefs)] <- -log10(this.coefs[, grep("P", colnames(this.coefs))])
    beta.mat[i, rownames(this.coefs)] <- this.coefs[, "Estimate"]
  }
  
  adj.mat <- beta.mat
  adj.mat[adj.mat < 0] <- -1
  adj.mat[adj.mat > 0] <- 1
  
  
  adj.freq <- apply(adj.mat, 2, function(this) sum(this != 0))/nrow(adj.mat)
  rm.feat.idxs  <- c(which(adj.freq < freq.co))
  rm.row.idxs   <- unique(unlist(apply(adj.mat[, rm.feat.idxs, drop = FALSE], 2, function(this) which(this != 0))))
  
  if (any(rm.feat.idxs)) {
    print(paste0("Removing ", length(rm.feat.idxs), 
                 " features because they didnt have large enough frequency, and the models they are in (", 
                 length(rm.row.idxs), ")" ))
    ens.mat  <- ens.mat[-rm.row.idxs, -rm.feat.idxs, drop = FALSE]
    beta.mat <- beta.mat[-rm.row.idxs, -rm.feat.idxs, drop = FALSE]
    adj.mat  <- adj.mat[-rm.row.idxs, -rm.feat.idxs, drop = FALSE]
  }
  
  #removing null models
  feat.cts <- apply(abs(adj.mat), 1, sum)
  null.mods <- which(feat.cts == 0)
  
  if (any(null.mods) & rem.null.mods == TRUE) {
    print(paste0("Removing ", length(null.mods), " null models"))
    ens.mat  <- ens.mat[-null.mods, , drop = FALSE]
    beta.mat <- beta.mat[-null.mods, , drop = FALSE]
    adj.mat  <- adj.mat[-null.mods, , drop = FALSE]
  }
  
  pval.sc <- sparseClust(x = ens.mat, makePlots = FALSE, less.sparse = FALSE)
  
  
  
  cols = rainbow_hcl(n = pval.sc$K)
  #palette(terrain.colors(ncols))
  
  names(cols) = names(pval.sc$centers)
  
  
  csl = list(cluster = cols)
  ha.df           = data.frame(cluster = pval.sc$classification)
  rownames(ha.df) = names(pval.sc$classification)
  
  
  
  ha = HeatmapAnnotation(df = ha.df, col = csl, which = "row", show_legend = FALSE)
  
  pdf(file = fname, height = 8.5, width = 11)  
  
  rd = hclust(dist(ens.mat), method = "ward.D2")
  rc = hclust(dist(t(ens.mat)), method = "ward.D2")
  
  #calculate max column name height based on fontsize and max number of characters in column name
  cm.per.point <- 0.0352777778 
  mnh <- ceiling(feat.font.size*cm.per.point*(max(sapply(colnames(ens.mat), nchar)) + 2))
  #the '+2' above is padding
  
  print(ha + Heatmap(ens.mat, column_names_gp = gpar(fontsize = feat.font.size), name = "-log10 p-value", show_row_names = FALSE,
                     column_title = "Heatmap of Features in the Ensemble (-log10 of p-value)", 
                     column_title_gp = gpar(fontsize = feat.font.size, fontface = "bold"),
                     cluster_columns = rc, cluster_rows = rd, 
                     show_row_dend = show.rd, column_names_max_height = unit(mnh, "cm")))
  dev.off()
  
  if (0) {
    #ToDo: prettify the plots, make row and column order same as plots above
    print(Heatmap(beta.mat, column_names_gp = gpar(fontsize = 8), name = "Beta Estimate",
                  column_title = expression(paste("Heatmap of Features in the Ensemble (Estimated beta - ", hat(beta), ")")), #))#,
                  clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2"))
    
    print(Heatmap(adj.mat, column_names_gp = gpar(fontsize = 8), name = "Heatmap of Features in the Ensemble (Discritized Beta)",
                  column_title = expression(paste(hat(beta), " in [-1, 0, 1])")), #))#,
                  clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2"))
  }
}

#plot a single model
#ToDo - extend to multiple factors, 
#     - change rsq calculation to not rely on rebuilding a model
#     - once the ensemble code is changed such that the data frame is not kept for every
#       single model, we need to change the code below to put the data matrix back 
#       into the lm object
#plotSingleModel <- function(model, this.dat, p.co = .1, plot.sig.only = FALSE, out.file = "./plots.pdf") {
plotSingleModel <- function(model, p.co = .1, plot.sig.only = FALSE, out.file = "./plots.pdf",
                            theme.size = 10) {
  this.dat <- model$refit.mod$data
  #model$refit.mod$data <- this.dat
  
  classes  <- as.character(dplyr:::summarise_each(this.dat, funs(class(.))))
  
  this.lm  <- model$refit.mod
  
  mod.summ <- model$mod.summ$coefficients
  sig.only.coefs <- mod.summ[mod.summ[, 4] <= p.co, , drop = FALSE]
  
  #mapping maes between the factors that are fitted and their names in the
  #design matrix
  fac.map <- mapFactorNames(colnames(this.dat), rownames(mod.summ)[-1], classes)
  
  rmse    <- round(model$mod.summ$rmse, 3)
  adj.rsq <- round(model$mod.summ$adj.r.squared, 3)
  
  myTheme <- ttheme_default(colhead = list(fg_params = list(parse = TRUE, cex = .8)),
                            rowhead = list(fg_params = list(parse = TRUE, cex = .5)),
                            core = list(fg_params = list(cex = .75)))
  
  
  pdf(out.file)
  myTable = xtable(round(data.frame(mod.summ, check.names = FALSE), 6))
  #g1 <- grid.table(myTable, theme = myTheme)
  grid.arrange(tableGrob(myTable, theme = myTheme), nrow = 1)
  
  myTable.short = xtable(round(data.frame(sig.only.coefs, check.names = FALSE), 6))
  if (nrow(myTable.short) > 0)
    grid.arrange(tableGrob(myTable.short, theme = myTheme), nrow = 1, newpage = TRUE)
  
  #g2 <- grid.table(myTable.short, theme = myTheme)
  #grid.arrange(g1, g2, nrow = 1, newpage = FALSE)
  
  
  resp.var = as.character(this.lm$formula)[2]
  
  lmr = fortify(this.lm)
  
  
  #ggtitle(paste0("Fitted vs Actual values of ",  resp.var, "adjusted R^2: ", 
  #               adj.rsq, "; RMSE: ", rmse)) +
  
  plot.title = paste0("Fitted vs Actual values of ",  resp.var)
  plot.subtitle = paste0("adjusted R^2: ", adj.rsq, "; RMSE: ", rmse)
  
  print(ggplot(lmr, aes_string(x = ".fitted", y = resp.var)) + 
          geom_point(size = 4) +
          xlab(paste0("Fitted values of ", resp.var)) +
          ylab(paste0("Actual values of ", resp.var)) +
          ggtitle(bquote(atop(.(plot.title), atop(italic(.(plot.subtitle)), "")))) +
          theme_bw(theme.size)) 
  
  
  
  gghist(this.dat[, resp.var], vname = resp.var, nbins = 30, larger.font = TRUE)
  
  ylab <- paste0("Fitted values of ", resp.var)
  
  
  #plot each individual varialbe
  pvs <- names(lmr)
  pvs <- pvs[-c(which(pvs == resp.var), grep("\\.", pvs))]
  
  if (plot.sig.only == TRUE & nrow(sig.only.coefs) >= 1) {
    sig.only.names <- rownames(sig.only.coefs)
    
    #convert name in summary table for a factor back to the factor name
    match.names <- match(sig.only.names, names(fac.map))
    if (any(!is.na(match.names))) {
      map.idxs <- which(!is.na(match.names))
      sig.only.names[map.idxs] <- fac.map[match.names[map.idxs]]
    }
    pvs <- intersect(pvs, sig.only.names)
  }
  
  pv.classes <- summarise_each(this.dat[, pvs, drop = FALSE], funs(class(.)))
  
  for (v in pvs) {
    visreg(this.lm, xvar = v, points = list(cex = 2.5), ylab = ylab) 
  }
  
  for (v in pvs) {
    if (class(lmr[, v]) == "factor") {
      print(ggplot(lmr, aes_string(x = v, y = resp.var)) + 
              geom_boxplot(size = 4) + 
              #facet_grid(as.formula(paste(". ~", fac.var))) +
              #geom_smooth(method = "lm", se = TRUE)+
              theme_bw(theme.size))
    } else {
      print(ggplot(lmr, aes_string(x = v, y = resp.var)) + 
              geom_point(size = 4) + 
              #facet_grid(as.formula(paste(". ~", fac.var))) +
              geom_smooth(method = "lm", se = TRUE)+
              theme_bw(theme.size))
    }
  }
  
  if (any(pv.classes == "factor")) {
    facs <- names(pv.classes)[pv.classes == "factor"]
    
    pvs <- pvs[-which(pvs %in% facs)]
    
    for (f in facs) {
      for (v in pvs) {
        visreg(this.lm, xvar = v, by = f, points = list(cex = 2.5), ylab = ylab) 
      }
    }
    
    for (f in facs) {
      for (v in pvs) {
        print(ggplot(lmr, aes_string(x = v, y = resp.var)) + 
                geom_point(size = 4) + 
                facet_grid(as.formula(paste(". ~", f))) +
                geom_smooth(method = "lm", se = TRUE)+
                theme_bw(theme.size))
      }
    }
  }
  
  print(autoplot(this.lm, which = 1:6, colour = 'dodgerblue3',
                 smooth.colour = 'black', smooth.linetype = 'dashed',
                 ad.colour = 'blue',
                 label.size = 3, label.n = 5, label.colour = 'blue',
                 ncol = 3))
  
  dev.off()
}


plotSingleModel_enet <- function(model, enet.mod, dset, out.file = "./bestMod_enet_plots.pdf",
                                 theme.size = 10, sparse.coef = FALSE) {
  this.dat <- model$refit.mod$data
  #model$refit.mod$data <- this.dat
  resp.var <- as.character(model$refit.mod$formula)[2]
  
  if (class(dset) == "data.frame")
    dset <- makeModMatrix(dset)
  
  colnames(dset) <- make.names(colnames(dset))
  
  best.s <- ifelse(sparse.coef == TRUE, model$cv.stats["lambda", "sparse"], model$cv.stats["lambda", "min"])
  
  betas <- coef(enet.mod, s = best.s)
  betas <- betas[betas[, 1] != 0, ]
  betas <- data.frame(Estimate = betas)
  rownames(betas) <- make.names(rownames(betas))
  
  y.hat <- predict(enet.mod, best.s, newx = dset)
  y     <- this.dat[, resp.var]
  
  #ToDO: make the two lines below not ugly
  ym  <- data.frame(y); colnames(ym) <- resp.var
  yhm <- data.frame(y.hat); colnames(yhm) <- ".fitted"
  
  
  #calculate rmse, r-squared and adjusted r-squared
  mod.stats <- getModRefitErrors_vec(y, y.hat, nrow(betas) - 1)
  
  rmse    <- round(mod.stats$rmse, 3)
  adj.rsq <- round(mod.stats$adj.r.squared, 3)
  
  #classes  <- as.character(dplyr:::summarise_each(this.dat, funs(class(.))))
  
  #this.lm  <- model$refit.mod
  
  #mod.summ <- model$mod.summ$coefficients
  #sig.only.coefs <- mod.summ[mod.summ[, 4] <= p.co, , drop = FALSE]
  
  #mapping maes between the factors that are fitted and their names in the
  #design matrix
  
  myTheme <- ttheme_default(colhead = list(fg_params = list(parse = TRUE, cex = .8)),
                            rowhead = list(fg_params = list(parse = TRUE, cex = .6)),
                            core = list(fg_params = list(cex = .75)))
  
  
  pdf(out.file)
  myTable = xtable(round(data.frame(betas, check.names = FALSE), 6))
  #g1 <- grid.table(myTable, theme = myTheme)
  grid.arrange(tableGrob(myTable, theme = myTheme), nrow = 1)
  
  
  
  
  #ggtitle(paste0("Fitted vs Actual values of ",  resp.var, "adjusted R^2: ", 
  #               adj.rsq, "; RMSE: ", rmse)) +
  lmr <- cbind(yhm, dset[, rownames(betas)[-1]], ym)
  
  plot.title = paste0("Fitted vs Actual values of ",  resp.var)
  plot.subtitle = paste0("adjusted R^2: ", adj.rsq, "; RMSE: ", rmse)
  
  print(ggplot(lmr, aes_string(x = ".fitted", y = resp.var)) + 
          geom_point(size = 4) +
          xlab(paste0("Fitted values of ", resp.var)) +
          ylab(paste0("Actual values of ", resp.var)) +
          ggtitle(bquote(atop(.(plot.title), atop(italic(.(plot.subtitle)), "")))) +
          theme_bw(theme.size)) 
  
  
  
  gghist(lmr[, resp.var], vname = resp.var, nbins = 30, larger.font = TRUE)
  
  ylab <- paste0("Fitted values of ", resp.var)
  
  
  #plot each individual varialbe
  pvs <- names(dplyr:::select(lmr, -one_of(".fitted", resp.var)))
  pv.classes <- sapply(lmr[, pvs], class)
  
  for (v in pvs) {
    if (class(lmr[, v]) == "factor" | class(lmr[, v]) == "character") {
      print(ggplot(lmr, aes_string(x = v, y = resp.var)) + 
              geom_boxplot(size = 4) + 
              #facet_grid(as.formula(paste(". ~", fac.var))) +
              #geom_smooth(method = "lm", se = TRUE)+
              theme_bw(theme.size))
    } else {
      print(ggplot(lmr, aes_string(x = v, y = resp.var)) + 
              geom_point(size = 4) + 
              #facet_grid(as.formula(paste(". ~", fac.var))) +
              geom_smooth(method = "lm", se = TRUE)+
              theme_bw(theme.size))
    }
  }
  
  if (any(pv.classes == "factor")) {
    facs <- names(pv.classes)[pv.classes == "factor" | pv.classes == "character"]
    
    pvs <- pvs[-which(pvs %in% facs)]
    
    for (f in facs) {
      for (v in pvs) {
        print(ggplot(lmr, aes_string(x = v, y = resp.var)) + 
                geom_point(size = 4) + 
                facet_grid(as.formula(paste(". ~", f))) +
                geom_smooth(method = "lm", se = TRUE)+
                theme_bw(theme.size))
      }
    }
  }
  
  dev.off()
}

cvCurvePlot <- function(cv.stats, lambdas, nzeros, err, err.sd, err.name, plot.title, f.size = 16) {
  p <- list()
  anno <- cv.stats %>% 
    t() %>%
    data.frame(stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column(var = "criterion") %>%
    mutate(lambda = log(lambda))
  
  d <- data.frame('lambda' = log(lambdas), 'cvm' = err,
                  'cvup' = err + err.sd, 
                  'cvlo' = err - err.sd,
                  'nz' = nzeros)
  
  xlab <- 'log(Lambda)'
  ylab <- err.name
  
  plot.data        <- d
  plot.data$label  <- rep(max(d$cvup), nrow(plot.data))
  
  
  indexer <- seq(1, nrow(plot.data), length = 12)
  
  p[[1]] <- ggplot2::ggplot(plot.data) +
    geom_point(aes_string('lambda', 'cvm')) +
    geom_errorbar(aes_string(x = 'lambda', ymin = 'cvlo', ymax = 'cvup')) +
    geom_text(data = plot.data[indexer, ], aes(x = lambda, y = label, label = nz)) +
    geom_point(data = anno, 
               mapping = aes(x = lambda,  y = meanCvError, col = criterion, shape = criterion),
               size = 4, alpha = .75) +
    geom_vline(data = anno, aes(xintercept = lambda), linetype = "dashed") +
    xlab('log(Lambda)') +
    ylab(ylab) +
    labs(title = plot.title) +
    theme_bw(f.size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

repeatedCVSummaryPlot <- function(full.model, models, f.size = 16, family) {
  all.lambda.types <- colnames(models[[1]]$cv.stats)
  
  #making data frame that captures errors across each repeat run
  errDist <- lapply(all.lambda.types, function(et) {
    sapply(models, function(x) x$cv.stats["meanCvError", et])
  })
  names(errDist) <- all.lambda.types
  
  err.plot <- data.frame(errDist, stringsAsFactors = FALSE) %>% 
    gather(key = "Criterion", value = "meanCvError")
  
  #making data frame that captures number of non-zero coefficients (including intercept)
  #across each repeated run
  nzDist <- lapply(all.lambda.types, function(et) {
    sapply(models, function(x) x$cv.stats["numNonZero", et])
  })
  names(nzDist) <- all.lambda.types
  nz.plot       <- data.frame(nzDist, stringsAsFactors = ) %>% 
    gather(key = "Criterion", value = "nnz") 
  
  
  p <- NULL
  
  #plotting doesn't work multinomial...not sure why, not time to debug autoplot
  #going to ignore for now
  # if (family == "multinomial") {
  #   p[[1]] <- NA
  # } else {
  #   p[[1]] <- autoplot(full.model, label = FALSE, main = "Elastic Net Selection Paths") + 
  #               guides(colour = FALSE) +
  #               theme_bw(f.size) +
  #               theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  # }
  # 
  p[[1]] <- ggplot(err.plot) +
    geom_density(aes(x = meanCvError, y = ..scaled.., fill = Criterion), 
                 alpha = .5) +
    xlab(models[[1]]$err.name) + 
    ylab("Scaled Density") + 
    labs(title = paste0("Mean CV-Error Over ", length(models), " Repeated CV Runs"),
         subtitle = "Distribution for Each Criterion") + 
    theme_bw(f.size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = .5))
  
  p[[2]] <- ggplot(nz.plot) +
    geom_density(aes(x = nnz, y = ..scaled.., fill = Criterion), 
                 alpha = .5) +
    xlab("Number of Non-zero Coefficients") + 
    ylab("Scaled Density") + 
    labs(title = paste0("Number of Non-zero Coefficients Across ", length(models), " Repeated CV Runs"),
         subtitle = "Distribution for Each Criterion") + 
    theme_bw(f.size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = .5))
  return(p)
}
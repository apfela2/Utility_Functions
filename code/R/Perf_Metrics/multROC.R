# Make ROC curves from multiple models on same plot

# Needs "submodules/P01229_IO_TMB_Surrogates/code/R/eval_models.R"

# preds <- all_preds
# response <- PARAMS$response
# data <- PARAMS[["dat"]]


multROC <- function(data, preds, response, plot.title) {
  plot.dat <- vector("list", length(preds))
  for (i in 1:length(preds)) {
    this.preds <- preds[[i]]$pred.mat
    
    # if (PARAMS$recalibrate)
    #   this.preds <- apply(this.preds, 2, function(this) recal(this.y = dat.all$MutLoad, this.y_hat = this))
    
    this.preds <- data.frame(this.preds)
    # agg_y_hat  <- apply(this.preds, 1, mean)
    # agg_y_hat  <- apply(this.preds, 1, median)
    
    num_y <- data[[i]][[response]]
    levels(num_y) <- (1:length(levels(num_y))) 
    num_y <- as.numeric(num_y) - 1   # Note: 1 is CR/PR and 0 is PD/SD/NE
    cutoffs <- unique(unlist(round(this.preds, 3)))
    # cutoffs <- round(this.preds, 3)
    cutoffs <- cutoffs[order(cutoffs)]
    
    compBinPerf <- function(actual, pred) {
      tp = sum(pred == 1 & actual == 1) #true positive
      fp = sum(pred == 1 & actual == 0) #false positive
      tn = sum(pred == 0 & actual == 0) #true negative
      fn = sum(pred == 0 & actual == 1) #false negative
      
      ppv <- tp/(tp + fp)
      npv <- tn/(tn + fn)
      sens <- tp/(tp + fn)
      spec <- tn/(tn + fp)
      return(c(tp = tp, fp = fp, tn = tn, fn = fn, ppv = ppv, npv = npv, sens = sens, spec = spec))
    }
    
    
    my.pred <- vector("list", ncol(this.preds))
    for (j in 1:ncol(this.preds)) {
      my.pred[[j]] <- data.frame(do.call("rbind", lapply(cutoffs,
                                                         function(co) {
                                                           
                                                           return(c("Cutoff" = co, compBinPerf(actual = num_y, pred = ifelse(this.preds[,j] >= co, 1, 0)), rep = j))
                                                         })), stringsAsFactors = FALSE)
      
      # my.pred[[j]] <- data.frame(do.call("rbind", return(c("Cutoff" = cutoffs, compBinPerf(actual = num_y, pred = ifelse(this.preds[,j] >= cutoffs, 1, 0)), rep = j))
      # ), stringsAsFactors = FALSE)
      
    }
    
    my.pred <- do.call("rbind", my.pred)  
    
    my.pred_g <- group_by(my.pred, Cutoff)
    med_pred <- summarise(my.pred_g, ppv = median(ppv), npv = median(npv), sens = median(sens), spec = median(spec))
    
    bin.perf  <- lapply(this.preds, function(this) evalBinPred(y = num_y, y_hat = this))
    all.aucs  <- data.frame(auc = sapply(bin.perf, function(x) x$auc))
    #m.auc     <- median(all.aucs$auc)
    
    # bin.stats <- do.call("rbind", lapply(seq_len(length(bin.perf)), 
    #                                      function(idx){
    #                                        this.dset      <- bin.perf[[idx]]$pred.mat
    #                                        this.dset$rep  <- rep(idx, nrow(this.dset))
    #                                        this.dset$rank <- seq_len(nrow(this.dset))
    #                                        #this.dset$Prevalence <- sapply(as.vector(this.dset$Cutoff), function(x) (sum(this.dset$Cutoff >= x) - 1)/(nrow(this.dset) - 1))
    #                                        #this -1's in the line above account for the Inf in Cutoff
    #                                        return(this.dset)
    #                                      }
    # ))
    
    #make plots
    outer.loop.str <- paste0(PARAMS$outer_cv_nfolds, "-fold CV, repeated ", PARAMS$num_outer_rep, " times")
    
    title <- paste0(plot.title, "\n(", outer.loop.str, ")")
    plot.dat[[i]] <- mutate(med_pred, model = deparse(preds[[i]]$formula), AUC = median(all.aucs$auc))
    # plot.dat[[i]] <- plot.dat[[i]][plot.dat[[i]]$rep == which(near(all.aucs, median(all.aucs$auc)))[1],]
  }
  plot.dat <- do.call("rbind", plot.dat)
  
  model_levels <- rep(NA, length(preds))
  for (i in 1:length(preds)) {
    model_levels[i] <- deparse(preds[[i]]$formula)
  }
  plot.dat$model <- factor(plot.dat$model, levels = model_levels )
  
  # plots[[length(plots) + 1]] <-
  perf <- list()
  
  # Make ROC plot with curve for each model
  perf[["ROC_plot"]] <- ggplot(plot.dat, aes(y = sens, x = 1 - spec, color = model)) + 
    # geom_line(aes(group = rep), alpha = .1, size = 1) +
    # geom_line() + 
    geom_point(aes(shape=model), size = 2) +
    scale_shape_manual(values = seq(65, length.out = length(unique(plot.dat$model)))) +
    guides(shape = guide_legend(override.aes = list(size=5))) +
    # scale_x_reverse() +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), size = .8, col = "grey") +
    ylab("Sensitivity") +
    xlab("1 - Specificity") +
    labs(title = title) +
    coord_fixed() +
    theme_bw(PARAMS$font_size) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = rel(.9)),
          axis.title.x = element_text(size = rel(.9), face = "bold"),
          axis.title.y = element_text(size = rel(.9), face = "bold"))
  
  # Make PPV/NPV plot
  perf[["PPV_plot"]] <- ggplot(plot.dat, aes(y = ppv, x = 1 - npv, color = model)) + 
    # geom_line(aes(group = rep), alpha = .1, size = 1) +
    # geom_line() +
    geom_point(aes(shape=model), size = 2) +
    scale_shape_manual(values = seq(65, length.out = length(unique(plot.dat$model)))) +
    guides(shape = guide_legend(override.aes = list(size=5))) +
    # scale_x_reverse() +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), size = .8, col = "grey") +
    ylab("PPV") +
    xlab("1 - NPV") +
    labs(title = title) +
    coord_fixed() +
    theme_bw(PARAMS$font_size) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = rel(.9)),
          axis.title.x = element_text(size = rel(.9), face = "bold"),
          axis.title.y = element_text(size = rel(.9), face = "bold"))
  
  # Calculate median AUC for each model
  plot.dat_g <- group_by(plot.dat, model)
  perf[["AUC"]] <- summarise(plot.dat_g, AUC = round(mean(AUC), 2))
  
  return(perf)
}


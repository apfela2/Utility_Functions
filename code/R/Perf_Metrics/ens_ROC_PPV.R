# Depends on evalBinPred

## Abraham Apfel
## December 29, 2018


# Make ROC curve from ens prediction matrix


ens_ROC_PPV <- function(this.preds, y, subtitle = "Repeated CV", title, outer_folds,
                        outer_rep, limx = NULL, labx = "Gene Score", font = 18, ens = T){
  
  
  this.preds <- data.frame(this.preds)
  # agg_y_hat  <- apply(this.preds, 1, mean)
  agg_y_hat  <- apply(this.preds, 1, median)
  
  
  bin.perf  <- lapply(this.preds, function(this) evalBinPred(y = y, y_hat = this))
  all.aucs  <- data.frame(auc = sapply(bin.perf, function(x) x$auc))
  #m.auc     <- median(all.aucs$auc)
  
  bin.stats <- do.call("rbind", lapply(seq_len(length(bin.perf)), 
                                       function(idx){
                                         this.dset      <- bin.perf[[idx]]$pred.mat
                                         this.dset$rep  <- rep(idx, nrow(this.dset))
                                         this.dset$rank <- seq_len(nrow(this.dset))
                                         #this.dset$Prevalence <- sapply(as.vector(this.dset$Cutoff), function(x) (sum(this.dset$Cutoff >= x) - 1)/(nrow(this.dset) - 1))
                                         #this -1's in the line above account for the Inf in Cutoff
                                         return(this.dset)
                                       }
  ))
  #adding column for prevelance
  
  #make plots
  plots <- list()
  
  # plot.subtitle <- paste0("model type: ", mod.type, ", feature set: ", feature.set)
  
  # if (!is.na(sel.alpha))
  #   plot.subtitle <- paste0(plot.subtitle, " (", num.genes.full.mod, " selected), alpha: ", sel.alpha)
  
  outer.loop.str <- paste0(outer_folds, "-fold CV, repeated ", outer_rep, " times")
  
  #ToDO:
  #1) AUC CURVE ACROSS ALL REPEATS
  #2) SENS AND SPEC VS CUTOFF ACROSS ALL REPEATS
  #3) NPV, PVV, ACC VS CUTOFF ACROSS ALL REPEATS
  
  if(ens){
    plot.title    <- paste0(title, "\n(", outer.loop.str, ")\nmedian AUC = ", 
                            round(median(all.aucs$auc), 2), ", sd = ", round(sd(all.aucs$auc), 2))
  } else {
    plot.title    <- paste0(title, "\n AUC = ", 
                            round(median(all.aucs$auc), 2))
  }
  
  plot.dat <- bin.stats %>% dplyr::select(True.positive.rate, Specificity, Cutoff, rank, rep)
  
  if(ens){
    plots[["ROC"]] <- 
      ggplot(plot.dat, aes(y = True.positive.rate, x = 1 - Specificity)) + 
      geom_line(aes(group = rep), alpha = .1, size = 1) +
      geom_smooth(col = "black") + 
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), size = .8, col = "grey") +
      ylab("Sensitivity") +
      xlab("1 - Specificity") +
      labs(title = plot.title, subtitle = subtitle) +
      coord_fixed() +
      theme_bw(font) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, size = rel(.9)),
            axis.title.x = element_text(size = rel(.9), face = "bold"),
            axis.title.y = element_text(size = rel(.9), face = "bold"))
  } else {
    plots[["ROC"]] <- 
      ggplot(plot.dat, aes(y = True.positive.rate, x = 1 - Specificity)) + 
      geom_line(aes(group = rep), size = 1) +
      # geom_smooth(col = "black") + 
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), size = 0.8, col = "grey") +
      ylab("Sensitivity") +
      xlab("1 - Specificity") +
      labs(title = plot.title, subtitle = subtitle) +
      coord_fixed() +
      theme_bw(font) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, size = rel(.9)),
            axis.title.x = element_text(size = rel(.9), face = "bold"),
            axis.title.y = element_text(size = rel(.9), face = "bold"))
    
  }
  
  rm(plot.title)
  
  
  #AGGREGATING DATA TO MAKE PLOTS
  comb.dat <- bin.stats %>%
    dplyr::select(-one_of("Cutoff", "rep", "rank")) %>%
    tidyr::gather(key = "Metric", value = "Performance")
  #plot.dat$Cutoff <- rep(bin.stats$Cutoff, times = length(unique(plot.dat$Metric)))
  
  comb.dat <- cbind(comb.dat,
                    do.call("rbind", lapply(seq_len(length(unique(comb.dat$Metric))), 
                                            function(i) return(bin.stats[, c("rep", "rank", "Cutoff")]))))
  q975 <- function(x) { return(quantile(x, .975, na.rm = T))}
  q025 <- function(x) { return(quantile(x, .025, na.rm = T))}
  mn <- function(x) {return(mean(x, na.rm = T))}
  med <- function(x) {return(median(x, na.rm = T))}
  
  
  
  agg.dat <- comb.dat %>%
    dplyr::select(-rep) %>%
    group_by(Metric, rank) %>%
    summarise_all(funs(mn, med, q975, q025)) %>%
    data.frame(stringsAsFactors = FALSE) %>%
    dplyr::mutate(Metric = ifelse(Metric == "True.positive.rate", "Sensitivity", Metric))
  agg.dat$Metric <- gsub("\\.", " ", agg.dat$Metric)
  
  x.min.breaks <-  seq(round(min(bin.stats$Cutoff, na.rm = T), 1),
                       round(max(bin.stats$Cutoff[bin.stats$Cutoff < Inf], na.rm = T), 1), 0.1)
  y.min.breaks <- seq(0, 1, 0.05)
  
  # 
  #sens and spec vs cutoff
  # plot.title <- paste0("Performance vs Predicted OR Cutoff (", outer.loop.str, ")")
  # plot.dat   <- agg.dat %>%
  #   dplyr::filter(Metric %in% c("Sensitivity", "Specificity", "Prevalence"))
  # 
  # # 
  # plots[[length(plots) + 1]] <-  
  #   ggplot(plot.dat, aes(y = Performance_med, x = Cutoff_med, col = Metric, fill = Metric)) + 
  #   geom_ribbon(aes(ymin = Performance_q025, ymax = Performance_q975), alpha = .5, colour=NA) +
  #   geom_line(size = 1, alpha = .5) +
  #   ylab("Performance") +
  #   xlab("Gene Score") +
  #   #geom_smooth() + 
  #   labs(title = plot.title, subtitle = plot.subtitle) +
  #   scale_x_continuous(minor_breaks = x.min.breaks) +
  #   scale_y_continuous(minor_breaks = y.min.breaks) +
  #   theme_bw(PARAMS$font_size) +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         plot.subtitle = element_text(hjust = 0.5, size = rel(.9)))
  
  rm(plot.dat) 
  
  
  ##ppv, npv, acc vs cutoff
  
  
  plot.dat <- agg.dat %>%
    dplyr::filter(Metric %in% c("Accuracy", "Positive predictive value", "Negative predictive value", "Prevalence"))
  
  if (ens) {
    plot.title <- paste0(title, "\n(", outer.loop.str, ")")
  } else {
    plot.title <- paste0(title)
  }
  plots[["PPV"]] <- 
    ggplot(plot.dat, aes(y = Performance_med, x = Cutoff_med, col = Metric, fill = Metric)) + 
    geom_ribbon(aes(ymin = Performance_q025, ymax = Performance_q975), alpha = .5, colour=NA) +
    geom_line(size = 1) +
    scale_x_continuous(minor_breaks = x.min.breaks) +
    scale_y_continuous(minor_breaks = y.min.breaks) +
    ylab("Performance") +
    xlab(labx) +
    #stat_smooth(se = FALSE, level = .5) + 
    labs(title = plot.title, subtitle = subtitle) +
    theme_bw(font) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = rel(.9)))
  
  if(!is.null(limx)){
    plots[["PPV"]] <- plots[["PPV"]] + xlim(limx)
  }
  
  # cv.perf.list[[mn]] <- list(plots = plots, binary.perf = bin.perf, continuous.perf = cont.perf, aucs = all.aucs,
  #                            featureSet = feature.set, modelType = mod.type, alpha = sel.alpha, 
  #                            agg.y.hat = agg_y_hat, agg.dat = agg.dat, bin.stats = bin.stats, num.genes.full.mod = num.genes.full.mod)
  
  cv.perf <- list(plots = plots, binary.perf = bin.perf, aucs = all.aucs,
                  agg.y.hat = agg_y_hat, agg.dat = agg.dat, bin.stats = bin.stats)
  
  return(cv.perf)
}

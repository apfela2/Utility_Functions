## Make plot of AUC vs time


# predict.time is vector of times you want to plot ROC curves for
# Stime is survival time
# status is censoring indicator where 1 means event and 0 is censored
# marker is your matrix of predicted biomarker scores
# data is your data object


# Stime <- PARAMS$response
# status <- PARAMS$status
# marker <- all_preds
# data <- PARAMS$dat
# plot.title <- "TEST"
# predict.time = seq(30.4375, 365.25, 30.4375)

multTimeAUC <- function(predict.time = seq(30.4375, 365.25, 30.4375), Stime,
                        status, marker, data, plot.title) {
  
  plot.dat <- vector("list", length(marker))
  
  for (i in 1:length(marker)) {
    this.preds <- marker[[i]]$pred.mat
    
    survivalROC_data <- list()
    
    for (j in 1:ncol(this.preds)) {
      survivalROC_helper <- function(t) {
        survivalROC(Stime        = data[[i]][[Stime]],
                    status       = data[[i]][[status]],
                    marker       = this.preds[,j],
                    predict.time = t,
                    method       = "NNE",
                    cut.values = unique(as.vector(round(this.preds, 2))),
                    span = 0.25 * nrow(data[[i]])^(-0.20))
      }
      
      
      ## Evaluate every Month
      survivalROC_data[[j]] <- data_frame(t = predict.time) %>%
        mutate(survivalROC = purrr::map(t, survivalROC_helper),
               ## Extract scalar AUC
               auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
               rep = j,
               ## Put cut off dependent values in a data_frame
               df_survivalROC = purrr::map(survivalROC, function(obj) {
                 as_data_frame(obj[c("cut.values","TP","FP")])
               })) %>%
        dplyr::select(-survivalROC) %>%
        unnest() %>%
        arrange(t, FP, TP)
    }
    
    # Collapse list of predictions into dataframe, one row for each repeat
    survivalROC_data <- do.call("rbind", survivalROC_data)
    
    # Calculate median TPR/FPR/AUC across all repeats
    my.pred_g <- group_by(survivalROC_data, t, cut.values)
    med_pred <- summarise(my.pred_g, TPR = median(TP, na.rm = T), FPR = median(FP, na.rm = T), AUC = median(auc, na.rm = T))
    med_pred <- ungroup(med_pred)
    
    
    plot.dat[[i]] <- mutate(med_pred, model = deparse(marker[[i]]$formula))
    
    outer.loop.str <- paste0(PARAMS$outer_cv_nfolds, "-fold CV, repeated ", PARAMS$num_outer_rep, " times")
    
    title <- paste0(plot.title, "\n(", outer.loop.str, ")")
  }
  
  plot.dat <- do.call("rbind", plot.dat)
  
  model_levels <- rep(NA, length(marker))
  for (i in 1:length(marker)) {
    model_levels[i] <- deparse(marker[[i]]$formula)
  }
  plot.dat$model <- factor(plot.dat$model, levels = model_levels )
  
  # Create list to store plots
  perf <- list()  
  
  # Plot
  
  perf[["timeAUC_plot"]] <- ggplot(plot.dat, mapping = aes(x = t, y = AUC, color = model)) +
    geom_line(show.legend = F) +
    geom_point(aes(shape=model), size = rel(3)) +
    scale_shape_manual(values = seq(65, length.out = length(unique(plot.dat$model)))) +
    guides(shape = guide_legend(override.aes = list(size=5))) +
    theme_bw(18) +
    theme(axis.text.x = element_text(vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank()) +
    xlab("Months") +
    ylab("Time-dependent AUC") +
    labs(title = plot.title) +
    scale_x_continuous(breaks = seq(30.4375, max(predict.time), 30.4375),
                       labels = seq(1, max(predict.time)/30.4375)) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    ylim(0,1)
  
  return(perf[["timeAUC_plot"]])
}

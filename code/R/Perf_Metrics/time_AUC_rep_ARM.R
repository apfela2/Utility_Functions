## Make plot of AUC vs time for matrix of predictions


time_AUC_rep <- function(predict.time = seq(30.4375, 730.5, 30.4375), Stime, status, marker, data, plot.title = paste("Time vs AUC Plot")) {
 
  
  
  multTimeAUC_arm <- function(predict.time = seq(30.4375, 365.25, 30.4375), Stime,
                              status, arm, marker, data, plot.title) {
    
    nivo.dat <- vector("list", length(marker))
    combo.dat <- vector("list", length(marker))
    plot.dat <- vector("list", length(marker))
    
    for (i in 1:length(marker)) {
      this.preds <- marker[[i]]$pred.mat
      this.preds <- marker[[i]]$pred.mat
      this.preds <- as.data.frame(this.preds)
      this.preds$Stime <- data[[i]][, Stime]
      this.preds$status <- data[[i]][, status]
      this.preds$ARM <- data[[i]][, arm]
      this.preds$ARM <- with(this.preds, ifelse(ARM == "Nivolumab 3 mg/kg + Ipilimumab 1 mg/kg", "Combo",
                                                ifelse(ARM == "Nivolumab 3 mg/kg + Ipilimumab Placebo", "Nivo", NA)))
      
      combo <- filter(this.preds, ARM == "Combo")
      combo <- dplyr::select(combo, -Stime, -status, -ARM)
      combo <- as.matrix(combo)
      
      survivalROC_data <- list()
      
      for (j in 1:ncol(combo)) {
        survivalROC_helper <- function(t) {
          survivalROC(Stime        = this.preds[which(this.preds$ARM == "Combo"),"Stime"],
                      status       = this.preds[which(this.preds$ARM == "Combo"),"status"],
                      marker       = combo[,j],
                      predict.time = t,
                      method       = "NNE",
                      cut.values = unique(as.vector(round(combo, 2))),
                      span = 0.25 * nrow(combo)^(-0.20))
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
      
      
      combo.dat[[i]] <- mutate(med_pred, model = deparse(marker[[i]]$formula), ARM = "Combo")
      
      outer.loop.str <- paste0(PARAMS$outer_cv_nfolds, "-fold CV, repeated ", PARAMS$num_outer_rep, " times")
      
      title <- paste0(plot.title, "\n(", outer.loop.str, ")")
      
      ##################################################################
      
      nivo <- filter(this.preds, ARM == "Nivo")
      nivo <- dplyr::select(nivo, -Stime, -status, -ARM)
      nivo <- as.matrix(nivo)
      
      survivalROC_data <- list()
      
      for (j in 1:ncol(nivo)) {
        survivalROC_helper <- function(t) {
          survivalROC(Stime        = this.preds[which(this.preds$ARM == "Nivo"),"Stime"],
                      status       = this.preds[which(this.preds$ARM == "Nivo"),"status"],
                      marker       = nivo[,j],
                      predict.time = t,
                      method       = "NNE",
                      cut.values = unique(as.vector(round(nivo, 2))),
                      span = 0.25 * nrow(nivo)^(-0.20))
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
      
      
      nivo.dat[[i]] <- mutate(med_pred, model = deparse(marker[[i]]$formula), ARM = "Nivo")
      
      plot.dat[[i]] <- rbind(combo.dat[[i]], nivo.dat[[i]])
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
      scale_shape_manual(values = seq(49, length.out = length(unique(plot.dat$model)))) +
      guides(shape = guide_legend(override.aes = list(size=5))) +
      facet_wrap(~ ARM) +
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
  
  
   
  survivalROC_data <- list()
  
  for (i in 1:ncol(marker)) {
    survivalROC_helper <- function(t) {
      survivalROC(Stime        = data[[Stime]],
                  status       = data[[status]],
                  marker       = marker[,i],
                  predict.time = t,
                  method       = "NNE",
                  span = 0.25 * nrow(data)^(-0.20))
    }
    
    
    ## Evaluate at each
    survivalROC_data[[i]] <- tibble(t = predict.time) %>%
      mutate(survivalROC = purrr::map(t, survivalROC_helper),
             ## Extract scalar AUC
             auc = purrr::map_dbl(survivalROC, magrittr::extract2, "AUC"),
             ## Put cut off dependent values in a tibble
             df_survivalROC = purrr::map(survivalROC, function(obj) {
               as_tibble(obj[c("cut.values","TP","FP")])
             })) %>%
      dplyr::select(-survivalROC) %>%
      unnest() %>%
      arrange(t, FP, TP)
  }
  
  
  
  ## Change labels to months
  # time_labels <- c('91.3' = "3 Months", '182.6' = "6 Months", '365.25' = "12 Months",
  # '547.8' = "18 Months", '730.5' = "24 Months")
  
  
  ## Plot
  survivalROC_data_long <- bind_rows(survivalROC_data, .id = "rep")
  
  survivalROC_data_long_g <- group_by(survivalROC_data_long, t)
  med_auc <- summarise(survivalROC_data_long_g, med_auc = median(auc))
  
  
  ggplot(data = survivalROC_data_long, mapping = aes(x = t, y = auc)) +
    geom_point(aes(group = rep), alpha = .1, size = .3, col = "grey") +
    geom_line(aes(group = rep), alpha = .3, size = .3, col = "grey") +
    geom_line(data = med_auc, mapping = aes(x = t, y = med_auc), col = "black") + 
    theme_bw(PARAMS$font_size) +
    theme(axis.text.x = element_text(vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank()) +
    xlab("Landmark Time (Months)") +
    ylab("Time-dependent AUC") +
    labs(title = plot.title) +
    scale_x_continuous(breaks = predict.time,
                       labels = seq(1, length(predict.time))) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    ylim(0,1)
}

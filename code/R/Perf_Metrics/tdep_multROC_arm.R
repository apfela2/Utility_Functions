# Make time-dependent ROC curves for multiple models by treatment


# predict.time is vector of times you want to plot ROC curves for
# Stime is survival time
# status is censoring indicator where 1 means event and 0 is censored
# marker is your matrix of predicted biomarker scores
# data is your data object

# Stime <- PARAMS$response
# status <- PARAMS$status
# marker <- all_preds
# data <- PARAMS$dat
# arm <- "ACTARM"
# plot.title <- "TEST"
# predict.time = c(182.6, 273.9, 365.25)
# time_labels = c('182.6' = "6 Months", '273.9' = "9 Months", '365.25' = "12 Months", 'Combo' = "Combo", "Nivo" = "Nivo")

tdep_multROC_arm <- function(predict.time = c(182.6, 273.9, 365.25),
                             time_labels = c('182.6' = "6 Months",'273.9' = "9 Months",
                                             '365.25' = "12 Months"),
                             Stime, status, arm, marker, data, plot.title) {
  
  nivo.dat <- vector("list", length(marker))
  combo.dat <- vector("list", length(marker))
  plot.dat <- vector("list", length(marker))
  
  for (i in 1:length(marker)) {
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
      
      ## Evaluate at each time point
      survivalROC_data[[j]] <- data_frame(t = predict.time)%>%
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
    
    # Calculate C-index
    
    all.cindex <- sapply(as.data.frame(combo), function (x) {
      calc <- concordance.index (x, surv.time = this.preds[which(this.preds$ARM == "Combo"), "Stime"],
                                 surv.event = this.preds[which(this.preds$ARM == "Combo"),"status"])
      return(calc$c.index)
    })
    
    combo.dat[[i]] <- mutate(med_pred, model = deparse(marker[[i]]$formula), Cindex = median(all.cindex), ARM = "Combo")
    
    ## Plot
    
    # plot.subtitle <- paste("Mean C-index =", round(mean(cindex), 2), "Std =", round(sd(cindex), 2), "across", ncol(marker), "repeats")
    outer.loop.str <- paste0(PARAMS$outer_cv_nfolds, "-fold CV, repeated ", PARAMS$num_outer_rep, " times")
    
    title <- paste0(plot.title, "\n(", outer.loop.str, ")")
    
    #########################################################
    
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
      
      ## Evaluate at each time point
      survivalROC_data[[j]] <- data_frame(t = predict.time)%>%
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
    
    # Calculate C-index
    
    all.cindex <- sapply(as.data.frame(nivo), function (x) {
      calc <- concordance.index (x, surv.time = this.preds[which(this.preds$ARM == "Nivo"), "Stime"],
                                 surv.event = this.preds[which(this.preds$ARM == "Nivo"),"status"])
      return(calc$c.index)
    })
    
    nivo.dat[[i]] <- mutate(med_pred, model = deparse(marker[[i]]$formula), Cindex = median(all.cindex), ARM = "Nivo")
    
    plot.dat[[i]] <- rbind(combo.dat[[i]], nivo.dat[[i]])
  }
  
  plot.dat <- do.call("rbind", plot.dat)
  
  model_levels <- rep(NA, length(marker))
  for (i in 1:length(marker)) {
    model_levels[i] <- deparse(marker[[i]]$formula)
  }
  plot.dat$model <- factor(plot.dat$model, levels = model_levels )
  
  # Alternate x-axis location for annotations
  xtext <- rep(c(0.75, 0.85), each = length(predict.time), times = length(marker))
  xtext <- xtext[1:(length(xtext)/2)]
  
  # Create list to store plots
  perf <- list()  
  
  plot.dat_g <- group_by(plot.dat, ARM, model)
  Cindex.dat <- summarise(plot.dat_g, Cindex = round(mean(Cindex), 2)) %>%
    ungroup()
  
  # plot.dat$model <- factor(plot.dat$model, levels = )
  
  perf[["ROC_plot"]] <- ggplot(plot.dat, mapping = aes(x = FPR, y = TPR, color = model)) +
    # geom_line() +
    # scale_linetype() +
    scale_color_discrete() +
    coord_fixed() +
    geom_label(data = plot.dat %>% dplyr::select(t, model, ARM, AUC) %>% unique %>%
                 arrange(ARM, model, t) %>%
                 mutate(y = rep(seq(0.4, 0.05, length.out = length(marker)),
                                times = length(levels(as.factor(plot.dat$ARM))),
                                each = length(predict.time)),
                        x = 0.75),
               mapping = aes(x = x, y = y, label = sprintf("%.3f", AUC)), size = rel(3), show.legend = F) +
    annotate("text", x=0.65, y=0.5, label="Med AUC", hjust=0, size=rel(4)) +
    geom_point(aes(shape=model), size = rel(2)) +
    scale_shape_manual(values = seq(49, length.out = length(unique(plot.dat$model)))) +
    guides(shape = guide_legend(override.aes = list(size=5))) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), size = .8, col = "grey") +
    facet_grid(ARM ~ t, labeller = as_labeller(time_labels)) +
    labs(title = title) +
    theme_bw(PARAMS$font_size) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          # plot.subtitle = element_text(hjust = 0.5, size = rel(.9)),
          strip.background = element_blank())
  
  # g = ggplotGrob(p)
  # g$layout$clip[g$layout$name == "panel"] = "off"
  # perf[["ROC_plot"]] <- grid.draw(g)
  
  perf[["Cindex_Combo"]] <- filter(Cindex.dat, ARM == "Combo")
  perf[["Cindex_Combo"]] <- dplyr::select(perf[["Cindex_Combo"]], -ARM)
  perf[["Cindex_Nivo"]] <- filter(Cindex.dat, ARM == "Nivo")
  perf[["Cindex_Nivo"]] <- dplyr::select(perf[["Cindex_Nivo"]], -ARM)
  
  return(perf)
}
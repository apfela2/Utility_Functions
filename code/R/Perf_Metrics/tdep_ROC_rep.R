# Function to plot time-dependent ROC curves for matrix of predictions


# predict.time is vector of times you want to plot ROC curves for
# Stime is survival time
# status is censoring indicator where 1 means event and 0 is censored
# marker is your matrix of predicted biomarker scores
# data is your data object

tdep_ROC_rep <- function(predict.time = c(91.3, 182.6, 365.25, 547.8, 730.5), time_labels = c('91.3' = "3 Months", '182.6' = "6 Months", '365.25' = "12 Months", '547.8' = "18 Months", '730.5' = "24 Months"),
                         Stime, status, marker, data, plot.title) {
  
  
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
    
    
    ## Evaluate at each time point
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
  
  
  # Calculate C-index
  
  cindex <- sapply(as.data.frame(marker), function (x) {
    calc <- concordance.index (x, surv.time = data[[Stime]], surv.event = data[[status]])
    return(calc$c.index)
  })
  
  
  ## Plot
  survivalROC_data_long <- bind_rows(survivalROC_data, .id = "rep")
  
  survivalROC_data_long_g <- group_by(survivalROC_data_long, t)
  med_auc <- summarise(survivalROC_data_long_g, med_auc = median(auc))
  
  
  
  plot.subtitle <- paste("Mean C-index =", round(mean(cindex), 2), "Std =", round(sd(cindex), 2), "across", ncol(marker), "repeats")
  
  survivalROC_data_long %>%
    ggplot( mapping = aes(x = FP, y = TP)) +
    geom_point(aes(group = rep), alpha = .1, size = .3) +
    geom_line(aes(group = rep), alpha = .1, size = .3) +
    geom_smooth(col = "red") + 
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), size = .8, col = "grey") +
    geom_label(data = med_auc %>% dplyr::select(t,med_auc) %>% unique,
               mapping = aes(label = paste("Median AUC = ", sprintf("%.3f", med_auc))), x = 0.75, y = 0.125) +
    facet_wrap( ~ t, labeller = as_labeller(time_labels)) +
    labs(title = plot.title, subtitle = plot.subtitle) +
    theme_bw(PARAMS$font_size) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = rel(.9)),
          strip.background = element_blank())
}



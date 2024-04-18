## Make plot of AUC vs time for matrix of predictions


time_AUC_rep <- function(predict.time = seq(30.4375, 730.5, 30.4375), Stime, status, marker, data, plot.title = paste("Time vs AUC Plot")) {
  
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


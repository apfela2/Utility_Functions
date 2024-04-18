# Function to plot time-dependent ROC curves for vector of predictions


# predict.time is vector of times you want to plot ROC curves for
# Stime is survival time
# status is censoring indicator where 1 means event and 0 is censored
# marker is your vector of predicted biomarker score
# data is your data object

tdep_ROC <- function(predict.time = c(91.3, 182.6, 365.25, 547.8, 730.5), time_labels = c('91.3' = "3 Months", '182.6' = "6 Months", '365.25' = "12 Months", '547.8' = "18 Months", '730.5' = "24 Months"),
                     Stime, status, marker, data, plot.title) {
  
  survivalROC_helper <- function(t) {
    survivalROC(Stime        = data[[Stime]],
                status       = data[[status]],
                marker       = marker,
                predict.time = t,
                method       = "NNE",
                span = 0.25 * nrow(data)^(-0.20))
  }
  
  
  ## Evaluate every 20 days
  survivalROC_data <- tibble(t = predict.time) %>%
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
  
  # Calculate C-index
  
  cindex <- concordance.index (marker, surv.time = data[[Stime]], surv.event = data[[status]])
  cindex <- cindex$c.index
  
  plot.subtitle <- paste("C-index =", round(cindex, 2))
  
  ## Plot
  survivalROC_data %>%
    ggplot( mapping = aes(x = FP, y = TP)) +
    # geom_point() +
    geom_line() +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1)) +
    geom_label(data = survivalROC_data %>% dplyr::select(t,auc) %>% unique,
               mapping = aes(label = paste("AUC = ", sprintf("%.3f", auc))), x = 0.75, y = 0.125) +
    facet_wrap( ~ t, labeller = as_labeller(time_labels)) +
    theme_bw(18) +
    labs(title = plot.title, subtitle = plot.subtitle) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = rel(.9)),
          strip.background = element_blank())
}

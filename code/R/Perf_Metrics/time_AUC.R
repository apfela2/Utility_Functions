## Make plot of AUC vs time

time_AUC <- function(predict.time = seq(30.4375, 730.5, 30.4375), Stime, status, marker, data, plot.title) {
  
  survivalROC_helper <- function(t) {
    survivalROC(Stime        = data[[Stime]],
                status       = data[[status]],
                marker       = marker,
                predict.time = t,
                method       = "NNE",
                span = 0.25 * nrow(data)^(-0.20))
  }
  
  
  ## Evaluate every Month
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
  
  # Plot
  survivalROC_data %>%
    ggplot( mapping = aes(x = t, y = auc)) +
    geom_point() +
    geom_line() +
    theme_bw(18) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank()) +
    xlab("Months") +
    ylab("Time-dependent AUC") +
    labs(title = plot.title) +
    scale_x_continuous(breaks = predict.time,
                       labels = seq(1, length(predict.time))) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    ylim(0,1)
}

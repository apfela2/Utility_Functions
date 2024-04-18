tdep_PPV_trt_biom <- function(predict.time = c(91.3, 182.6, 365.25, 547.8, 730.5), time_labels = c('91.3' = "3 Months", '182.6' = "6 Months", '365.25' = "12 Months", '547.8' = "18 Months", '730.5' = "24 Months"),
                              Stime, status, data, plot.title, cut_length = 100, xaxislab, trt = ACTARM, biom) {
  
  marker <- data[[biom]]
  
  
  
  Surv <- Surv(data[[Stime]], data[[status]])
  model <- coxph(Surv ~ data[[biom]], data = data)
  negative <- ( coef(model)[1] < 0 )
  
  if (negative) {
    marker <- -marker
  }
  
  # Change output from SeSpPPVNPV function to only return prevalence column instead of dataframe from Stats to aid in future looping
  SeSpPPVNPV2 <- function (cutpoint, T, delta, marker, other_markers = NULL, cause, 
                           weighting = "marginal", times, iid = FALSE) 
  {
    if (length(delta) != length(T) | length(marker) != length(T) | 
        length(delta) != length(T)) {
      stop("lengths of vector T, delta and marker have to be equal\n")
    }
    if (missing(times)) {
      stop("Choose at least one time for computing the time-dependent AUC\n")
    }
    if (!weighting %in% c("marginal", "cox", "aalen")) {
      stop("the weighting argument must be marginal (default), cox or aalen.\n")
    }
    if (weighting %in% c("cox", "aalen") & !missing(other_markers) & 
        !class(other_markers) == "matrix") {
      stop("argument other_markers must be a matrix\n")
    }
    if (weighting %in% c("cox", "aalen") & !missing(other_markers)) {
      if (!nrow(other_markers) == length(marker)) 
        stop("lengths of vector T, delta, marker and number of rows of other_markers have to be equal\n")
    }
    if (weighting != "marginal" & iid) {
      stop("Weighting must be marginal for computing the iid representation \n Choose iid=FALSE or weighting=marginal in the input arguments")
    }
    if (weighting %in% c("cox", "aalen") & !missing(other_markers)) {
      is_not_na <- as.logical(apply(!is.na(cbind(T, delta, 
                                                 marker, other_markers)), 1, prod))
      T <- T[is_not_na]
      delta <- delta[is_not_na]
      marker <- marker[is_not_na]
      other_markers <- as.matrix(other_markers[is_not_na, ])
    }
    else {
      is_not_na <- as.logical(apply(!is.na(cbind(T, delta, 
                                                 marker)), 1, prod))
      T <- T[is_not_na]
      delta <- delta[is_not_na]
      marker <- marker[is_not_na]
    }
    start_computation_time <- Sys.time()
    n <- length(T)
    n_times <- length(times)
    if (n_times == 1) {
      times <- c(0, times)
      n_times <- 2
    }
    times <- times[order(times)]
    times_names <- paste("t=", times, sep = "")
    CumInci <- rep(NA, n_times)
    surv <- rep(NA, n_times)
    names(CumInci) <- times_names
    names(surv) <- times_names
    Stats <- matrix(NA, nrow = n_times, ncol = 6)
    colnames(Stats) <- c("Cases", "survivor at t", "Other events at t", 
                         "Censored at t", "Positive (X>c)", "Negative (X<=c)")
    rownames(Stats) <- times_names
    Stats[, c("Positive (X>c)", "Negative (X<=c)")] <- matrix(c(sum(marker > 
                                                                      cutpoint), sum(marker <= cutpoint)), byrow = TRUE, ncol = 2, 
                                                              nrow = nrow(Stats))
    order_T <- order(T)
    T <- T[order_T]
    delta <- delta[order_T]
    marker <- marker[order_T]
    if (weighting == "marginal") {
      weights <- pec::ipcw(Surv(failure_time, status) ~ 1, 
                           data = data.frame(failure_time = T, status = as.numeric(delta != 
                                                                                     0)), method = "marginal", times = times, subjectTimes = T, 
                           subjectTimesLag = 1)
    }
    if (weighting == "cox") {
      if (missing(other_markers)) {
        marker_censoring <- marker
      }
      other_markers <- other_markers[order_T, ]
      marker_censoring <- cbind(marker, other_markers)
      colnames(marker_censoring) <- paste("X", 1:ncol(marker_censoring), 
                                          sep = "")
      fmla <- as.formula(paste("Surv(T,status)~", paste(paste("X", 
                                                              1:ncol(marker_censoring), sep = ""), collapse = "+")))
      data_weight <- as.data.frame(cbind(data.frame(T = T, 
                                                    status = as.numeric(delta != 0)), marker_censoring))
      weights <- pec::ipcw(fmla, data = data_weight, method = "cox", 
                           times = as.matrix(times), subjectTimes = data_weight[, 
                                                                                "T"], subjectTimesLag = 1)
    }
    if (weighting == "aalen") {
      if (missing(other_markers)) {
        marker_censoring <- marker
      }
      other_markers <- other_markers[order_T, ]
      marker_censoring <- cbind(marker, other_markers)
      colnames(marker_censoring) <- paste("X", 1:ncol(marker_censoring), 
                                          sep = "")
      fmla <- as.formula(paste("Surv(T,status)~", paste(paste("X", 
                                                              1:ncol(marker_censoring), sep = ""), collapse = "+")))
      data_weight <- as.data.frame(cbind(data.frame(T = T, 
                                                    status = as.numeric(delta != 0)), marker_censoring))
      weights <- pec::ipcw(fmla, data = data_weight, method = "aalen", 
                           times = as.matrix(times), subjectTimes = data_weight[, 
                                                                                "T"], subjectTimesLag = 1)
    }
    Mat.data <- cbind(T, delta, marker)
    colnames(Mat.data) <- c("T", "delta", "marker")
    Weights_cases_all <- 1/(weights$IPCW.subjectTimes * n)
    Weights_cases_all <- Weights_cases_all
    FP_1 <- rep(NA, n_times)
    TP <- rep(NA, n_times)
    FP_2 <- rep(NA, n_times)
    PPV <- rep(NA, n_times)
    NPV_1 <- rep(NA, n_times)
    NPV_2 <- rep(NA, n_times)
    names(FP_1) <- times_names
    names(TP) <- times_names
    names(FP_2) <- times_names
    names(PPV) <- times_names
    names(NPV_1) <- times_names
    names(NPV_2) <- times_names
    if (iid) {
      MatInt0TcidhatMCksurEff <- Compute.iid.KM(T, delta != 
                                                  0)
      epsiloni.Se <- matrix(NA, nrow = n, ncol = n_times)
      epsiloni.Sp1 <- matrix(NA, nrow = n, ncol = n_times)
      epsiloni.Sp2 <- matrix(NA, nrow = n, ncol = n_times)
      epsiloni.PPV <- matrix(NA, nrow = n, ncol = n_times)
      epsiloni.NPV1 <- matrix(NA, nrow = n, ncol = n_times)
      epsiloni.NPV2 <- matrix(NA, nrow = n, ncol = n_times)
      se.Se <- rep(NA, n_times)
      se.Sp1 <- rep(NA, n_times)
      se.Sp2 <- rep(NA, n_times)
      names(se.Sp1) <- times_names
      names(se.Sp2) <- times_names
      names(se.Se) <- times_names
      colnames(epsiloni.Se) <- times_names
      colnames(epsiloni.Sp1) <- times_names
      colnames(epsiloni.Sp2) <- times_names
      colnames(epsiloni.PPV) <- times_names
      colnames(epsiloni.NPV1) <- times_names
      colnames(epsiloni.NPV2) <- times_names
    }
    for (t in 1:n_times) {
      Cases <- (Mat.data[, "T"] < times[t] & Mat.data[, "delta"] == 
                  cause)
      Controls_1 <- (Mat.data[, "T"] > times[t])
      Controls_2 <- (Mat.data[, "T"] < times[t] & Mat.data[, 
                                                           "delta"] != cause & Mat.data[, "delta"] != 0)
      if (weights$method != "marginal") {
        Weights_controls_1 <- 1/(weights$IPCW.times[, t] * 
                                   n)
      }
      else {
        Weights_controls_1 <- rep(1/(weights$IPCW.times[t] * 
                                       n), times = n)
      }
      Weights_controls_1 <- Weights_controls_1
      Weights_cases <- Weights_cases_all
      Weights_controls_2 <- Weights_cases_all
      Weights_cases[!Cases] <- 0
      Weights_controls_1[!Controls_1] <- 0
      Weights_controls_2[!Controls_2] <- 0
      den_TP_t <- sum(Weights_cases)
      den_FP_1_t <- sum(Weights_controls_1)
      den_FP_2_t <- sum(Weights_controls_2) + sum(Weights_controls_1)
      den_PPV_t <- sum(Weights_cases[which(Mat.data[, "marker"] > 
                                             cutpoint)] + Weights_controls_1[which(Mat.data[, 
                                                                                            "marker"] > cutpoint)] + Weights_controls_2[which(Mat.data[, 
                                                                                                                                                       "marker"] > cutpoint)])
      den_NPV_t <- sum(Weights_cases[which(Mat.data[, "marker"] <= 
                                             cutpoint)] + Weights_controls_1[which(Mat.data[, 
                                                                                            "marker"] <= cutpoint)] + Weights_controls_2[which(Mat.data[, 
                                                                                                                                                        "marker"] <= cutpoint)])
      if (den_TP_t != 0) {
        TP[t] <- sum(Weights_cases[which(Mat.data[, "marker"] > 
                                           cutpoint)])/den_TP_t
      }
      if (den_FP_1_t != 0) {
        FP_1[t] <- sum(Weights_controls_1[which(Mat.data[, 
                                                         "marker"] > cutpoint)])/den_FP_1_t
      }
      if (den_FP_2_t != 0) {
        FP_2[t] <- sum(Weights_controls_1[which(Mat.data[, 
                                                         "marker"] > cutpoint)] + Weights_controls_2[which(Mat.data[, 
                                                                                                                    "marker"] > cutpoint)])/den_FP_2_t
        NPV_2[t] <- ((1 - FP_2[t]) * den_FP_2_t)/den_NPV_t
      }
      CumInci[t] <- c(den_TP_t)
      surv[t] <- c(den_FP_1_t)
      Stats[t, 1:4] <- c(sum(Cases), sum(Controls_1), sum(Controls_2), 
                         n - sum(Cases) - sum(Controls_1) - sum(Controls_2))
      PPV[t] <- (TP[t] * CumInci[t])/den_PPV_t
      NPV_1[t] <- ((1 - FP_1[t]) * surv[t])/den_NPV_t
      if (iid) {
        Int0tdMCsurEffARisk <- NA
        Int0tdMCsurEffARisk <- MatInt0TcidhatMCksurEff[max(which(Mat.data[, 
                                                                          "T"] <= times[t])), ]
        Weights_cases <- Weights_cases
        Weights_controls_2 <- Weights_controls_2
        Weights_controls_1 <- Weights_controls_1
        epsiloni.Se[, t] <- Weights_cases * n * (as.numeric(Mat.data[, 
                                                                     "marker"] > cutpoint) - TP[t])/CumInci[t] + colMeans(MatInt0TcidhatMCksurEff * 
                                                                                                                            (Weights_cases * n * (as.numeric(Mat.data[, "marker"] > 
                                                                                                                                                               cutpoint) - TP[t])/CumInci[t]))
        epsiloni.Sp1[, t] <- (n/sum(Controls_1)) * as.numeric(Mat.data[, 
                                                                       "T"] > times[t]) * (as.numeric(Mat.data[, "marker"] <= 
                                                                                                        cutpoint) - (1 - FP_1[t]))
        epsiloni.Sp2[, t] <- (Weights_controls_2 * n + Weights_controls_1 * 
                                n) * (as.numeric(Mat.data[, "marker"] <= cutpoint) - 
                                        (1 - FP_2[t]))/(1 - CumInci[t]) + colMeans(MatInt0TcidhatMCksurEff * 
                                                                                     ((Weights_controls_2 * n) * (as.numeric(Mat.data[, 
                                                                                                                                      "marker"] <= cutpoint) - (1 - FP_2[t]))/(1 - 
                                                                                                                                                                                 CumInci[t]))) + mean((Weights_controls_1 * 
                                                                                                                                                                                                         n) * (as.numeric(Mat.data[, "marker"] <= cutpoint) - 
                                                                                                                                                                                                                 (1 - FP_2[t]))/(1 - CumInci[t])) * Int0tdMCsurEffARisk
        epsiloni.PPV[, t] <- (as.numeric(Mat.data[, "marker"] > 
                                           cutpoint)/den_PPV_t) * ((Weights_cases * n + 
                                                                      Weights_controls_2 * n) * (as.numeric(Mat.data[, 
                                                                                                                     "delta"] == cause) - as.numeric(Mat.data[, "delta"] != 
                                                                                                                                                       0) * PPV[t]) - Weights_controls_1 * n * PPV[t]) + 
          colMeans(MatInt0TcidhatMCksurEff * ((as.numeric(Mat.data[, 
                                                                   "marker"] > cutpoint)/den_PPV_t) * (Weights_cases * 
                                                                                                         n + Weights_controls_2 * n) * (as.numeric(Mat.data[, 
                                                                                                                                                            "delta"] == cause) - as.numeric(Mat.data[, 
                                                                                                                                                                                                     "delta"] != 0) * PPV[t]))) - mean((as.numeric(Mat.data[, 
                                                                                                                                                                                                                                                            "marker"] > cutpoint)/den_PPV_t) * Weights_controls_1 * 
                                                                                                                                                                                                                                         n * PPV[t]) * Int0tdMCsurEffARisk
        epsiloni.NPV2[, t] <- (as.numeric(Mat.data[, "marker"] <= 
                                            cutpoint)/den_NPV_t) * ((Weights_cases * n + 
                                                                       Weights_controls_2 * n) * (as.numeric(Mat.data[, 
                                                                                                                      "delta"] != 0 & Mat.data[, "delta"] != cause) - 
                                                                                                    as.numeric(Mat.data[, "delta"] != 0) * NPV_2[t]) + 
                                                                      Weights_controls_1 * n * (1 - NPV_2[t])) + colMeans(MatInt0TcidhatMCksurEff * 
                                                                                                                            ((as.numeric(Mat.data[, "marker"] <= cutpoint)/den_NPV_t) * 
                                                                                                                               (Weights_cases * n + Weights_controls_2 * n) * 
                                                                                                                               (as.numeric(Mat.data[, "delta"] != 0 & Mat.data[, 
                                                                                                                                                                               "delta"] != cause) - as.numeric(Mat.data[, 
                                                                                                                                                                                                                        "delta"] != 0) * NPV_2[t]))) + mean((as.numeric(Mat.data[, 
                                                                                                                                                                                                                                                                                 "marker"] <= cutpoint)/den_NPV_t) * Weights_controls_1 * 
                                                                                                                                                                                                                                                              n * (1 - NPV_2[t])) * Int0tdMCsurEffARisk
        epsiloni.NPV1[, t] <- (as.numeric(Mat.data[, "marker"] <= 
                                            cutpoint)/den_NPV_t) * ((Weights_cases * n + 
                                                                       Weights_controls_2 * n) * (-as.numeric(Mat.data[, 
                                                                                                                       "delta"] != 0) * NPV_1[t]) + Weights_controls_1 * 
                                                                      n * (1 - NPV_1[t])) + colMeans(MatInt0TcidhatMCksurEff * 
                                                                                                       ((as.numeric(Mat.data[, "marker"] <= cutpoint)/den_NPV_t) * 
                                                                                                          (Weights_cases * n + Weights_controls_2 * n) * 
                                                                                                          (-as.numeric(Mat.data[, "delta"] != 0) * NPV_1[t]))) + 
          mean((as.numeric(Mat.data[, "marker"] <= cutpoint)/den_NPV_t) * 
                 Weights_controls_1 * n * (1 - NPV_1[t])) * 
          Int0tdMCsurEffARisk
      }
    }
    if (iid) {
      se.Se <- apply(epsiloni.Se, 2, sd)/sqrt(n)
      se.Sp1 <- apply(epsiloni.Sp1, 2, sd)/sqrt(n)
      se.Sp2 <- apply(epsiloni.Sp2, 2, sd)/sqrt(n)
      se.PPV <- apply(epsiloni.PPV, 2, sd)/sqrt(n)
      se.NPV1 <- apply(epsiloni.NPV1, 2, sd)/sqrt(n)
      se.NPV2 <- apply(epsiloni.NPV2, 2, sd)/sqrt(n)
    }
    if (iid) {
      inference <- list(mat_iid_rep_Se = epsiloni.Se, mat_iid_rep_Sp1 = epsiloni.Sp1, 
                        mat_iid_rep_Sp2 = epsiloni.Sp2, mat_iid_rep_PPV = epsiloni.PPV, 
                        mat_iid_rep_NPV1 = epsiloni.NPV1, mat_iid_rep_NPV2 = epsiloni.NPV2, 
                        vect_se_Se = se.Se, vect_se_Sp1 = se.Sp1, vect_se_Sp2 = se.Sp2, 
                        vect_se_PPV = se.PPV, vect_se_NPV1 = se.NPV1, vect_se_NPV2 = se.NPV2)
    }
    else {
      inference <- NA
    }
    stop_computation_time <- Sys.time()
    if (max(Stats[, 3]) == 0) {
      out <- list(TP = TP, FP = FP_1, PPV = PPV, NPV = NPV_1, 
                  Stats = Stats[, 5], inference = inference, computation_time = difftime(stop_computation_time, 
                                                                                         start_computation_time, units = "secs"), iid = iid, 
                  n = n, times = times, weights = weights, cutpoint = cutpoint)
      class(out) <- "ipcwsurvivalSeSpPPVNPV"
      out
    }
    else {
      out <- list(TP = TP, FP_1 = FP_1, FP_2 = FP_2, PPV = PPV, 
                  NPV_1 = NPV_1, NPV_2 = NPV_2, Stats = Stats, inference = inference, 
                  computation_time = difftime(stop_computation_time, 
                                              start_computation_time, units = "secs"), iid = iid, 
                  n = n, times = times, weights = weights, cutpoint = cutpoint)
      class(out) <- "ipcwcompetingrisksSeSpPPVNPV"
      out
    }
  }
  
  
  
  
  PPV_data <- list()
  
  PPV_helper <- function(cut) {
    SeSpPPVNPV2( cutpoint = cut,
                 T        = data[[Stime]],
                 delta    = data[[status]],
                 marker   = marker,
                 cause    = 1,
                 times    = predict.time,
                 iid      = F)
  }
  
  
  ## Evaluate at chosen time points
  PPV_data <- tibble(cut = seq(min(marker), max(marker), length.out = cut_length)) %>%
    mutate(PPV = purrr::map(cut, PPV_helper),
           ## Extract scalar AUC
           # PPV = purrr::map_dbl(PPV, magrittr::extract2, "PPV"),
           # NPV = purrr::map_dbl(PPV, magrittr::extract2, "NPV"),
           # Prevalence = purrr::map_dbl(PPV, function(obj){
           #   magrittr::extract2("Stats")[((length(times)*3)+1):(length(times)*4)]),
           # Prevalence1 = as.vector(purrr::map_dfc(PPV, magrittr::extract2, "Stats")[2,4]),
           # times = purrr::map_dbl(PPV, magrittr::extract2, "times"),
           # cutpoint = purrr::map_dbl(PPV, magrittr::extract2, "cutpoint")) %>%
           df_PPV = purrr::map(PPV, function(obj) {
             # as_tibble(obj[c("PPV", "NPV", "times")])
             as_tibble(obj[c("PPV", "NPV", "times", Prevalence = "Stats")])
           })) %>%
    # Stats = purrr::map_at(PPV, c(((length(times)*3)+1):(length(times)*4)), function(obj) {
    #    as_tibble(obj["Stats"])
    # })) %>%
    dplyr::select(-PPV) %>%
    unnest() %>%
    # arrange(times, cut, PPV) %>%
    rename(Prevalence = Stats)
  
  
  
  ## Set up data for plot
  
  # Make prevalence a proportion
  PPV_data$Prevalence <- PPV_data$Prevalence/nrow(data)
  
  if (negative) {
    PPV_data$cut <- -PPV_data$cut
    PPV_data$Prevalence <- 1 - PPV_data$Prevalence
  }
  
  
  PPV_data <- tidyr::gather(PPV_data, key = "Metric", value = "Performance", Prevalence)
  
  agg.dat <- PPV_data %>%
    group_by(times, Metric, cut) %>%
    data.frame(stringsAsFactors = FALSE)
  
  drops <- c("PPV", "NPV")
  agg.dat <- agg.dat[, !(names(agg.dat)) %in% drops, drop = FALSE]
  
  ##########################################
  
  
  
  
  ONE <- data[as.numeric(data[[trt]]) == 1,]
  TWO <- data[as.numeric(data[[trt]]) == 2,]
  
  TRT_labels <- c(levels(data[[trt]]))
  
  
  marker <- ONE[[biom]]
  
  Surv <- Surv(ONE[[Stime]], ONE[[status]])
  
  model <- coxph(Surv ~ ONE[[biom]], data = ONE)
  
  negative <- ( coef(model)[1] < 0 )
  
  if (negative) {
    marker <- -marker
  }
  
  
  PPV_data <- list()
  
  PPV_helper <- function(cut) {
    SeSpPPVNPV2( cutpoint = cut,
                 T        = ONE[[Stime]],
                 delta    = ONE[[status]],
                 marker   = marker,
                 cause    = 1,
                 times    = predict.time,
                 iid      = F)
  }
  
  
  ## Evaluate at chosen time points
  PPV_data <- tibble(cut = seq(min(marker), max(marker), length.out = cut_length)) %>%
    mutate(PPV = purrr::map(cut, PPV_helper),
           ## Extract scalar AUC
           # PPV = purrr::map_dbl(PPV, magrittr::extract2, "PPV"),
           # NPV = purrr::map_dbl(PPV, magrittr::extract2, "NPV"),
           # Prevalence = purrr::map_dbl(PPV, function(obj){
           #   magrittr::extract2("Stats")[((length(times)*3)+1):(length(times)*4)]),
           # Prevalence1 = as.vector(purrr::map_dfc(PPV, magrittr::extract2, "Stats")[2,4]),
           # times = purrr::map_dbl(PPV, magrittr::extract2, "times"),
           # cutpoint = purrr::map_dbl(PPV, magrittr::extract2, "cutpoint")) %>%
           df_PPV = purrr::map(PPV, function(obj) {
             # as_tibble(obj[c("PPV", "NPV", "times")])
             as_tibble(obj[c(PPV_ONE = "PPV", NPV_ONE = "NPV", "times")])
           })) %>%
    # Stats = purrr::map_at(PPV, c(((length(times)*3)+1):(length(times)*4)), function(obj) {
    #    as_tibble(obj["Stats"])
    # })) %>%
    dplyr::select(-PPV) %>%
    unnest() %>%
    # arrange(times, cut, PPV) %>%
    rename(PPV_ONE = PPV, NPV_ONE = NPV)
  
  
  
  ## Set up data for plot
  
  if (negative) {
    PPV_data$cut <- -PPV_data$cut
  }
  
  
  
  PPV_data <- tidyr::gather(PPV_data, key = "Metric", value = "Performance", PPV_ONE, NPV_ONE)
  
  PPV_data <- dplyr::mutate(.data = PPV_data, Metric = ifelse(Metric == "PPV_ONE", "NPV_ONE",
                                                              ifelse(Metric == "NPV_ONE", "PPV_ONE", Metric)))
  
  
  
  
  agg.dat_ONE <- PPV_data %>%
    group_by(times, Metric, cut) %>%
    data.frame(stringsAsFactors = FALSE)
  
  
  ##############################################################################
  
  marker <- TWO[[biom]]
  
  Surv <- Surv(TWO[[Stime]], TWO[[status]])
  
  model <- coxph(Surv ~ TWO[[biom]], data = TWO)
  
  negative <- ( coef(model)[1] < 0 )
  
  if (negative) {
    marker <- -marker
  }
  
  
  PPV_data <- list()
  
  PPV_helper <- function(cut) {
    SeSpPPVNPV2( cutpoint = cut,
                 T        = TWO[[Stime]],
                 delta    = TWO[[status]],
                 marker   = marker,
                 cause    = 1,
                 times    = predict.time,
                 iid      = F)
  }
  
  
  ## Evaluate at chosen time points
  PPV_data <- tibble(cut = seq(min(marker), max(marker), length.out = cut_length)) %>%
    mutate(PPV = purrr::map(cut, PPV_helper),
           ## Extract scalar AUC
           # PPV = purrr::map_dbl(PPV, magrittr::extract2, "PPV"),
           # NPV = purrr::map_dbl(PPV, magrittr::extract2, "NPV"),
           # Prevalence = purrr::map_dbl(PPV, function(obj){
           #   magrittr::extract2("Stats")[((length(times)*3)+1):(length(times)*4)]),
           # Prevalence1 = as.vector(purrr::map_dfc(PPV, magrittr::extract2, "Stats")[2,4]),
           # times = purrr::map_dbl(PPV, magrittr::extract2, "times"),
           # cutpoint = purrr::map_dbl(PPV, magrittr::extract2, "cutpoint")) %>%
           df_PPV = purrr::map(PPV, function(obj) {
             # as_tibble(obj[c("PPV", "NPV", "times")])
             as_tibble(obj[c(PPV_TWO = "PPV", NPV_TWO = "NPV", "times")])
           })) %>%
    # Stats = purrr::map_at(PPV, c(((length(times)*3)+1):(length(times)*4)), function(obj) {
    #    as_tibble(obj["Stats"])
    # })) %>%
    dplyr::select(-PPV) %>%
    unnest() %>%
    # arrange(times, cut, PPV) %>%
    rename(PPV_TWO = PPV, NPV_TWO = NPV)
  
  
  
  ## Set up data for plot
  
  if (negative) {
    PPV_data$cut <- -PPV_data$cut
  }
  
  
  
  PPV_data <- tidyr::gather(PPV_data, key = "Metric", value = "Performance", PPV_TWO, NPV_TWO)
  
  PPV_data <- dplyr::mutate(.data = PPV_data, Metric = ifelse(Metric == "PPV_TWO", "NPV_TWO",
                                                              ifelse(Metric == "NPV_TWO", "PPV_TWO", Metric)))
  
  
  
  
  agg.dat_TWO <- PPV_data %>%
    group_by(times, Metric, cut) %>%
    data.frame(stringsAsFactors = FALSE)
  
  agg.dat <- rbind(agg.dat, agg.dat_ONE, agg.dat_TWO)
  
  
  x.min.breaks <-  seq(round(min(agg.dat$cut, na.rm = T), 1),
                       round(max(agg.dat$cut, na.rm = T), 1), 0.2)
  y.min.breaks <- seq(0, 1, 0.05)
  
  
  
  
  # Make plot
  ggplot(agg.dat, aes(y = Performance, x = cut, col = Metric, fill = Metric)) + 
    geom_line(size = 1) +
    # geom_line(aes(y = Upper, col = Metric),size = 1, linetype = "dotted") +
    # geom_line(aes(y = Lower, col = Metric),size = 1, linetype = "dotted") +
    # 
    scale_x_continuous(minor_breaks = x.min.breaks) +
    scale_y_continuous(minor_breaks = y.min.breaks) +
    scale_color_discrete(breaks = c("PPV_ONE", "NPV_ONE", "PPV_TWO", "NPV_TWO", "Prevalence"),
                         labels = c(paste("PPV :", TRT_labels[1]), paste("NPV :", TRT_labels[1]), paste("PPV :", TRT_labels[2]), paste("NPV :", TRT_labels[2]), "Prevalence")) +
    ylab("Performance") +
    xlab(xaxislab) +
    #stat_smooth(se = FALSE, level = .5) + 
    labs(title = plot.title) +
    theme_bw(18) +
    facet_wrap( ~ times, labeller = as_labeller(time_labels)) +
    theme(plot.title = element_text(hjust = 0.5))
  #           plot.subtitle = element_text(hjust = 0.5, size = rel(.9)))
}


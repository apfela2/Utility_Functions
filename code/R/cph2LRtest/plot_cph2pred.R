
# x = object from cph2predLRtest
# ylimA = c(lower_limit, upper_limit) for 3-way interaction plot for aVar faceted by bVar
# ylimB = c(lower_limit, upper_limit) for 3-way interaction plot for bVar faceted by aVar
# confInt = Should plot include confidence intervals

plot.cph2pred <- function(x, ylimA = NULL, ylimB = NULL, confInt = TRUE){
  if(!is(x, "cph2pred")) {
    stop("x must be a 'cph2pred' object.")
  }
  
  ##########################################
  # Set up data for partial effects plots
  ##########################################
  
  ### A. set up dummy values/reference levels: 
  
  ### use median for continuous values and 0 for factors
  data <- x$argumentList$data
  var_mod <- colnames(data)[colnames(data) %in% c(x$argumentList$aVar,
     x$argumentList$bVar, x$argumentList$covars) ]
  
  # Identify which variables are numeric or categorical
  isnum <- unlist(lapply(data[, var_mod], is.numeric))
  var_num <- var_mod[isnum]
  var_cat <- var_mod[!isnum]
  # var_num <- var_mod[sapply(data[,var_mod], is.numeric)]
  # var_cat <- var_mod[!sapply(data[,var_mod], is.numeric)]
  
  if(length(x$argumentList$trtVar) == 0) {
    var_trt <- "SINGLE_ARM"
    data[[var_trt]] <- factor("DUMMY")
  } else {
    var_trt <- x$argumentList$trtVar  
  }
  
  
  var_ab <- c(x$argumentList$aVar, x$argumentList$bVar)
  var_cov <- x$argumentList$covars[x$argumentList$covars %in% var_num]
  
  
  {
#######################
    # Rewritten: dplyr
    #     # ## Get median of all numerical vars
    # meds=data %>% dplyr::select(all_of(var_num))%>%
    #   summarise(across(var_num, list(med=median), na.rm=T))
    # colnames(meds)=sub('_med', '', colnames(meds))

    # Rewritten: dplyr
    
###################
    # Select numeric variables from data
    num_vars <- data[, var_num]
    # Calculate median for each numeric variable
    meds <- apply(num_vars, 2, median, na.rm = TRUE)
    # Convert meds to a dataframe and rename columns
    meds <- data.frame(meds)
    meds <- as.data.frame(t(meds))
######################
    # dplyr
    # # Calculate quantiles for aVar and bVar
    # quants <- data %>% dplyr::select(all_of(var_ab))
    # quants <- lapply(quants, quantile, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
    # quants <- data.frame(quants)
    # quants <- dplyr::select(meds, -all_of(var_ab)) %>% bind_cols(quants)
    # dplyr
    #########################
    
    # Calculate quantiles for aVar and bVar
    quants <- data[, var_ab]
    quants <- lapply(quants, quantile, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
    quants <- as.data.frame(quants)
    quants <- cbind(meds[, -which(names(meds) %in% var_ab)], quants)
    #############################
    # dplyr
    # ## Set up categorical vars to later pull reference levels
    # ref_c=lapply(data %>% dplyr::select(all_of(var_cat)) , 
    #              function(x)levels(x))
    # ref_c=do.call(cbind, ref_c) %>% data.frame()
    # for (i in colnames(ref_c))
    #   ref_c[[i]]=factor(ref_c[[i]], levels=levels(data[[i]]))
    # dplyr
    
    ###############################
    # Get reference levels for categorical variables
    ref_c <- lapply(data[, var_cat, drop = F], function(x) levels(x))
    ref_c <- do.call(cbind, ref_c)
    ref_c <- as.data.frame(ref_c)
    # Set reference levels for each categorical variable
    for (i in colnames(ref_c)) {
      ref_c[[i]] <- factor(ref_c[[i]], levels = levels(data[[i]]))
    }
    ################################
  }
  
  ### parameters: 
  nintp<-10  ## number of interpolated numbers 
  pred_dat <- list()
  pred_dat_names <- list()
  
  
  # Capture number of treatment arms
    narms <- length(levels(as.factor(data[[var_trt]])))
  
    
    # For each BM we create 2 (one for each ARM) sequences from 10% to 90% followed by median for remaining BM
    pred_dat<-NULL
    for(i in var_num)
    {
      ra <- quantile(data[[i]], probs = c(0.1, 0.9))
      breaks <- seq(ra[1], ra[2], length.out=nintp)
      
      
      ## new_data_matrix for imputing var i
      
    
      mat_i <- data.frame(biom=i, TRT=rep(levels(data[[var_trt]]), each=nintp), rep(breaks, narms))
      colnames(mat_i)[ncol(mat_i)] <- i
      
      if(x$argumentList$withInteractions){
        ### add default value for other bioms, or set reference level for categorical vars
        if (i %in% var_ab){
          #####################
          # # dplyr
          # mat_i_dummy <- quants %>% bind_cols(ref_c[1, , drop = F] )%>% dplyr::select(-any_of(colnames(mat_i)))
          # mat_i_dummy <- lapply(mat_i_dummy, rep, times = nintp*narms)
          # mat_i <- lapply(mat_i, rep, each = 5)
          # mat_i=mat_i %>% bind_cols(mat_i_dummy)
          # # dplyr
          #######################
          
          # Combine quantiles and reference levels
          mat_i_dummy <- cbind(quants, ref_c[1, , drop = F])
          mat_i_dummy <- mat_i_dummy[, !(colnames(mat_i_dummy) %in% colnames(mat_i)), drop = F]
          
          # Repeat columns of mat_i_dummy and mat_i
          mat_i_dummy <- data.frame(lapply(mat_i_dummy, rep, times = nintp * narms))
          mat_i <- data.frame(lapply(mat_i, rep, each = 5))
          
          # Combine mat_i and mat_i_dummy
          mat_i <- cbind(mat_i, mat_i_dummy)
          ##############################
        } else {
        ############################
          # # dplyr
          # mat_i_dummy=meds %>% bind_cols(ref_c[1, , drop = F] )%>% dplyr::select(-any_of(colnames(mat_i)))
          # mat_i=mat_i %>% bind_cols(mat_i_dummy)
          # # dplyr
          ##############################
          # Combine meds and reference levels
          mat_i_dummy <- cbind(meds, ref_c[1, , drop = F])
          mat_i_dummy <- mat_i_dummy[, !(colnames(mat_i_dummy) %in% colnames(mat_i))]
        
            # Combine mat_i and mat_i_dummy
          mat_i <- cbind(mat_i, mat_i_dummy)
          #################################
        }
      } else {
        ##################################
        # # dplyr
        # mat_i_dummy=meds %>% bind_cols(ref_c[1, , drop = F] )%>% dplyr::select(-any_of(colnames(mat_i)))
        # mat_i=mat_i %>% bind_cols(mat_i_dummy)
        # # dplyr
        ####################################
        # Combine meds and reference levels
        mat_i_dummy <- cbind(meds, ref_c[1, , drop = F])
        mat_i_dummy <- mat_i_dummy[, !(colnames(mat_i_dummy) %in% colnames(mat_i))]
       
        # Combine mat_i and mat_i_dummy
        mat_i <- cbind(mat_i, mat_i_dummy)
        ###############################
      }
      
      
      if(is.null(pred_dat)){
        pred_dat=mat_i
      } else {
        ############################
        # dplyr
        # pred_dat=pred_dat %>% bind_rows(mat_i)
        # dplyr
        ############################
        # Combine pred_dat and mat_i
        pred_dat <- rbind(pred_dat, mat_i)
        ###############################
      }
    }
    
  
    
    ## for categorical vars
    for(i in var_cat)
    {
      ra= levels(data[[i]])
      
      ## new_data_matrix for imputing var i
      mat_i <- data.frame(biom=i, TRT=rep(levels(data[[var_trt]]), each=length(ra)) ,rep(ra, narms))
      colnames(mat_i)[ncol(mat_i)] <- i
      
      #############################################

      # # dplyr
      # ### add default value for other bioms, or set reference level for categorical vars
      # mat_i_dummy=meds %>% bind_cols(ref_c[1, , drop = F] )%>% dplyr::select(-any_of(colnames(mat_i)))
      # mat_i=mat_i %>% bind_cols(mat_i_dummy)
      # # dplyr
      ###############################################
      # Combine meds and reference levels
      mat_i_dummy <- cbind(meds, ref_c[1, , drop = F])
      mat_i_dummy <- mat_i_dummy[, !(colnames(mat_i_dummy) %in% colnames(mat_i))]
      
      # Combine mat_i and mat_i_dummy
      mat_i <- cbind(mat_i, mat_i_dummy)
      ################################################
      
      if(is.null(pred_dat)){
        pred_dat <- mat_i
      } else {
        ##########################
        # # dplyr
        # pred_dat=pred_dat %>% bind_rows(mat_i)
        # # dplyr
        ############################
        # Combine pred_dat and mat_i
        pred_dat <- rbind(pred_dat, mat_i)
        ###############################
      }
    }
    
    
    
    
    ### B: 
    # add basis functions for splines transformation of continuous variables
    # Note to use the knots from the training data 
    if(x$argumentList$useRCS) {
      nk <- ifelse(!is.null(x$argumentList$knotarg), x$argumentList$knotarg, 4)
      
      pred_cont <- lapply(var_num, function(i)
      {
        ## extract knots from original imputed data
        z <- rcs(data[[i]],nk)
        z_new <- rcs(pred_dat[[i]], attributes(z)$parms)
        
        z <- matrix(nrow = nrow(z_new), ncol = ncol(z_new))
        for(j in 1:ncol(z_new)){
          z[,j] <- z_new[,j]
        }
        z <- data.frame(z)
        # z=data.frame(as.numeric(z_new[,1]), as.numeric(z_new[,2]), check.names = F)
        znames <- vector(length = ncol(z))
        for (j in 1:ncol(z)){
          if (j==1) {
            znames[j] <-  i
          } else {
            ############################
            # stringr
            # znames[j] <-  paste0(i, str_c(rep(".p", times = j-1), sep = "", collapse = ""))
            # stringr
            ############################
            # Create a vector of strings
            vec <- paste0(rep(".p", times = j-1), collapse = "")
            # Combine i and vec
            znames[j] <- paste0(i, vec)
            ###############################
          }
        }
        colnames(z) <- znames
        z
      })
      
      pred_cont <- do.call(cbind, pred_cont)
      
      
      #########################################
      # dplyr
      # t=pred_cont %>% bind_cols(pred_dat %>% dplyr::select(-var_num))
      # dplyr
      ##########################################
      # Select all columns from pred_dat except var_num
      
      non_num <- pred_dat[, -which(names(pred_dat) %in% var_num)]
      # Combine pred_cont and selected_cols
      t <- cbind(pred_cont, non_num)
      ########################################
    } else t <- data.frame(pred_dat)
    
    # t$y= sample(data[[x$argumentList$survVar]], nrow(t), replace = T) ## add a pseudo y response
    
    t <- t[, union(setdiff(colnames(t), "TRT"), "TRT"), drop = F]  ### put TRT to the last column, so the model.matrix output will have the same order of labeling the interaction term: e.g. CD8:TRT instead of TRT:CD8
    
    # rownames(t)=1:nrow(t)
    
    for ( i in colnames(ref_c))
      t[[i]] <- factor(t[[i]], levels=levels(ref_c[[i]]))
    
    if(x$argumentList$withInteractions){
      
      # Create all 2-way interactions between aVar and bVar
      #####################################################
      # # dplyr
      # a_df <- dplyr::select(t, var_ab[1], str_subset(colnames(t), var_ab[2]))
      # a_df <- mutate_at(a_df, vars(contains(var_ab[2])), ~.*a_df[,1])
      # a_df <- rename_at(a_df, vars(contains(var_ab[2])), ~ paste0(var_ab[1], "_x_", .))
      # a_df <- dplyr::select(a_df, -var_ab[1])
      # # dplyr
      ####################################################
      # Select columns
      a_df <- t[, c(var_ab[1], grep(var_ab[2], colnames(t), value = TRUE)), drop = F]
      # Multiply columns
      cols_to_multiply <- grep(var_ab[2], colnames(a_df), value = TRUE)
      a_df[, cols_to_multiply] <- a_df[, cols_to_multiply, drop = F] * a_df[, var_ab[1]]
      # Rename columns
      cols_to_rename <- grep(var_ab[2], colnames(a_df), value = TRUE)
      new_col_names <- paste0(var_ab[1], "_x_", cols_to_rename)
      colnames(a_df) <- c(var_ab[1], new_col_names)
      # Remove column
      a_df <- a_df[, -which(names(a_df) == var_ab[1]), drop = F]
      #######################################################
      # # dplyr
      # b_df <- dplyr::select(t, var_ab[2], str_subset(colnames(t), var_ab[1]))
      # b_df <- mutate_at(b_df, vars(contains(var_ab[1])), ~.*b_df[,1])
      # b_df <- rename_at(b_df, vars(contains(var_ab[1])), ~ paste0(., "_x_", var_ab[2]))
      # b_df <- dplyr::select(b_df, -colnames(b_df[,1:2]))
      # # dplyr
      ##########################################################
      # Select columns
      b_df <- t[, c(var_ab[2], grep(var_ab[1], colnames(t), value = TRUE)), drop = F]
      # Multiply columns
      cols_to_multiply <- grep(var_ab[1], colnames(b_df), value = TRUE)
      b_df[, cols_to_multiply] <- b_df[, cols_to_multiply, drop = F] * b_df[, var_ab[2]]
      # Rename columns
      cols_to_rename <- grep(var_ab[1], colnames(b_df), value = TRUE)
      new_col_names <- paste0(cols_to_rename, "_x_", var_ab[2])
      colnames(b_df) <- c(var_ab[2], new_col_names)
      # Remove column
      b_df <- b_df[, -which(names(b_df) == var_ab[2]), drop = F]
      b_df <- b_df[, -which(names(b_df) == paste0(var_ab[1], "_x_", var_ab[2])), drop = F]
      #############################################################
      # # dplyr
      # t <- bind_cols(t, a_df)
      # t <- bind_cols(t, b_df)
      # # dplyr
      #########################################################
      t <- cbind(t, a_df)
      t <- cbind(t, b_df)
    }
    
    
    
    
    # Rename TRT as actual treatment variable name
    colnames(t)[colnames(t) == "TRT"] <- var_trt
    
    
   
  # Make Predictions
  pred_Full <- data.frame(predict(x$fitFull, t, conf.int = 0.95))
  colnames(pred_Full) <- paste0(colnames(pred_Full), "_predFull")
  
  pred_A <- data.frame(predict(x$fitA, t, conf.int = 0.95))
  colnames(pred_A) <- paste0(colnames(pred_A), "_predA")
  
  pred_B <- data.frame(predict(x$fitB, t, conf.int = 0.95))
  colnames(pred_B) <- paste0(colnames(pred_B), "_predB")
  
  ################################
  # # dplyr
  # preds <- bind_cols(pred_Full, pred_A, pred_B)
  # # dplyr
  #################################
  preds <- cbind(pred_Full, pred_A, pred_B)  
  ##################################
  
  name_preds <- colnames(preds)
  #################################
  # # dplyr
  # t <- bind_cols(t, preds)
  # # dplyr
  #################################
  t <- cbind(t, preds)
  ##################################
  
  
  #######################################
  # # dplyr
  # preds_cont <- filter(t, biom %in% var_num)
  # # dplyr
  ######################################
  preds_cont <- t[t$biom %in% var_num, , drop = F]
  ####################################
  
  if(x$argumentList$withInteractions){
    ##################################
    # # dplyr
    # preds_cont <- mutate(preds_cont, id = ifelse(biom %in% var_ab,
    #                                              c(1:(nintp*narms*5)),
    #                                              c(rep(1:(nintp*narms), times = length(var_cov)))))
    # # dplyr
    #################################
    preds_cont$id <- ifelse(preds_cont$biom %in% var_ab, 
                            c(1:(nintp*narms*5)),
                            c(rep(1:(nintp*narms), times = length(var_cov))))
    #####################################
  } else {
    ######################################
    # # dplyr
    # preds_cont <- mutate(preds_cont, id = c(rep(1:(nintp*narms), times = length(var_num))))
    # # dplyr
    #######################################
    preds_cont$id <- c(rep(1:(nintp*narms), times = length(var_num)))
    ###################################
  }
  
  ###########################################
  # # dplyr
  # preds_cont <- mutate(preds_cont, biom_id = str_c(biom, id))
  # 
  # preds_cont <- dplyr::select(preds_cont, name_preds, biom_id)
  # # dplyr
  ###########################################
  preds_cont$biom_id <- paste0(preds_cont$biom, preds_cont$id)
  
  preds_cont <- preds_cont[, c(name_preds, "biom_id"), drop = F]
  ############################################
  
  # Set up data for plotting
######################################  
  # # dplyr
  # plot_dat_cont <- dplyr::select(t, var_num, var_trt, biom)
  # plot_dat_cont <- dplyr::filter(plot_dat_cont, biom %in% var_num )
  # # dplyr
  ###################################
  plot_dat_cont <- t[,c(var_num, var_trt, "biom"), drop = F]
  plot_dat_cont <- plot_dat_cont[plot_dat_cont$biom %in% var_num, , drop = F]
  #####################################
  
  if(x$argumentList$withInteractions){
  ######################################
    # # dplyr
    #   plot_dat_cont <- mutate(plot_dat_cont, id = ifelse(biom %in% var_ab, 
    #                                                    c(1:(nintp*narms*5)),
    #                                                    c(rep(1:(nintp*narms), times = length(var_cov)))))
    # # dplyr
    #################################
      plot_dat_cont$id <- ifelse(plot_dat_cont$biom %in% var_ab,
                                 c(1:(nintp*narms*5)),
                                 c(rep(1:(nintp*narms), times = length(var_cov))))
  } else {
    ###############################
    # # dplyr
    # plot_dat_cont <- mutate(plot_dat_cont, id = c(rep(1:(nintp*narms), times = length(var_num))))
    # # dplyr
    ##############################
    plot_dat_cont$id <- c(rep(1:(nintp*narms), times = length(var_num)))
    #################################
  }
  
  ###################################
  # # dplyr
  # plot_dat_cont <- dplyr::mutate(plot_dat_cont, biom_id = stringr::str_c(biom, id))
  # # dplyr
  ###################################
  plot_dat_cont$biom_id <- paste0(plot_dat_cont$biom, plot_dat_cont$id)
  ###################################
  
  
  bioms <- list()
  
  if(x$argumentList$withInteractions){
    bioms_ab <- list()
    
    for (i in var_num) {
      if(i %in% var_ab) {
        #################################
        # # dplyr
        # bioms_ab[[i]] <- filter(plot_dat_cont, biom == i) 
        # bioms_ab[[i]] <- bioms_ab[[i]][,c(var_ab, var_trt)]
        # bioms_ab[[i]] <- mutate(bioms_ab[[i]], id = 1:nrow(bioms_ab[[i]]))
        # bioms_ab[[i]] <- mutate(bioms_ab[[i]], biom_id = str_c(i, id))
        # # dplyr
        #################################
        bioms_ab[[i]] <- plot_dat_cont[plot_dat_cont$biom == i,] 
        bioms_ab[[i]] <- bioms_ab[[i]][,c(var_ab, var_trt)]
        bioms_ab[[i]]$id <- 1:nrow(bioms_ab[[i]])
        bioms_ab[[i]]$biom_id <- paste0(i, bioms_ab[[i]]$id)
        #################################
      } else {
        ################################
        # # dplyr
        # bioms[[i]] <- filter(plot_dat_cont, biom == i) 
        # bioms[[i]] <- bioms[[i]][,c(i, var_trt)]
        # bioms[[i]] <- mutate(bioms[[i]], id = 1:nrow(bioms[[i]]))
        # # bioms[[i]] <- mutate(bioms[[i]], biom_id = str_c(i, id))
        # # dplyr
        ####################################
        bioms[[i]] <- plot_dat_cont[plot_dat_cont$biom == i,]
        bioms[[i]] <- bioms[[i]][,c(i, var_trt)]
        bioms[[i]]$id <- 1:nrow(bioms[[i]])
        ###################################
      }
    }
    
    if(sum(bioms_ab[[var_ab[1]]][[var_trt]] != bioms_ab[[var_ab[2]]][[var_trt]]) > 0){
      stop("Predictions not alligned correctly. Bug in code.")
    }
    ################################
    # # dplyr
    # if(!is_empty(var_cov)){
    #   # dplyr
    ##################################
      if(length(var_cov) > 0){
    #################################
        
      for(i in length(var_cov)) {
        if(i == 1){
          if(sum(bioms[[var_cov[i]]][[var_trt]] != bioms[[var_cov[i+1]]][[var_trt]]) > 0){
            stop("Predictions not alligned correctly. Bug in code.")
          }
        } else {
          if(sum(bioms[[var_cov[i]]][[var_trt]] != bioms[[var_cov[i-1]]][[var_trt]]) >0){
            stop("Predictions not alligned correctly. Bug in code.")
          }
        }
      }
    }
    
    plot_a <- bioms_ab[[var_ab[1]]]
    
    
    ########################################
    # # dplyr
    # plot_a <- dplyr::select(plot_a, -id)
    # plot_a <- tidyr::pivot_longer(plot_a, -c(var_trt, var_ab[2], biom_id), names_to = "Biomarker", values_to = "Value")
    # # dplyr
    #########################################
    plot_a <- plot_a[, -which(names(plot_a) == "id")]
    plot_a <- reshape(plot_a, varying = var_ab[1], idvar = "biom_id",
                      v.names = "Value", timevar = "Biomarker",
                      times = var_ab[1], direction = "long")
    ########################################
    
    plot_b <- bioms_ab[[var_ab[2]]]
    
    #########################################
    # # dplyr
    # plot_b <- dplyr::select(plot_b, -id)
    # plot_b <- pivot_longer(plot_b, -c(var_trt, var_ab[1], biom_id), names_to = "Biomarker", values_to = "Value")
    # # dplyr
    ##########################################
    plot_b <- plot_b[, -which(names(plot_b) == "id")]
    plot_b <- reshape(plot_b, varying = var_ab[2], idvar = "biom_id",
                      v.names = "Value", timevar = "Biomarker",
                      times = var_ab[2], direction = "long")
    #########################################
    
    # Merge with Predictions
    
    ###################################
    # # dplyr
    # plot_a <- inner_join(plot_a, preds_cont)
    # plot_b <- inner_join(plot_b, preds_cont)
    # # dplyr
    ##########################################################
    plot_a <- merge(plot_a, preds_cont, by = "biom_id", all = FALSE)
    plot_b <- merge(plot_b, preds_cont, by = "biom_id", all = FALSE)
    #########################################################
    plot_a[[var_trt]] <- as.factor(as.character(plot_a[[var_trt]]))
    plot_a[[var_ab[2]]] <- factor(as.character(round(plot_a[[var_ab[2]]], 3)), levels = c(round(quants[[var_ab[2]]], 3)), labels = c("10%", "25%", "50%", "75%", "90%"))
    # levels(plot_a[[var_ab[2]]]) <- paste(sort(as.numeric(levels(plot_a[[var_ab[2]]]))))
    
    plot_b[[var_trt]] <- as.factor(as.character(plot_b[[var_trt]]))
    # plot_b[[var_ab[1]]] <- as.factor(as.character(round(plot_b[[var_ab[1]]], 3)))
    # levels(plot_b[[var_ab[1]]]) <- paste(sort(as.numeric(levels(plot_b[[var_ab[1]]]))))
    plot_b[[var_ab[1]]] <- factor(as.character(round(plot_b[[var_ab[1]]], 3)), levels = c(round(quants[[var_ab[1]]], 3)), labels = c("10%", "25%", "50%", "75%", "90%"))
    
    #######################################
    # # dplyr
    # if(!is_empty(var_cov)){
    #   for (i in 2:length(bioms)){
    #     bioms[[i]] <- dplyr::select(bioms[[i]], -var_trt)
    #   }
    #   plot_dat_red <- Reduce(full_join, bioms)
    #   plot_dat_red <- pivot_longer(plot_dat_red, -c(var_trt, id), names_to = "Biomarker", values_to = "Value")
    #   plot_dat_red <- mutate(plot_dat_red, biom_id = str_c(Biomarker, id))
    #   
    #   # Merge Predictions with plotting data
    #   plot_dat_red <- inner_join(plot_dat_red, preds_cont)
    #   
    #   
    #   # plot_dat_cont <- dplyr::select(plot_dat_red, biom_id, name_preds)
    #   
    #   
    #   plot_dat <- arrange(plot_dat_red, Biomarker, var_trt, Value)
    # # dplyr
      ##############################
      if(length(var_cov) > 0) {
        for (i in 2:length(bioms)){
          bioms[[i]] <- bioms[[i]][, -which(names(bioms[[i]]) == var_trt)]
        }
        
        plot_dat_red <- bioms[[1]]
        
        for(i in 2:length(bioms)){
          plot_dat_red <- merge(plot_dat_red, bioms[[i]], all = TRUE)
        }
       plot_dat_red <- reshape(plot_dat_red, varying = var_cov, v.names = "Value",
                               timevar = "Biomarker", times = var_cov, direction = "long")
       plot_dat_red$biom_id <- paste0(plot_dat_red$Biomarker, plot_dat_red$id)
       
       # Merge Predictions with plotting data
       plot_dat_red <- merge(plot_dat_red, preds_cont, by = "biom_id", all = F)
       
       plot_dat <- plot_dat_red[order(plot_dat_red$Biomarker, plot_dat_red[, var_trt], plot_dat_red$Value), ]
      ##############################
     
      plot_dat[[var_trt]] <- as.factor(as.character(plot_dat[[var_trt]]))
    }
  } else {
    for (i in var_num) {
      
      ##################################
      # # dplyr
      # bioms[[i]] <- filter(plot_dat_cont, biom == i) 
      # bioms[[i]] <- bioms[[i]][,c(i, var_trt)]
      # bioms[[i]] <- mutate(bioms[[i]], id = 1:nrow(bioms[[i]]))  
      # # dplyr
      #####################################
      bioms[[i]] <- plot_dat_cont[plot_dat_cont$biom == i,] 
      bioms[[i]] <- bioms[[i]][,c(i, var_trt)]
      bioms[[i]]$id <- 1:nrow(bioms[[i]])
      #####################################
    }
    
    for (i in 1:length(bioms)){
      if(i ==1){
        
        if(sum(bioms[[i]][[var_trt]] != bioms[[i+1]][[var_trt]]) > 0){
          stop("Predictions not alligned correctly. Bug in code.")
        }
      } else {
        if(sum(bioms[[i]][[var_trt]] != bioms[[i-1]][[var_trt]]) >0){
          stop("Predictions not alligned correctly. Bug in code.")
        }
      }
    }
    
    
    
      ########################################
    #   # dplyr
    # for (i in 2:length(bioms)){
    #   bioms[[i]] <- dplyr::select(bioms[[i]], -var_trt)
    # }
    #   plot_dat_red <- Reduce(full_join, bioms)
    #   plot_dat_red <- pivot_longer(plot_dat_red, -c(var_trt, id), names_to = "Biomarker", values_to = "Value")
    #   plot_dat_red <- mutate(plot_dat_red, biom_id = str_c(Biomarker, id))
    #   
    #   # Merge Predictions with plotting data
    #   plot_dat_red <- inner_join(plot_dat_red, preds_cont)
    #   
    #   
    #   # plot_dat_cont <- dplyr::select(plot_dat_red, biom_id, name_preds)
    #   
    #   
    #   plot_dat <- arrange(plot_dat_red, Biomarker, var_trt, Value)
    #   # dplyr
      #########################################
      for (i in 2:length(bioms)){
      bioms[[i]] <- bioms[[i]][, -which(names(bioms[[i]]) == var_trt)]
      }
      
      plot_dat_red <- bioms[[1]]
      
      for(i in 2:length(bioms)){
        plot_dat_red <- merge(plot_dat_red, bioms[[i]], all = TRUE)
      }
      plot_dat_red <- reshape(plot_dat_red, varying = var_num, v.names = "Value",
                              timevar = "Biomarker", times = var_num, direction = "long")
      plot_dat_red$biom_id <- paste0(plot_dat_red$Biomarker, plot_dat_red$id)
      
      # Merge Predictions with plotting data
      plot_dat_red <- merge(plot_dat_red, preds_cont, by = "biom_id", all = F)
      
      plot_dat <- plot_dat_red[order(plot_dat_red$Biomarker, plot_dat_red[, var_trt], plot_dat_red$Value), ]
      ##########################################
    
      plot_dat[[var_trt]] <- as.factor(as.character(plot_dat[[var_trt]]))
    
  }
  
  
  ###############################################
  # # dplyr
  # if(is_empty(x$argumentList$trtVar)){
  # # dplyr
  ############################################
    if(length(x$argumentList$trtVar) == 0) {
  ###########################################
    if(x$argumentList$withInteractions) {
      
      # Make partial effects plots
      PE_FullA <- ggplot(plot_a, aes(x = Value, y = linear.predictors_predFull)) +
        # geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull), alpha = 0.1, linetype = "dotted") +
        # scale_color_manual(values = c("red", "blue")) +
        geom_line() +
        # ylim(-3, 2) +
        facet_grid(cols = vars(get(var_ab[2]))) +
        ylab("Relative log Hazard") +
        xlab(var_ab[1]) +
        ggtitle(paste0("Full Model Fit: Faceted by ", var_ab[2])) +
        theme_bw(18)
      
      
      
      PE_FullB <- ggplot(plot_b, aes(x = Value, y = linear.predictors_predFull)) +
        # geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull), alpha = 0.1, linetype = "dotted") +
        # scale_color_manual(values = c("red", "blue")) +
        geom_line() +
        # ylim(-3, 2) +
        facet_grid(cols = vars(get(var_ab[1]))) +
        ylab("Relative log Hazard") +
        xlab(var_ab[2]) +
        ggtitle(paste0("Full Model Fit: Faceted by ", var_ab[1])) +
        theme_bw(18)
      
      
      
      PE_A <- ggplot(plot_a, aes(x = Value, y = linear.predictors_predA)) +
        # geom_ribbon(aes(ymin = lower_predA, ymax = upper_predA), alpha = 0.1, linetype = "dotted") +
        # scale_color_manual(values = c("red", "blue")) +
        geom_line() +
        # ylim(-3, 2) +
        # facet_wrap(~Biomarker) +
        ylab("Relative log Hazard") +
        xlab(var_ab[1]) +
        ggtitle(paste0("Model Fit: ", var_ab[1])) +
        theme_bw(18)
      
      
      
      PE_B <- ggplot(plot_b, aes(x = Value, y = linear.predictors_predB)) +
        # geom_ribbon(aes(ymin = lower_predB, ymax = upper_predB), alpha = 0.1, linetype = "dotted") +
        # scale_color_manual(values = c("red", "blue")) +
        geom_line() +
        # ylim(-3, 2) +
        # facet_wrap(~Biomarker) +
        ylab("Relative log Hazard") +
        xlab(var_ab[2]) +
        ggtitle(paste0("Model Fit: ", var_ab[2])) +
        theme_bw(18)
      
      if(!is.null(ylimA)){
        PE_FullA <- PE_FullA + coord_cartesian(ylim = ylimA)
      }  
      
      if(!is.null(ylimB)){
        PE_FullB <- PE_FullB + coord_cartesian(ylim = ylimB)
      }   
      
      if(confInt) {
        PE_FullA <- PE_FullA +
          geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull), alpha = 0.1, linetype = "dotted")
        PE_FullB <- PE_FullB +
          geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull), alpha = 0.1, linetype = "dotted")
        PE_A <- PE_A +
          geom_ribbon(aes(ymin = lower_predA, ymax = upper_predA), alpha = 0.1, linetype = "dotted")
        PE_B <- PE_B +
          geom_ribbon(aes(ymin = lower_predB, ymax = upper_predB), alpha = 0.1, linetype = "dotted")
      }
      
      #################################
      # # dplyr
      # if(!is_empty(var_cov)){
      #   # dplyr
      ###################################
        if(length(var_cov) > 0) {
      ###################################
        # Make partial effects plots
        PE_Full_cov <- ggplot(plot_dat, aes(x = Value, y = linear.predictors_predFull)) +
          # geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull), alpha = 0.1, linetype = "dotted") +
          # scale_color_manual(values = c("red", "blue")) +
          geom_line() +
          # ylim(-3, 2) +
          facet_wrap(~Biomarker) +
          ylab("Relative log Hazard") +
          xlab("Predictor Value") +
          ggtitle("Full Model Fit") +
          theme_bw(18)
        
        
        
        PE_A_cov <- ggplot(plot_dat, aes(x = Value, y = linear.predictors_predA)) +
          # geom_ribbon(aes(ymin = lower_predA, ymax = upper_predA), alpha = 0.1, linetype = "dotted") +
          # scale_color_manual(values = c("red", "blue")) +
          geom_line() +
          # ylim(-3, 2) +
          facet_wrap(~Biomarker) +
          ylab("Relative log Hazard") +
          xlab("Predictor Value") +
          ggtitle(paste0("Model Fit: ", var_ab[1])) +
          theme_bw(18)
        
        
        
        PE_B_cov <- ggplot(plot_dat, aes(x = Value, y = linear.predictors_predB)) +
          # geom_ribbon(aes(ymin = lower_predB, ymax = upper_predB), alpha = 0.1, linetype = "dotted") +
          # scale_color_manual(values = c("red", "blue")) +
          geom_line() +
          # ylim(-3, 2) +
          facet_wrap(~Biomarker) +
          ylab("Relative log Hazard") +
          xlab("Predictor Value") +
          ggtitle(paste0("Model Fit: ", var_ab[2])) +
          theme_bw(18)
        
        
        
        if(confInt) {
          PE_Full_cov <- PE_Full_cov +
            geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull), alpha = 0.1, linetype = "dotted")
          PE_A_cov <- PE_A_cov +
            geom_ribbon(aes(ymin = lower_predA, ymax = upper_predA), alpha = 0.1, linetype = "dotted")
          PE_B_cov <- PE_B_cov +
            geom_ribbon(aes(ymin = lower_predB, ymax = upper_predB), alpha = 0.1, linetype = "dotted")
        }
        
        PE_Full <- list(A = PE_FullA, B = PE_FullB, Cov = PE_Full_cov)
        PE_A <- list(A = PE_A, Cov = PE_A_cov)
        PE_B <- list(B = PE_B, Cov = PE_B_cov)
      } else {
        
        
        PE_Full <- list(A = PE_FullA, B = PE_FullB)
        PE_A <- list(A = PE_A)
        PE_B <- list(B = PE_B)
      }
      
      
      
      
    } else {
      
      # Make partial effects plots
      PE_Full <- ggplot(plot_dat, aes(x = Value, y = linear.predictors_predFull)) +
        # geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull), alpha = 0.1, linetype = "dotted") +
        # scale_color_manual(values = c("red", "blue")) +
        geom_line() +
        # ylim(-3, 2) +
        facet_wrap(~Biomarker) +
        ylab("Relative log Hazard") +
        xlab("Predictor Value") +
        ggtitle("Full Model Fit") +
        theme_bw(18)
      
      
      PE_A <- ggplot(plot_dat, aes(x = Value, y = linear.predictors_predA)) +
        # geom_ribbon(aes(ymin = lower_predA, ymax = upper_predA), alpha = 0.1, linetype = "dotted") +
        # scale_color_manual(values = c("red", "blue")) +
        geom_line() +
        # ylim(-3, 2) +
        facet_wrap(~Biomarker) +
        ylab("Relative log Hazard") +
        xlab("Predictor Value") +
        ggtitle(paste0("Model Fit: ", var_ab[1])) +
        theme_bw(18)
      
      PE_B <- ggplot(plot_dat, aes(x = Value, y = linear.predictors_predB)) +
        # geom_ribbon(aes(ymin = lower_predB, ymax = upper_predB), alpha = 0.1, linetype = "dotted") +
        # scale_color_manual(values = c("red", "blue")) +
        geom_line() +
        # ylim(-3, 2) +
        facet_wrap(~Biomarker) +
        ylab("Relative log Hazard") +
        xlab("Predictor Value") +
        ggtitle(paste0("Model Fit: ", var_ab[2])) +
        theme_bw(18)
    
      if(confInt) {
        PE_Full <- PE_Full +
          geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull), alpha = 0.1, linetype = "dotted")
        PE_A <- PE_A +
          geom_ribbon(aes(ymin = lower_predA, ymax = upper_predA), alpha = 0.1, linetype = "dotted")
        PE_B <- PE_B +
          geom_ribbon(aes(ymin = lower_predB, ymax = upper_predB), alpha = 0.1, linetype = "dotted")
      }
    }
  } else {
    if(x$argumentList$withInteractions) {
      
      # Make partial effects plots
      PE_FullA <- ggplot(plot_a, aes(x = Value, y = linear.predictors_predFull, color = get(var_trt), linetype = get(var_trt))) +
        # geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull, group = get(var_trt)), alpha = 0.1, linetype = "dotted") +
        # scale_color_manual(values = c("red", "blue")) +
        geom_line() +
        # ylim(-3, 2) +
        facet_grid(cols = vars(get(var_ab[2]))) +
        ylab("Relative log Hazard") +
        xlab(var_ab[1]) +
        ggtitle(paste0("Full Model Fit: Faceted by ", var_ab[2])) +
        theme_bw(18)
      
     
      
      PE_FullB <- ggplot(plot_b, aes(x = Value, y = linear.predictors_predFull, color = get(var_trt), linetype = get(var_trt))) +
        # geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull, group = get(var_trt)), alpha = 0.1, linetype = "dotted") +
        # scale_color_manual(values = c("red", "blue")) +
        geom_line() +
        # ylim(-3, 2) +
        facet_grid(cols = vars(get(var_ab[1]))) +
        ylab("Relative log Hazard") +
        xlab(var_ab[2]) +
        ggtitle(paste0("Full Model Fit: Faceted by ", var_ab[1])) +
        theme_bw(18)
      
      
      PE_A <- ggplot(plot_a, aes(x = Value, y = linear.predictors_predA, color = get(var_trt), linetype = get(var_trt))) +
        # geom_ribbon(aes(ymin = lower_predA, ymax = upper_predA, group = get(var_trt)), alpha = 0.1, linetype = "dotted") +
        # scale_color_manual(values = c("red", "blue")) +
        geom_line() +
        # ylim(-3, 2) +
        # facet_wrap(~Biomarker) +
        ylab("Relative log Hazard") +
        xlab(var_ab[1]) +
        ggtitle(paste0("Model Fit: ", var_ab[1])) +
        theme_bw(18)
      
      PE_B <- ggplot(plot_b, aes(x = Value, y = linear.predictors_predB, color = get(var_trt), linetype = get(var_trt))) +
        # geom_ribbon(aes(ymin = lower_predB, ymax = upper_predB, group = get(var_trt)), alpha = 0.1, linetype = "dotted") +
        # scale_color_manual(values = c("red", "blue")) +
        geom_line() +
        # ylim(-3, 2) +
        # facet_wrap(~Biomarker) +
        ylab("Relative log Hazard") +
        xlab(var_ab[2]) +
        ggtitle(paste0("Model Fit: ", var_ab[2])) +
        theme_bw(18)
      
      if(!is.null(ylimA)){
        PE_FullA <- PE_FullA + coord_cartesian(ylim = ylimA)
      }  
      
      if(!is.null(ylimB)){
        PE_FullB <- PE_FullB + coord_cartesian(ylim = ylimB)
      }   
      
      if(confInt) {
        PE_FullA <- PE_FullA +
          geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull, group = get(var_trt)), alpha = 0.1, linetype = "dotted")
        PE_FullB <- PE_FullB +
          geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull, group = get(var_trt)), alpha = 0.1, linetype = "dotted")
        PE_A <- PE_A +
          geom_ribbon(aes(ymin = lower_predA, ymax = upper_predA, group = get(var_trt)), alpha = 0.1, linetype = "dotted")
        PE_B <- PE_B +
          geom_ribbon(aes(ymin = lower_predB, ymax = upper_predB, group = get(var_trt)), alpha = 0.1, linetype = "dotted")
      }
      
      #################################
      # # dplyr
      # if(!is_empty(var_cov)){
      #   # dplyr
      ###############################
        if(length(var_cov) > 0) {
      #############################
        # Make partial effects plots
        PE_Full_cov <- ggplot(plot_dat, aes(x = Value, y = linear.predictors_predFull, color = get(var_trt), linetype = get(var_trt))) +
          # geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull, group = get(var_trt)), alpha = 0.1, linetype = "dotted") +
          # scale_color_manual(values = c("red", "blue")) +
          geom_line() +
          # ylim(-3, 2) +
          facet_wrap(~Biomarker) +
          ylab("Relative log Hazard") +
          xlab("Predictor Value") +
          ggtitle("Full Model Fit") +
          theme_bw(18)
        
        PE_A_cov <- ggplot(plot_dat, aes(x = Value, y = linear.predictors_predA, color = get(var_trt), linetype = get(var_trt))) +
          # geom_ribbon(aes(ymin = lower_predA, ymax = upper_predA, group = get(var_trt)), alpha = 0.1, linetype = "dotted") +
          # scale_color_manual(values = c("red", "blue")) +
          geom_line() +
          # ylim(-3, 2) +
          facet_wrap(~Biomarker) +
          ylab("Relative log Hazard") +
          xlab("Predictor Value") +
          ggtitle(paste0("Model Fit: ", var_ab[1])) +
          theme_bw(18)
        
        PE_B_cov <- ggplot(plot_dat, aes(x = Value, y = linear.predictors_predB, color = get(var_trt), linetype = get(var_trt))) +
          # geom_ribbon(aes(ymin = lower_predB, ymax = upper_predB, group = get(var_trt)), alpha = 0.1, linetype = "dotted") +
          # scale_color_manual(values = c("red", "blue")) +
          geom_line() +
          # ylim(-3, 2) +
          facet_wrap(~Biomarker) +
          ylab("Relative log Hazard") +
          xlab("Predictor Value") +
          ggtitle(paste0("Model Fit: ", var_ab[2])) +
          theme_bw(18)
       
        
        if(confInt) {
          PE_Full_cov <- PE_Full_cov +
            geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull, group = get(var_trt)), alpha = 0.1, linetype = "dotted")
          PE_A_cov <- PE_A_cov +
            geom_ribbon(aes(ymin = lower_predA, ymax = upper_predA, group = get(var_trt)), alpha = 0.1, linetype = "dotted")
          PE_B_cov <- PE_B_cov +
            geom_ribbon(aes(ymin = lower_predB, ymax = upper_predB, group = get(var_trt)), alpha = 0.1, linetype = "dotted")
        }
        
        PE_Full <- list(A = PE_FullA, B = PE_FullB, Cov = PE_Full_cov)
        PE_A <- list(A = PE_A, Cov = PE_A_cov)
        PE_B <- list(B = PE_B, Cov = PE_B_cov)
      } else {
        
        PE_Full <- list(A = PE_FullA, B = PE_FullB)
        PE_A <- list(A = PE_A)
        PE_B <- list(B = PE_B)
      }
      
      
      
      
    } else {
      
      # Make partial effects plots
      PE_Full <- ggplot(plot_dat, aes(x = Value, y = linear.predictors_predFull, color = get(var_trt), linetype = get(var_trt))) +
        # geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull, group = get(var_trt)), alpha = 0.1, linetype = "dotted") +
        # scale_color_manual(values = c("red", "blue")) +
        geom_line() +
        # ylim(-3, 2) +
        facet_wrap(~Biomarker) +
        ylab("Relative log Hazard") +
        xlab("Predictor Value") +
        ggtitle("Full Model Fit") +
        theme_bw(18)
      
      PE_A <- ggplot(plot_dat, aes(x = Value, y = linear.predictors_predA, color = get(var_trt), linetype = get(var_trt))) +
        # geom_ribbon(aes(ymin = lower_predA, ymax = upper_predA, group = get(var_trt)), alpha = 0.1, linetype = "dotted") +
        # scale_color_manual(values = c("red", "blue")) +
        geom_line() +
        # ylim(-3, 2) +
        facet_wrap(~Biomarker) +
        ylab("Relative log Hazard") +
        xlab("Predictor Value") +
        ggtitle(paste0("Model Fit: ", var_ab[1])) +
        theme_bw(18)
      
      PE_B <- ggplot(plot_dat, aes(x = Value, y = linear.predictors_predB, color = get(var_trt), linetype = get(var_trt))) +
        # geom_ribbon(aes(ymin = lower_predB, ymax = upper_predB, group = get(var_trt)), alpha = 0.1, linetype = "dotted") +
        # scale_color_manual(values = c("red", "blue")) +
        geom_line() +
        # ylim(-3, 2) +
        facet_wrap(~Biomarker) +
        ylab("Relative log Hazard") +
        xlab("Predictor Value") +
        ggtitle(paste0("Model Fit: ", var_ab[2])) +
        theme_bw(18)
      
      if(confInt) {
        PE_Full <- PE_Full +
          geom_ribbon(aes(ymin = lower_predFull, ymax = upper_predFull, group = get(var_trt)), alpha = 0.1, linetype = "dotted")
        PE_A <- PE_A +
          geom_ribbon(aes(ymin = lower_predA, ymax = upper_predA, group = get(var_trt)), alpha = 0.1, linetype = "dotted")
        PE_B <- PE_B +
          geom_ribbon(aes(ymin = lower_predB, ymax = upper_predB, group = get(var_trt)), alpha = 0.1, linetype = "dotted")
      }
    }
    
  }
  
  
  
  
  
  out <- list(Full = PE_Full, FitA = PE_A, FitB = PE_B)  
  out
  
}

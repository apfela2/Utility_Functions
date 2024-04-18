## Abraham Apfel
## December 28, 2018


# Function to output maineffects coefficients table from glinternet model

# Table of main effects coefficients mapped to their corresponding gene

# Depends on knitr, dplyr

# coef.out is output from coef function on single lambda value from glinternet output
# X is matrix of predictors supplied in glinternet model
# cat_var is vector of categorical variables in X
# time_effect if used discrete method with time_effect set to T


mainTable <- function(coef.out, X, cat_var = NULL, time_effect = F) {
  
  if (!is.null(cat_var)) {
    main_cat_index <- unlist(coef.out[[1]][["mainEffects"]][["cat"]])
    main_cont_index <- unlist(coef.out[[1]][["mainEffects"]][["cont"]])
    cont_X <- dplyr::select(as.data.frame(X), -cat_var)
    cat_X <- dplyr::select(as.data.frame(X), cat_var)
    mainGene_cat <- colnames(cat_X[main_cat_index])
    cat_levels <- max(as.data.frame(X)[cat_var]) + 1
    
    if(time_effect) {
      cont_X$Time <- NA
    }
    tmp_list <- vector("list", length(mainGene_cat))
    for (i in 1:length(mainGene_cat)) {
      tmp_list[[i]] <- rep(mainGene_cat[i], cat_levels[i])
      for (j in 1:length(tmp_list[[i]])) {
        tmp_list[[i]][j] <- paste0(tmp_list[[i]][j], j - 1)
      }
    }
    mainGene_cat <- do.call("c", tmp_list)
    
    mainGene_cont <- colnames(cont_X[main_cont_index])
    mainCoef_cat <- unlist(coef.out[[1]][["mainEffectsCoef"]][["cat"]])
    mainCoef_cont <- unlist(coef.out[[1]][["mainEffectsCoef"]][["cont"]])
    main_cat <- cbind(Biomarker = mainGene_cat, Coefficients = round(mainCoef_cat, 3))
    main_cont <- cbind(Biomarker = mainGene_cont, Coefficients = round(mainCoef_cont, 3))
    
    main <- rbind(main_cat, main_cont)
  } else {
    main_cont_index <- unlist(coef.out[[1]][["mainEffects"]][["cont"]])
    cont_X <- as.data.frame(X)
    mainGene_cont <- colnames(cont_X[main_cont_index])
    mainCoef_cont <- unlist(coef.out[[1]][["mainEffectsCoef"]][["cont"]])
    main_cont <- cbind(Biomarker = mainGene_cont, Coefficients = round(mainCoef_cont, 3))
  }
  
  # knitr::kable(main) 
}

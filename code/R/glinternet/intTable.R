## Abraham Apfel
## December 29, 2018


# Function to output interactions coefficients table from glinternet model

# Table of interaction coefficients mapped to their corresponding biomarkers

# Depends on knitr, dplyr

# coef.out is output from coef function on single lambda value from glinternet output
# X is matrix of predictors supplied in glinternet model
# cat_var is vector of categorical variables in X


intTable <- function(coef.out, X, cat_var = NULL) {
  
  if (!is.null(cat_var)) {
    cont_X <- dplyr::select(as.data.frame(X), -cat_var)
    cat_X <- dplyr::select(as.data.frame(X), cat_var)
    cat_levels <- max(as.data.frame(X)[cat_var]) + 1
    
    # Table of interactions
    int_cont_index <- unlist(coef.out[[1]][["interactions"]][["contcont"]])
    int_cat_index <- unlist(coef.out[[1]][["interactions"]][["catcont"]])
    
    intGene_cat <- rep(NA, nrow(int_cat_index))
    for (i in 1:nrow(int_cat_index)) {
      intGene_cat[i] <- colnames(cat_X[int_cat_index[i,1]])  
    }
    
    
    tmp_list <- vector("list", length(intGene_cat))
    for (i in 1:length(intGene_cat)) {
      tmp_list[[i]] <- rep(intGene_cat[i], cat_levels[int_cat_index[i, 1]])
      for (j in 1:length(tmp_list[[i]])) {
        tmp_list[[i]][j] <- paste0(tmp_list[[i]][j], j - 1)
      }
    }
    intGene_cat <- do.call("c", tmp_list)
    
    intGene_cat_cont <- rep(NA, nrow(int_cat_index))
    for (i in 1:nrow(int_cat_index)) {
      intGene_cat_cont[i] <- colnames(cont_X[int_cat_index[i,2]])  
    }
    tmp_list <- vector("list", length(intGene_cat_cont))
    for (i in 1:length(intGene_cat_cont)) {
      tmp_list[[i]] <- rep(intGene_cat_cont[i], cat_levels[int_cat_index[i, 1]])
    }
    intGene_cat_cont <- do.call("c", tmp_list)
    
    intCoef_cat <- unlist(coef.out[[1]][["interactionsCoef"]][["catcont"]])
    int_cat <- cbind(Biomarker1 = intGene_cat, Biomarker2 = intGene_cat_cont, Coefficients = round(intCoef_cat, 3))
    
    if(!is.null(int_cont_index)){
      intGene_cont1 <- rep(NA, nrow(int_cont_index))
      for (i in 1:nrow(int_cont_index)) {
        intGene_cont1[i] <- colnames(cont_X[int_cont_index[i,1]])  
      }
      intGene_cont2 <- rep(NA, nrow(int_cont_index))
      for (i in 1:nrow(int_cont_index)) {
        intGene_cont2[i] <- colnames(cont_X[int_cont_index[i,2]])  
      }
      
      intCoef_cont <- unlist(coef.out[[1]][["interactionsCoef"]][["contcont"]])
      int_cont <- cbind(Biomarker1 = intGene_cont1, Biomarker2 = intGene_cont2, Coefficients = round(intCoef_cont, 3))
      
      int <- rbind(int_cat, int_cont)
    } else {
      int <- int_cat
    }
  } else {
    cont_X <- as.data.frame(X)
    int_cont_index <- unlist(coef.out[[1]][["interactions"]][["contcont"]])
    intGene_cont1 <- rep(NA, nrow(int_cont_index))
    
    for (i in 1:nrow(int_cont_index)) {
      intGene_cont1[i] <- colnames(cont_X[int_cont_index[i,1]])  
    }
    intGene_cont2 <- rep(NA, nrow(int_cont_index))
    for (i in 1:nrow(int_cont_index)) {
      intGene_cont2[i] <- colnames(cont_X[int_cont_index[i,2]])  
    }
    
    intCoef_cont <- unlist(coef.out[[1]][["interactionsCoef"]][["contcont"]])
    int_cont <- cbind(Biomarker1 = intGene_cont1, Biomarker2 = intGene_cont2, Coefficients = round(intCoef_cont, 3))
  }
  
  # kable(int)
}  

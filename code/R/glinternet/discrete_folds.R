# Helper function to make CV on discrete data based off of folds defined by continuous data

# folds.list is the column with fold ID's in discrete data format
discrete_folds <- function(folds.list) {
  num.folds <- length(unique(folds.list))
  
  
  idxs.list <- list()
  # names(idxs.list) <- paste0(rep(gsub("  ", "_", names(folds.list)), each = num.folds), 
  #                            "_Fold_", 
  #                            rep(seq_len(num.folds), times = 1))
  
  
  for (f.idx in sort(unique(folds.list))) {
    tn              <- paste0("Run_1", "_Fold_", f.idx)
    idxs.list[[tn]] <- which(folds.list != f.idx)
  }
  
  
  return(idxs.list)
}
# plotHists
# gghist #function for making a gg plot histrogram
# make.multiplots
# multiplot
# fitted.vs.actual #returns ggplot object of fitted vs actural for a regression line



#some plotting functions
plotHists <- function(tmp, show.plots = FALSE, plots.fname = NULL, inc = 11, nc = 3) {	
  if (class(tmp) == "matrix") {
    cat("Input dataset must be a data frame, converting")
    tmp = data.frame(tmp, stringsAsFactors= FALSE)
  }
  
  plots = vector(mode = "list", length = ncol(tmp)); names(plots) = colnames(tmp)
  for (tv in names(tmp)) {
    plots[[tv]] = gghist(tmp[, tv], vname = tv, nbins = NULL, show.plot = FALSE) 
    #  ggplot(data = tmp) + 
    #  geom_histogram(aes_string(x = tv)) +
    #  theme_bw()
  }
  
  if (show.plots == TRUE)
    return(plots)
  
  make.multiplots(plots, inc = inc, nc = nc, out.file = plots.fname, silent = TRUE)
}

#function for making a gg plot histrogram
gghist = function(x, vname = "Variable", nbins = NULL, show.plot = TRUE, larger.font = FALSE) {
  if (is.numeric(x))
    x = data.frame(x)
  
  if (is.character(x))
    x = data.frame(x, stringsAsFactors = TRUE)
  
  if (ncol(x) > 1) {
    warning("this function makes a histogram for only one variable at time")
    colnames(x)[1] = "x"
  }
  
  if (is.null(nbins))
    nbins = floor(sqrt(nrow(x))*1.5) #multiplying by 1.5 to get slightly more bins
  
  title.size = 10
  axis.size = 8
  
  if (larger.font == TRUE) {
    title.size = 24
    axis.size = 20
  }
  
  if (is.numeric(x[, 1])) {
    p  = ggplot(data = x, aes(x = x)) + 
      geom_histogram(bins = nbins) +
      scale_x_continuous() + 
      xlab(vname) + 
      ggtitle(paste0("Distribution of ", vname)) + 
      theme_bw() +
      theme(plot.title = element_text(size = title.size),
            axis.title = element_text(size = axis.size))
    
  } else {
    p = ggplot(data = x, aes(x = x)) + 
      geom_bar() +
      xlab(vname) + 
      ggtitle(paste0("Distribution of ", vname)) +
      theme_bw() + 
      theme(plot.title = element_text(size = title.size),
            axis.title = element_text(size = axis.size))
    
  }
  
  if (show.plot == FALSE)
    return(p)
  
  print(p)
}


make.multiplots <- function(plots, inc = 11, nc = 3, out.file = NULL, silent = TRUE) {
  if (!is.null(out.file))
    pdf(out.file, width = 20, height =20) 
  
  plotted = rep(FALSE, length = length(plots))
  
  p.start = 1
  p.end   = min(p.start + inc, length(plots))
  
  while(any(!plotted)) {
    if (is.null(out.file) & silent == FALSE)
      quartz()
    cat(".")
    multiplot(plotlist = plots[p.start:p.end], cols = nc)
    plotted[p.start:p.end] = TRUE
    p.start = p.end + 1
    p.end   = min(p.start + inc, length(plots)) 
  }
  
  if (!is.null(out.file))
    dev.off()
}

multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#returns ggplot object of fitted vs actural for a regression line
fitted.vs.actual <- function(mod.fit, title = "Fitted Vs Actual") {
  yhat = fitted(mod.fit)
  y = mod.fit$model[, all.vars(lm.res$call)[1]]
  
  dat = data.frame(y = y, yhat = yhat)
  
  return(ggplot(data = dat, aes(x = yhat, y = y)) + geom_point())
}
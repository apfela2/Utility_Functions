## Abraham Apfel
## July 18, 2018


# Functions


# Function to transform continuous survival into discrete
disc <- function(dat, cens, surv.time, intervals = min(dat[[ surv.time ]]) : ((max(dat[[ surv.time ]])) + 1)) {
  
  # cens is binary censoring variable where 0 is censored and 1 is observed
  # surv.time is name of column in dat with continuous survival times
  # intervals is a vector of right interval borders for which the discrete survival times should be broken into. By default they increase by increments of 1 from minimum observed survival time to maximum
  
  #  colDat <- dat[[ surv.time ]]
  # intervals<-min(colDat) : (max(colDat)+1)
  dat <- dat[which(!is.na(dat[[surv.time]])),]
  
  if (any(dat[[ surv.time ]] < 0)) {
    stop("surival times must be non-negative - call to disc() with surv.time=", surv.time)
  }
  
  disc<-discSurv::contToDisc(dat, surv.time ,intervals)
  
  # Transform data to discrete time format
  disc_dat<-discSurv::dataLong(disc,"timeDisc",censColumn = cens)
}

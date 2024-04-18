cph2predLRtest <- function( aVar, bVar, survVar, censorVar, data,
   trtVar = character(0), covars = character(0), useRCS = TRUE, knotArg = 4,
   withInteractions = TRUE, keepData = TRUE, ... ) {
#
# NOTES:
  # useRCS: consider "a", "b", c("a", "b"), character(), instead of logical?
  #         Or, "a", "b", "both", "none"?
  #         Or, list(a = 4, b = 3), list(), list(b = NA),...,and remove knotArg?
  #         Or, c(a = 4, b = 3), c(), c(b = 3), ...?
# DATE CREATED:
#    2023-05-10 
#
# DESCRIPTION:
#    Fits nested Cox PH models to compare, via likelihood-ratio tests, a model
#    with two explanatory variables, A and B, to the reduced models with only
#    A or only B, and to compare each of those models to a "null" model with
#    neither A nor B.
#
#    Can optionally specify:
#    (1) to use restricted cubic splines for both A and B (currently, can use
#    rcs for only both or neither A and B). You can specify the number of knots
#    if you like but it must be the same for A and B.
#    (2) to include interactions between A and B in the full model (using %ia%
#    to exclude doubly-nonlinear interaction terms if using rcs for A and B).
#    Currently no ":" option to include doubly-nonlinear interaction terms. 
#    (3) to include a treatment-arm variable in the models, in which case
#    interactions between treatment arm and both A and B also will be included.
#    (Currently, you can NOT exclude the treatment interactions.)
#    (4) to include a set of covariates in the model. Currently these can be
#    included only additively; no interactions with A, B, treatment arm, or
#    each other are allowed.
#
# ARGUMENTS:
# aVar  character string giving the name of the column in data to use as
#    variable "A".
# bVar  character string giving the name of the column in data to use as
#    variable "B".
# survVar  character string giving the name of the column in data to use as the
#    time-to-event outcome.
# censorVar  character string giving the name of the column in data to use as
#    the censoring variable. NOTE: The censoring variable MUST contain only
#    values 1 (censored) and 0 (not censored)! BE SURE YOU ARE USING A
#    CENSORING INDICATOR HERE, NOT AN EVENT INDICATOR!!
# data  dataframe that includes at least all the necessary variables. Missing
#    values are NOT allowed in any of the variables used for fitting the models.
#    NOTE: give the dataframe name, NOT a character string.
# trtVar  character string giving the name of the column in data to use as the
#    treatment-arm variable. Should be a factor. The default, character(0),
#    specifies NO treatment arm in the models.
# covars  character vector giving the name(s) of column(s) in data to use as
#    covariates to include additively in the models. The default, character(0),
#    specifies NO covariates.
# useRCS  logical flag. If TRUE (the default), both A and B will be modeled
#    using restricted cubic splines.
# knotArg  numeric vector specifying the number of knots, or the knot locations,
#    to use if useRCS is TRUE. The default is 4. Note that this differs from
#    the ordinary default for rcs. You can, if you really want, use NA to tell
#    this function to use the usual default number of knots (typically 5, if
#    not overriden by options).
# withInteractions  logical flag. If TRUE (the default), the full model
#    includes interactions between A and B. NOTE: if length(trtVar) > 0 - that
#    is, you specify to include a treatment arm variable in the models - then
#    all interactions between treatment arm and both A and B
#    will be included in the models. THIS argument controls only whether or not
#    to include interactions between A and B.
# keepData  logical flag. If TRUE (the default), the object passed as the value
#    for argument "data" is included in the argumentList component of the
#    returned object. Otherwise NA is returned for that component. The default
#    is useful for e.g., getting plot methods to work without extra arguments.
#    But because the data object might be very large, we give you the option
#    to NOT store it as part of the returned object.
# ...  additional optional arguments to pass to function cph when fitting the
#    Cox PH models.
#
# VALUE:
#    a list with class "cph2pred" and with components:
# summary  a 5 by 3 numeric matrix giving the key likelihood-ratio (LR) test
#    results collected from the next 5 components. Just for convenience.
# testAgivenB  an "lrtest" object giving the results of a LR test for the
#    overall effect of variable A, given that variable B is in the model.
# testBgivenA  an "lrtest" object giving the results of a LR test for the
#    overall effect of variable B, given that variable A is in the model.
# testAandB  an "lrtest" object, or (if no covariates and no treatment arm
#    in the model) a list with single component "stats", giving the results
#    of a LR test for the overall effect of both variables A and B.
# testA  an "lrtest" object, or (if no covariates and no treatment arm
#    in the model) a list with single component "stats", giving the results
#    of a LR test for the overall effect of variable A (ignoring variable B).
# testB  an "lrtest" object, or (if no covariates and no treatment arm
#    in the model) a list with single component "stats", giving the results
#    of a LR test for the overall effect of variable B (ignoring variable A).
# formFull  formula constructed to fit the full model
# formA  formula constructed to fit the reduced model that excludes variable B
# formB  formula constructed to fit the reduced model that excludes variable A
# formNULL  formula constructed to fit the reduced (NULL) model that excludes
#    both variables A and B
# fitFull  a "cph" object giving the fitted full model
# fitA  a "cph" object giving the fitted reduced model that excludes variable B
# fitB  a "cph" object giving the fitted reduced model that excludes variable A
# fitNULL  a "cph" object giving the fitted reduced model that excludes both
#   variables A and B
# argumentList  a list giving the values of all the arguments in the call to
# this function, whether passed explicitly or by default. But see the keepData
#   argument for one caveat.
# call  a "call" object giving the current call to this function, with any
#   abbreviated argument names expanded to their full names.
#
#
# SEE ALSO:
#   rms::rcs, rms::rcspline.eval, rms::cph, rms::lrtest
#
# EXAMPLES:
# set.seed(699)
# n <- 1000
# dat <- data.frame( r = rnorm(n), theta = rnorm(n), OS = rexp(n),
#    OSC = rbinom(n, 1, 0.25), ACTARM = as.factor(rep(c("C", "T"), each = n/2)),
#    X1 = rnorm(n), X2 = as.factor(rbinom(n, 1, 0.5)))
###
# Single-arm trial (or ignore treatment arm):
###
# cph2predLRtest( aVar = "r", bVar = "theta", survVar = "OS", censorVar = "OSC",
#    data = dat )
# cph2predLRtest( aVar = "r", bVar = "theta", survVar = "OS", censorVar = "OSC",
#    data = dat, useRCS = FALSE )
# cph2predLRtest( aVar = "r", bVar = "theta", survVar = "OS", censorVar = "OSC",
#    data = dat, knotArg = 5 )
# cph2predLRtest( aVar = "r", bVar = "theta", survVar = "OS", censorVar = "OSC",
#    data = dat, knotArg = 5, withInteractions = FALSE )
# cph2predLRtest( aVar = "r", bVar = "theta", survVar = "OS", censorVar = "OSC",
#    data = dat, method = "efron" )
# cph2predLRtest( aVar = "r", bVar = "theta", survVar = "OS", censorVar = "OSC",
#    data = dat, covars = c("X1", "X2") )
# cph2predLRtest( aVar = "r", bVar = "theta", survVar = "OS", censorVar = "OSC",
#    data = dat, covars = c("rcs(X1, 4)", "X2") )
###
# Two-arm trial:
###
# cph2predLRtest( aVar = "r", bVar = "theta", survVar = "OS", censorVar = "OSC",
#    data = dat, trtVar = "ACTARM" )
# cph2predLRtest( aVar = "r", bVar = "theta", survVar = "OS", censorVar = "OSC",
#    data = dat, trtVar = "ACTARM", withInteractions = FALSE )
# cph2predLRtest( aVar = "r", bVar = "theta", survVar = "OS", censorVar = "OSC",
#    data = dat, trtVar = "ACTARM", useRCS = FALSE )
# cph2predLRtest( aVar = "r", bVar = "theta", survVar = "OS", censorVar = "OSC",
#    data = dat, trtVar = "ACTARM", withInteractions = FALSE, useRCS = FALSE )
###
# THIS CURRENTLY NOT ALLOWED. CANNOT USE ANY TRANSFORMATIONS, SUCH AS rcs, OF
# COVARIATES WHEN GIVE A VALUE FOR trtVar (that is, for multi-arm trials):
###
# cph2predLRtest( aVar = "r", bVar = "theta", survVar = "OS", censorVar = "OSC",
#    data = dat, trtVar = "ACTARM", covars = c("rcs(X1, 4)", "X2") )
###
# BUT THIS SHOULD WORK:
###
# cph2predLRtest( aVar = "r", bVar = "theta", survVar = "OS", censorVar = "OSC",
#    data = dat, trtVar = "ACTARM", covars = c("X1", "X2") )
#
###
# Check that all required var names are in data and none have missing values.
###
cnames <- colnames(data)
if (!any(cnames == aVar))
   stop("Value for argument aVar not found in colnames(data).")
if (!any(cnames == bVar))
   stop("Value for argument bVar not found in colnames(data).")
if (!any(cnames == survVar))
   stop("Value for argument survVar not found in colnames(data).")
if (!any(cnames == censorVar))
   stop("Value for argument censorVar not found in colnames(data).")
usecovars <- (length(covars) > 0)
usetrt <- (length(trtVar) > 0)
# Skip checking for covar names so that can allow transformations (e.g., rcs)
# of covars. For now, these are allowed only when length(trtVar) == 0.
# if (usecovars && !all(is.element(covars, cnames)))
#    stop("Value(s) for argument covars not found in colnames(data).")
if (usetrt && !any(cnames == trtVar))
   stop("Value for argument trtVar not found in colnames(data).")
# Also skip checking for missing values in covars. If present, the cph calls
# later will fail.
# thevars <- c(aVar, bVar, survVar, censorVar, trtVar, covars)
thevars <- c(aVar, bVar, survVar, censorVar, trtVar)
if (any(is.na(data[, thevars, drop = FALSE])))
   stop("Missing values not allowed in variables included in models.")
###
# Check that censorVar contains only values 0 and 1.
###
if (!all(is.element(data[[censorVar]], 0:1)))
   stop("data[[censorVar]] must contain only values 0 and 1.")
###
# Build left-hand-side for formula.
###
lhs <- paste( "Surv(", survVar, ", 1 -", censorVar, ")" )
###
# Add in rcs if desired. Currently only allows for both a and b or neither.
# I leave here 2 separate calls to is.na, even though ever so slightly
# inefficient, to make it easier to specify separate knotArg values for
# variables a and b, should we wish to do so later.
###
aVar.orig <- aVar
bVar.orig <- bVar
if (useRCS) {
   aVar <- paste0( "rcs( ", aVar )
   if (!is.na(knotArg)) aVar <- paste( aVar, ",", knotArg )
   aVar <- paste( aVar, ")" )
   bVar <- paste0( "rcs( ", bVar )
   if (!is.na(knotArg)) bVar <- paste( bVar, ",", knotArg )
   bVar <- paste( bVar, ")" )
}
###
# Build formula for full model and fit model.
###
#
# Options to explore when useRCS = TRUE, withInteractions = TRUE, and
# usetrt = TRUE.
# If useRCS = FALSE, no issue; all the 2-way and the 3-way (if withInteractions
# = TRUE) interactions would be fine.
# If withInteractions = FALSE, no issue; only 2-way interactions and they would
# be fine, with or without useRCS = TRUE.
# So, only a problem when both useRCS and withInteractions are TRUE! The
# problem centers on having a 3-way interaction in which 2 of the variables use
# rcs. Only affects the FULL model with more than one trt arm, since that is
# the only model here that can have a 3-way interaction.
# (1) trt + rcs(r):trt + rcs(theta):trt + rcs(r) %ia% rcs(theta) %ia% trt
# (2) trt + rcs(r):trt + rcs(theta):trt + rcs(r) %ia% rcs(theta) : trt
# (3) trt + rcs(r):trt + rcs(theta):trt + r:theta:trt
#
# (1) will not work; to quote Frank Harrell here:
# https://stats.stackexchange.com/questions/541916/\
# meaning-of-interaction-with-ia-in-rms-three-way-interaction
# "...%ia% does not extend to three-way interactions."
# (2) will not work; I tried and cph fails with a somewhat cryptic complaint
# (3) should work but is a reduced model compared to what I would LIKE to fit.
#
# Can I solve this using model.matrix (e.g. method = "model.matrix" in initial
# call to cph for full model) to build exactly the model matrix I want?
# Probably. But can I do this with a reasonable amount of programming effort?
# OK, yes, for some, possibly dubious, definitions of "reasonable".
#
rhsFull <- c( aVar, bVar )
if (withInteractions) {
    if (useRCS) {
       rhsFull <- c( rhsFull, paste( aVar, "%ia%", bVar ) )
    } else {
       rhsFull <- c( rhsFull, paste0( aVar, ":", bVar ) )
    }
}
if (!useRCS || !withInteractions || !usetrt) {
   if (usetrt) {
      trtia <- paste( trtVar, rhsFull, sep = ":" )
      rhsFull <- c( trtVar, rhsFull, trtia )
   }
   if (usecovars) rhsFull <- c( covars, rhsFull )
   rhsFull <- paste( rhsFull, collapse = " + " )
   formFull <- as.formula(paste( lhs, "~", rhsFull ))
   fitFull <- cph( formFull, data = data, na.action = na.fail, ... )
} else {
   # useRCS and withInteractions and usetrt all TRUE; need special approach for
   # 3-way interx.
   # Need initial single-arm formula here, without covars; then call cph with
   # method = "model.matrix".
   rhsFull <- paste( rhsFull, collapse = " + " )
   formFull <- as.formula(paste( lhs, "~", rhsFull ))
   mm <- cph( formFull, data = data, na.action = na.fail, method =
      "model.matrix" )
   # mm has the basis functions for A and B and the not-doubly-nonlinear
   # interaction variables.
   rhsFull <- colnames(mm)
   # Rename mm columns so will work in a formula
   rhsFull <- gsub("\'", ".p", rhsFull)
   rhsFull <- gsub(" \\* ", "_x_", rhsFull)
   colnames(mm) <- rhsFull
   trtia <- paste( trtVar, rhsFull, sep = ":" )
   rhsFull <- c( trtVar,  rhsFull, trtia )
   if (usecovars) {
      rhsFull <- c( covars, rhsFull )
      dat2 <- data[, c(survVar, censorVar, trtVar, covars), drop = FALSE]
   }
   else {
      dat2 <- data[, c(survVar, censorVar, trtVar), drop = FALSE]
   }
   rhsFull <- paste( rhsFull, collapse = " + " )
   formFull <- as.formula(paste( lhs, "~", rhsFull ))
   mm <- cbind(mm, dat2)
   fitFull <- cph( formFull, data = mm, na.action = na.fail, ... )
}
###
# Build formula for a-only model and fit model.
###
rhsA <- aVar
if (usetrt) {
    trtia <- paste( trtVar, rhsA, sep = ":" )
    rhsA <- c( trtVar, rhsA, trtia )
}
if (usecovars) rhsA <- c( covars, rhsA )
rhsA <- paste( rhsA, collapse = " + " )
formA <- as.formula(paste( lhs, "~", rhsA ))
fitA <- cph( formA, data = data, na.action = na.fail, ... )
###
# Build formula for b-only model and fit model.
###
rhsB <- bVar
if (usetrt) {
    trtia <- paste( trtVar, rhsB, sep = ":" )
    rhsB <- c( trtVar, rhsB, trtia )
}
if (usecovars) rhsB <- c( covars, rhsB )
rhsB <- paste( rhsB, collapse = " + " )
formB <- as.formula(paste( lhs, "~", rhsB ))
fitB <- cph( formB, data = data, na.action = na.fail, ... )
###
# Build formula for null model and fit model.
###
if (usetrt) {
   if (usecovars) {
      rhsNULL <- c( covars, trtVar )
      rhsNULL <- paste( rhsNULL, collapse = " + " )
   } else {
      rhsNULL <- trtVar
   }
} else {
   if (usecovars) {
      rhsNULL <- paste( covars, collapse = " + " )
   } else {
      # NOTE: coxph and anova would work with this but cph does not return
      # a useful result that could be passed to lrtest.
      # Will need a work-around for this special case.
      rhsNULL <- "1"
   }
} 
formNULL <- as.formula(paste( lhs, "~", rhsNULL ))
fitNULL <- cph( formNULL, data = data, na.action = na.fail, ... )
###
# Do the desired likelihood-ratio tests.
###
testAgivenB <- lrtest( fitFull, fitB )
testBgivenA <- lrtest( fitFull, fitA )
if (usetrt || usecovars) {
   testAandB <- lrtest( fitFull, fitNULL )
   testA <- lrtest( fitA, fitNULL )
   testB <- lrtest( fitB, fitNULL )
} else {
   # The work-around:
   testAandB <- list( stats = fitFull$stats[ c("Model L.R.", "d.f.", "P") ] )
   testA <- list( stats = fitA$stats[ c("Model L.R.", "d.f.", "P") ] )
   testB <- list( stats = fitB$stats[ c("Model L.R.", "d.f.", "P") ] )
   names(testAandB$stats)[1] <- "L.R. Chisq"
   names(testA$stats)[1] <- "L.R. Chisq"
   names(testB$stats)[1] <- "L.R. Chisq"
}
###
# Construct summary table.
# Perhaps this should be a separate "summary" function, but convenient for
# me to put it here for now...
###
statsMat <- array( NA, c(5, 3), list( c(
   "Overall Effect of A | B",
   "Overall Effect of B | A",
   "Overall Effect of both A and B",
   "Overall Effect of A",
   "Overall Effect of B"),
   c("L.R. Chisq", "d.f.", "P") ) )
statsMat[ "Overall Effect of A | B", ] <- testAgivenB$stats
statsMat[ "Overall Effect of B | A", ] <- testBgivenA$stats
statsMat[ "Overall Effect of both A and B", ] <- testAandB$stats
statsMat[ "Overall Effect of A", ] <- testA$stats
statsMat[ "Overall Effect of B", ] <- testB$stats
###
# Assemble argument list.
###
argumentList <- list(
   aVar = aVar.orig,
   bVar = bVar.orig,
   survVar = survVar,
   censorVar = censorVar,
   data = if (keepData) data else NA,
   trtVar = trtVar,
   covars = covars,
   useRCS = useRCS,
   knotArg = knotArg,
   withInteractions = withInteractions,
   keepData = keepData,
   ...) 
###
# Assemble return list and go home.
###
out <- list(
      summary = statsMat,
      testAgivenB = testAgivenB,
      testBgivenA = testBgivenA,
      testAandB = testAandB,
      testA = testA,
      testB = testB,
      formFull = formFull,
      formA = formA,
      formB = formB,
      formNULL = formNULL,
      fitFull = fitFull,
      fitA = fitA,
      fitB = fitB,
      fitNULL = fitNULL,
      argumentList = argumentList,
      call = match.call() )
oldClass(out) <- "cph2pred"
out
}

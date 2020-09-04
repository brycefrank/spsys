#' Constructs a function that produces the GREG estimate
#' of every variable defined in the mapping. The mapping defines
#' the covariates that will be used for each variable.
#' 
#' @param mapping a named list of attribute - covariate pairs such that
#' each named element is a vector of covariate names that will be used to fit
#' the greg for the named element.
#' @return a function that can be called on a SysFrame object that returns
#' a GREGOut object.
greg <- function(mapping) {
  # TODO look into heteroskedasticity (probably can use nlme::gls)
    # would require some expansion of the mapping...potentially pass some sort of wrapped nlme::gls object?
  
  greg_estimator <- function(sys_frame, pop) {
    # TODO a little redundant to have the mapping as a named list since the response variables are
    # defined on the left hand side, but this will do for now
    mods  <- list()
    estim <- list() # TODO this should just be similar to the VarOUT stuff, check there.
    
    n <- nrow(sys_frame@data)
    N <- nrow(pop)
    p <- length(names(mapping))
    
    resids <- data.frame(matrix(nrow=n, ncol=p))
    colnames(resids) <- names(mapping)
    
    for(att in names(mapping)) {
      form <- mapping[[att]]
      covt_names <- labels(terms(form))
      
      data_cols <- colnames(sys_frame@data)[colnames(sys_frame@data) %in% covt_names]
      data_cols <- c(att, data_cols)
      lm_data <- sys_frame@data[, data_cols]
      
      # Divide every column by the vector of inclusion probabilities
      # TODO double check that this is correct
      lm_data <- sweep(lm_data, 2, sys_frame@pi, FUN='/')
      
      mod <- lm(form, data=lm_data)
      mods[[att]] <- mod
      
      # Make the estimation of the mean for each attribute
      pop_mod <- pop[,data_cols]
      pop_mod <- sweep(pop_mod, 1, pop$pi, FUN='/')
      pop_pred <- predict(mod, newdata=pop_mod)
      mu <- sum(pop_pred) / N
      estim[[att]] <- mu
      
      # Aggregate the residuals into a data frame with attribute names
      resids[,att] <- resid(mod)
      
      browser()
      greg_out <- GREGOut(mods, estim, resids)
      return(greg_out)
    }
  }
  return(greg_estimator) 
}

setClass('GREGOut', slots = list(
  mods = 'list',
  est = 'list', # TODO make this a dataframe
  resid = 'data.frame'
))


GREGOut <- function(mods, est, resid) {
  greg_out <- new('GREGOut')
  greg_out@mods <- mods
  greg_out@est <- est
  greg_out@resid <- resid
  
  greg_out
}
#' Maximum-likelihood fitting of the CTP distribution
#'
#' @description
#' Maximum-likelihood fitting of the Complex Triparametric Pearson (CTP) distribution with parameters \eqn{a}, \eqn{b} and \eqn{\gamma}. Generic
#' methods are \code{print}, \code{summary}, \code{coef}, \code{logLik}, \code{AIC}, \code{BIC} and \code{plot}. 
#'
#' @usage
#' fitctp(x, astart = NULL, bstart = NULL, gammastart = NULL, 
#'           method = "L-BFGS-B", control = list(), ...)
#'        
#' @param x A numeric vector of length at least one containing only finite values.
#' @param astart A starting value for the parameter \eqn{a>0}; by default NULL.
#' @param bstart A starting value for the parameter \eqn{b}; by default NULL.
#' @param gammastart A starting value for the parameter \eqn{\gamma>max(0,2a)}; by default NULL.
#' @param method The method to be used in fitting the model. See 'Details'. 
#' @param control A list of parameters for controlling the fitting process.
#' @param ...  Additional parameters.
#'
#' @details If the starting values of the parameters \eqn{a}, \eqn{b} and \eqn{\gamma} are omitted (default option), 
#' they are computing by the method of moments (if possible; otherwise they must be entered).
#' 
#' The default method is \code{"L-BFGS-B"} (see details in \code{\link{optim}} function),
#' but non-linear minimization (\code{\link{nlm}}) or those included in the \code{optim} function (\code{"Nelder-Mead"}, 
#' \code{"BFGS"}, \code{"CG"} and \code{"SANN"}) may be used.
#' 
#' Standard error (SE) estimates for \eqn{a}, \eqn{b}
#' and \eqn{\gamma} are provided by the default method; otherwise, SE for \eqn{\gamma_0} where \eqn{\gamma=exp(\gamma_0)} is computed.
#'
#' @return An object of class \code{'fitCTP'} is a list containing the following components:
#'
#' \itemize{
#' \item \code{n}, the number of observations,
#' \item \code{initialValues}, a vector with the starting values used,
#' \item \code{coefficients}, the parameter ML estimates of the CTP distribution,
#' \item \code{se}, a vector of the standard error estimates,
#' \item \code{hessian}, a symmetric matrix giving an estimate of the Hessian at the solution found in the optimization of the log-likelihood function,
#' \item \code{cov}, an estimate of the covariance matrix of the model coefficients,
#' \item \code{corr}, an estimate of the correlation matrix of the model estimates,
#' \item \code{loglik}, the maximized log-likelihood,
#' \item \code{aic}, Akaike Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters,
#' \item \code{bic}, Bayesian Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters,
#' \item \code{code}, a code that indicates successful convergence of the fitter function used (see nlm and optim helps),
#' \item \code{converged},  logical value that indicates if the optimization algorithms succesfull,
#' \item \code{method}, the name of the fitter function used.
#' }
#' 
#' Generic functions:
#'
#' \itemize{
#' \item \code{print}: The print of a \code{'fitCTP'} object shows the ML parameter estimates and their standard errors in parenthesis.
#' \item \code{summary}: The summary provides the ML parameter estimates, their standard errors and the statistic and p-value of the Wald test to check if the parameters are significant.
#' This summary also shows the loglikelihood, AIC and BIC values, as well as the results for the chi-squared goodness-of-fit test and the Kolmogorov-Smirnov test for discrete variables. Finally, the correlation matrix between parameter estimates appears.
#' \item \code{coef}: It extracts the fitted coefficients from a \code{'fitCTP'} object.
#' \item \code{logLik}: It extracts the estimated log-likelihood from a \code{'fitCTP'} object.
#' \item \code{AIC}: It extracts the value of the Akaike Information Criterion from a \code{'fitCTP'} object.
#' \item \code{BIC}: It extracts the value of the Bayesian Information Criterion from a \code{'fitCTP'} object.
#' \item \code{plot}: It shows the plot of a \code{'fitCTP'} object. Observed and theoretical probabilities, empirical and theoretical cumulative distribution functions or empirical cumulative probabilities against theoretical cumulative probabilities are the three plot types.
#' }
#' 
#' @importFrom stats nlm optim coef runif
#' @export
#' @references 
#' 
#' \insertRef{RCSO2004}{cpd}
#' 
#' \insertRef{ROC2018}{cpd}
#' 
#' @seealso
#' Plot of observed and theoretical frequencies for a CTP fit: \code{\link{plot.fitCTP}}
#' 
#' Maximum-likelihood fitting for the CBP distribution: \code{\link{fitcbp}}.
#' 
#' Maximum-likelihood fitting for the EBW distribution: \code{\link{fitebw}}.
#'
#' @examples
#' set.seed(123)
#' x <- rctp(500, -0.5, 1, 2)
#' fitctp(x)
#' summary(fitctp(x, astart = 1, bstart = 1.1, gammastart = 3))

fitctp <- function(x, astart = NULL, bstart = NULL, gammastart = NULL, 
                   method = "L-BFGS-B", control = list(), ...)
{

  #Checking
  if ( mode(x) != "numeric")
    stop( paste ("non-numeric argument to mathematical function", sep = ""))

  if( is.null(control$maxit) )
    control$maxit <- 10000

  if( is.null(control$trace) )
    control$trace <- 0
  if( is.null(control$warn) )
    control$warn <- -2
  #set warning level
  defWarn <- getOption("warn")
  options(warn=control$warn)
  
  if ( !is.numeric(control$maxit) || control$maxit <= 0 )
    stop( "maximum number of iterations must be > 0" )
  
  if ( is.numeric(astart) && is.numeric(bstart) && is.numeric(gammastart) && (!is.nan(astart+bstart+gammastart))){
    moments<-FALSE
  }else{
    moments<-TRUE
  }
  if (moments){
    #Try to generate initial values using the method of moments
    n <- length(x)
    m1 <- mean(x)
    m2 <- sum(x ^ 2) / n
    m3 <- sum(x ^ 3) / n

    #System equations
    A <- matrix(c(1, 1 + m1, 1 + 2 * m1 + m2, m1, m1 + m2, m1 + 2 * m2 + m3, -m1, -m2, -m3), ncol = 3, nrow = 3)
    b <- matrix(c(0, -m2, -2 * m3 - m2), ncol = 1, nrow = 3)

    #Solve the system
    sol <- try(solve(A)%*%b,silent = TRUE)
    if ('try-error' %in% class(sol)){
      stop("The method of moments does not provide any estimates. Enter initial values for the parameters.")
    }
     
    #New initial values
    astart <- sol[2] / 2
    bstart <- sqrt(sol[1] - astart ^ 2)
    gammastart <- sol[3] + 1
    
    if (is.nan(bstart) || gammastart <= max(0, 2 * astart)){
      stop("The method of moments does not provide any estimates. Enter initial values for the parameters.")
    }

  }
  #Imaginary unit
  i<-sqrt(as.complex(-1))

  #Log-likelihood
  if (method != "L-BFGS-B") {
    logL <- function(p){
      a <- p[1]
      b <- p[2]
      gama <- exp(p[3])
      respuesta <- -sum(2*Re(complex_gamma(gama-a+b*i,log=TRUE)+log(complex_gamma(a+b*i+x))-log(complex_gamma(a+b*i)))-lgamma(gama-2*a)-lgamma(gama+x))
      return(respuesta)
    }

    pstart <- c(astart,bstart,log(gammastart))

  }else {
    logL <- function(p){
      a <- p[1]
      b <- p[2]
      gama <- p[3]
      respuesta <- -sum(2*Re(complex_gamma(gama-a+b*i,log=TRUE)+log(complex_gamma(a+b*i+x))-log(complex_gamma(a+b*i)))-lgamma(gama-2*a)-lgamma(gama+x))
      return(respuesta)
      }

      pstart <- c(astart,bstart,gammastart)
    }

  #Optimization process
  converged<-FALSE
  
  if (method == "nlm") {
    methodText = "nlm"
    fit <- try(nlm(logL, p = pstart, hessian = TRUE, iterlim = control$maxit, print.level = control$trace))
    if ('try-error' %in% class(fit)){
      if (control$trace > 0) {
        cat("Crashed 'nlm' initial fit", "\n")
      }
    }
    else {
      fit$value <- fit$minimum
      fit$par <- fit$estimate
      fit$convergence <- fit$code
      if (fit$convergence < 3)
        converged = TRUE
    }
  }
  else if (any(method == c("Nelder-Mead", "BFGS", "CG", "SANN"))) {
    methodText <- method
    fit <- try(optim(pstart, logL, method = method, hessian = TRUE,
                 control = list(maxit = control$maxit, trace = control$trace)))
    if ('try-error' %in% class(fit)){
      if (control$trace > 0) {
        cat(paste("Crashed '",method,"' initial fit",sep=","),"\n")
      }
    }
    else {
      if (fit$convergence == 0)
        converged = TRUE
    }
  }
  else if (any(method == c("L-BFGS-B"))) {
    methodText <- method
    fit <- try(optim(pstart, logL, method = method, lower = c(-Inf,0,0.0000001), upper = c(Inf,Inf,Inf), hessian = TRUE, control = list(maxit = control$maxit, trace = control$trace)))
    if ('try-error' %in% class(fit)){
      if (control$trace > 0) {
        cat(paste("Crashed '",method,"' initial fit",sep=","),"\n")
      }
    }
    else {
      if (fit$convergence == 0)
        converged<-TRUE
    }
  }else{
    stop("Incorrect method")
  }
  
  if (converged){
    if (method != "L-BFGS-B") {
      coef.table <- rbind(fit$par, deparse.level = 0)
      dimnames(coef.table) <- list("", c("a", "b", "log(gamma)"))
    }
    else {
      coef.table <- rbind(fit$par, deparse.level = 0)
      dimnames(coef.table) <- list("", c("a", "b", "gamma"))
	}
    
    #Results
    results <- list(
      x = x,
      n = length(x),
      loglik = - (fit$value + sum(lfactorial(x))),
      aic = 2 * (fit$value + sum(lfactorial(x))) + 6,
      bic = 2 * (fit$value + sum(lfactorial(x))) + 2 * log(length(x)),
      coefficients =  coef.table,
      code = fit$convergence,
      hessian = fit$hessian,
      cov = solve(fit$hessian),
      se = sqrt(diag(solve(fit$hessian))),
      corr = solve(fit$hessian) / (sqrt(diag(solve(fit$hessian))) %o%
                                     sqrt(diag(solve(fit$hessian)))),
      code = fit$convergence,
      converged = converged,
      initialValues = c(astart, bstart, gammastart),
      method = methodText
    )
  } else{
    results<-list()
    results$method<-methodText
    results$converged<-FALSE
    results$n<-nrow(x)
  }
  
  #restore warning level
  options(warn=defWarn)
  class(results) <- "fitCTP"
  results
}


#' Maximum-likelihood fitting of the CBP distribution
#'
#' @description
#' Maximum-likelihood fitting of the Complex Biparametric Pearson (CBP) distribution with parameters \eqn{b} and \eqn{\gamma}. Generic
#' methods are \code{print}, \code{summary}, \code{coef}, \code{logLik}, \code{AIC}, \code{BIC} and \code{plot}. 
#'
#' @usage
#' fitcbp(x, bstart = NULL, gammastart = NULL, method = "L-BFGS-B", control = list(), ...)
#'        
#'  
#' @param x A numeric vector of length at least one containing only finite values.
#' @param bstart A starting value for the parameter \eqn{b}; by default NULL.
#' @param gammastart A starting value for the parameter \eqn{\gamma>0}; by default NULL.
#' @param method The method to be used in fitting the model. See 'Details'. 
#' @param control A list of parameters for controlling the fitting process.
#' @param ...  Additional parameters.
#
#' @details If the starting values of the parameters \eqn{b} and \eqn{\gamma} are omitted (default option), 
#' they are computing by the method of moments (if possible; otherwise they must be entered).
#' 
#' The default method is \code{"L-BFGS-B"} (see details in \code{\link{optim}} function),
#' but non-linear minimization (\code{\link{nlm}}) or those included in the \code{optim} function (\code{"Nelder-Mead"}, 
#' \code{"BFGS"}, \code{"CG"} and \code{"SANN"}) may be used.
#' 
#' Standard error (SE) estimates for \eqn{b} and \eqn{\gamma} are provided by the default method; 
#' otherwise, SE for \eqn{\gamma_0} where \eqn{\gamma=exp{(\gamma_0})} is computed.
#' 
#' @return An object of class \code{'fitCBP'} is a list containing the following components:
#'
#' \itemize{
#' \item \code{n}, the number of observations,
#' \item \code{initialValues}, a vector with the starting values used,
#' \item \code{coefficients}, the parameter ML estimates of the CTP distribution,
#' \item \code{se}, a vector of the standard error estimates,
#' \item \code{hessian}, a symmetric matrix giving an estimate of the Hessian at the solution found in the optimization of the log-likelihood function,
#' \item \code{cov}, an estimate of the covariance matrix of the model coefficients,
#' \item \code{corr}, an estimate of the correlation matrix of the model estimates,
#' \item \code{loglik}, the maximized log-likelihood,
#' \item \code{aic}, Akaike Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters,
#' \item \code{bic}, Bayesian Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters,
#' \item \code{code}, a code that indicates successful convergence of the fitter function used (see nlm and optim helps),
#' \item \code{converged},  logical value that indicates if the optimization algorithms succesfull,
#' \item \code{method}, the name of the fitter function used.
#' }
#' 
#' Generic functions:
#'
#' \itemize{
#' \item \code{print}: The print of a \code{'fitCBP'} object shows the ML parameter estimates and their standard errors in parenthesis.
#' \item \code{summary}: The summary provides the ML parameter estimates, their standard errors and the statistic and p-value of the Wald test to check if the parameters are significant.
#' This summary also shows the loglikelihood, AIC and BIC values, as well as the results for the chi-squared goodness-of-fit test and the Kolmogorov-Smirnov test for discrete variables. Finally, the correlation matrix between parameter estimates appears.
#' \item \code{coef}: It extracts the fitted coefficients from a \code{'fitCBP'} object.
#' \item \code{logLik}: It extracts the estimated log-likelihood from a \code{'fitCBP'} object.
#' \item \code{AIC}: It extracts the value of the Akaike Information Criterion from a \code{'fitCBP'} object.
#' \item \code{BIC}: It extracts the value of the Bayesian Information Criterion from a \code{'fitCBP'} object.
#' \item \code{plot}: It shows the plot of a \code{'fitCBP'} object. Observed and theoretical probabilities, empirical and theoretical cumulative distribution functions or empirical cumulative probabilities against theoretical cumulative probabilities are the three plot types.
#' }
#' 
#' @importFrom stats nlm optim coef runif
#' @export
#' 
#' @references 
#' 
#' \insertRef{RCS2003}{cpd}
#' 
#' @seealso
#' Plot of observed and theoretical frequencies for a CBP fit: \code{\link{plot.fitCBP}}
#' 
#' Maximum-likelihood fitting for the CTP distribution: \code{\link{fitctp}}.
#' 
#' Maximum-likelihood fitting for the EBW distribution: \code{\link{fitebw}}.
#'
#' @examples
#' set.seed(123)
#' x <- rcbp(500, 1.75, 3.5)
#' fitcbp(x)
#' summary(fitcbp(x, bstart = 1.1, gammastart = 3))

fitcbp <- function(x, bstart = NULL, gammastart = NULL, 
                      method = "L-BFGS-B", control = list(), ...)
{
  #Checking
  if ( mode(x) != "numeric")
    stop( paste ("non-numeric argument to mathematical function", sep = ""))

  if( is.null(control$maxit) )
    control$maxit <- 10000

  if( is.null(control$trace) )
    control$trace <- 0
  if( is.null(control$warn) )
    control$warn <- -2
  #set warning level
  defWarn <- getOption("warn")
  options(warn=control$warn)
  
  if ( !is.numeric(control$maxit) || control$maxit <= 0 )
    stop( "maximum number of iterations must be > 0" )
  
  if ( is.numeric(bstart) && is.numeric(gammastart) && (!is.nan(bstart+gammastart)))
    moments<-FALSE
  else
    moments<-TRUE
  
  if (moments){
    #Try to generate initial values using the method of moments
    n <- length(x)
    m1 <- mean(x)
    m2 <- sum(x ^ 2) / n
    if (m2 - m1 ^ 2 -m1 <= 0){
      stop("The method of moments does not provide any estimates. Enter initial values for the parameters.")
    }else{
      #New initial values
      bstart <- sqrt(m1 * m2 / (m2 - m1  ^ 2 - m1))
      gammastart <- bstart ^ 2 / m1 +1
    }
  }
  
  if (gammastart<=0)
    stop("gammastart must be positive.")
  
  #Imaginary unit
    i<-sqrt(as.complex(-1))

  #Log-likelihood
    if (method != "L-BFGS-B") {
    logL <- function(p){
      b <- p[1]
      gama <- exp(p[2])
      respuesta <- -sum(2*Re(complex_gamma(gama+b*i,log=TRUE)+log(complex_gamma(b*i+x))-log(complex_gamma(b*i)))-lgamma(gama)-lgamma(gama+x))
      return(respuesta)
    }

    pstart <- c(bstart,log(gammastart))

  }else {
    logL <- function(p){
      b <- p[1]
      gama <- p[2]
      respuesta <- -sum(2*Re(complex_gamma(gama+b*i,log=TRUE)+log(complex_gamma(b*i+x))-log(complex_gamma(b*i)))-lgamma(gama)-lgamma(gama+x))
      return(respuesta)
    }

    pstart <- c(bstart,gammastart)
  }

  #Optimization process
  converged<-FALSE
  if (method == "nlm") {
    methodText = "nlm"
    fit <- try(nlm(logL, p = pstart, hessian = TRUE, iterlim = control$maxit, print.level = control$trace))
    if ('try-error' %in% class(fit)){
      if (control$trace > 0) {
        cat("Crashed 'nlm' initial fit", "\n")
      }
    }
    else {
      fit$value <- fit$minimum
      fit$par <- fit$estimate
      fit$convergence <- fit$code
      if (fit$convergence < 3)
        converged<-TRUE
    }
  }
  else if (any(method == c("Nelder-Mead", "BFGS", "CG", "SANN"))) {
    methodText <- method
    fit <- try(optim(pstart, logL, method = method, hessian = TRUE,
                 control = list(maxit = control$maxit, trace = control$trace)))
    if ('try-error' %in% class(fit)){
      if (control$trace > 0) {
        cat(paste("Crashed '",method,"' initial fit",sep=","),"\n")
      }
    }
    else {
      if (fit$convergence == 0)
        converged = TRUE
    }
  }
  else if (any(method == c("L-BFGS-B"))) {
    methodText <- method
    fit<-try(optim(pstart, logL, method = method, lower = c(0,0.0000001), upper = c(Inf,Inf), hessian = TRUE, control = list(maxit = control$maxit, trace = control$trace)))
    if ('try-error' %in% class(fit)){
      if (control$trace > 0) {
        cat(paste("Crashed '",method,"' initial fit",sep=","),"\n")
      }
    }
    else {
      if (fit$convergence == 0)
        converged<-TRUE
    }
  }else{
    stop("Incorrect method")
  }
  if (converged){
    if (method != "L-BFGS-B") {
      coef.table <- rbind(fit$par, deparse.level = 0)
      dimnames(coef.table) <- list("  ", c("b", "log(gamma)"))
    }
    else {
      coef.table <- rbind(fit$par, deparse.level = 0)
      dimnames(coef.table) <- list("  ", c("b", "gamma"))
    }
    
    #Results
    results <- list(
      x = x,
      n = length(x),
      loglik = - (fit$value + sum(lfactorial(x))),
      aic = 2 * (fit$value + sum(lfactorial(x))) + 4,
      bic = 2 * (fit$value + sum(lfactorial(x))) + 2 * log(length(x)),
      coefficients =  coef.table,
      code = fit$convergence,
      hessian = fit$hessian,
      cov = solve(fit$hessian),
      se = sqrt(diag(solve(fit$hessian))),
      corr = solve(fit$hessian) / (sqrt(diag(solve(fit$hessian))) %o%
                                     sqrt(diag(solve(fit$hessian)))),
      code = fit$convergence,
      converged = converged,
      initialValues = c(bstart, gammastart),
      method = methodText
    )
  }else{
    results<-list()
    results$method<-methodText
    results$converged<-FALSE
    results$n<-nrow(x)
  }
  #restore warning level
  options(warn=defWarn)
  class(results) <- "fitCBP"
  results
}

#' Maximum-likelihood fitting of the EBW distribution
#'
#' @description
#' Maximum-likelihood fitting of the Extended Biparametric Waring (EBW) distribution with parameters \eqn{\alpha}, \eqn{\rho} and \eqn{\gamma}. Generic
#' methods are \code{print}, \code{summary}, \code{coef}, \code{logLik}, \code{AIC}, \code{BIC} and \code{plot}. The method to be used in fitting the 
#' model is "L-BFGS-B" which allows constraints for each variable (see details in \code{\link{optim}} funtion). 
#'
#' @usage
#' fitebw(x, alphastart = NULL, rhostart = NULL, gammastart = NULL, 
#'           method = "L-BFGS-B", control = list(),...)
#'        
#' @param x A numeric vector of length at least one containing only finite values.
#' @param alphastart A starting value for the parameter \eqn{\alpha}; by default \code{NULL}.
#' @param gammastart A starting value for the parameter \eqn{\gamma>max(0,2\alpha)}; by default \code{NULL}.
#' @param rhostart A starting value for the parameter \eqn{\rho>0}; by default \code{NULL}.
#' @param method The method to be used in fitting the model. The default method is "L-BFGS-B" (optim).
#' @param control A list of parameters for controlling the fitting process.
#' @param ...  Additional parameters.
#'
#' @details If the starting value for \eqn{\alpha} is positive, the parameterization \eqn{(\alpha,\rho)} is used; 
#' otherwise, the parameterization \eqn{(\alpha,\gamma)} is used.
#' 
#' If the starting values of the parameters \eqn{\alpha}, \eqn{\gamma} or \eqn{\rho} are omitted (default option), 
#' they are computing by the method of moments (if possible; otherwise they must be entered).
#' 
#' The default method is \code{"L-BFGS-B"} (see details in \code{\link{optim}} function),
#' but non-linear minimization (\code{\link{nlm}}) or those included in the \code{optim} function 
#' (\code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"} and \code{"SANN"}) may be used.
#' 
#' Standard error (SE) estimates for \eqn{\alpha}, \eqn{\gamma} or \eqn{\rho} are provided by the default method; 
#' otherwise, SE for \eqn{\alpha_0} and \eqn{\gamma_0} where \eqn{\alpha=-exp(\alpha_0)} and \eqn{\gamma=exp(\gamma_0)}
#' (or for \eqn{\alpha_0} and \eqn{\rho_0} where \eqn{\alpha=exp(\alpha_0)} and \eqn{\rho=exp(\rho_0)}) are computed.
#'
#' @return An object of class \code{'fitEBW'} is a list containing the following components:
#'
#' \itemize{
#' \item \code{n}, the number of observations,
#' \item \code{initialValues}, a vector with the starting values used,
#' \item \code{coefficients}, the parameter ML estimates of the CTP distribution,
#' \item \code{se}, a vector of the standard error estimates,
#' \item \code{hessian}, a symmetric matrix giving an estimate of the Hessian at the solution found in the optimization of the log-likelihood function,
#' \item \code{cov}, an estimate of the covariance matrix of the model coefficients,
#' \item \code{corr}, an estimate of the correlation matrix of the model estimates,
#' \item \code{loglik}, the maximized log-likelihood,
#' \item \code{aic}, Akaike Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters,
#' \item \code{bic}, Bayesian Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters,
#' \item \code{code}, a code that indicates successful convergence of the fitter function used (see nlm and optim helps),
#' \item \code{converged},  logical value that indicates if the optimization algorithms succesfull.
#' \item \code{method}, the name of the fitter function used.
#' }
#' 
#' Generic functions:
#'
#' \itemize{
#' \item \code{print}: The print of a \code{'fitEBW'} object shows the ML parameter estimates and their standard errors in parenthesis.
#' \item \code{summary}: The summary provides the ML parameter estimates, their standard errors and the statistic and p-value of the Wald test to check if the parameters are significant.
#' This summary also shows the loglikelihood, AIC and BIC values, as well as the results for the chi-squared goodness-of-fit test and the Kolmogorov-Smirnov test for discrete variables. Finally, the correlation matrix between parameter estimates appears.
#' \item \code{coef}: It extracts the fitted coefficients from a \code{'fitEBW'} object.
#' \item \code{logLik}: It extracts the estimated log-likelihood from a \code{'fitEBW'} object.
#' \item \code{AIC}: It extracts the value of the Akaike Information Criterion from a \code{'fitEBW'} object.
#' \item \code{BIC}: It extracts the value of the Bayesian Information Criterion from a \code{'fitEBW'} object.
#' \item \code{plot}: It shows the plot of a \code{'fitEBW'} object. Observed and theoretical probabilities, empirical and theoretical cumulative distribution functions or empirical cumulative probabilities against theoretical cumulative probabilities are the three plot types.
#' }
#' 
#' @importFrom stats nlm optim coef runif var
#' @importFrom utils data
#' @export
#' @references 
#' 
#' \insertRef{COR2021}{cpd}
#'  
#' 
#' @seealso
#' Plot of observed and theoretical frequencies for a EBW fit: \code{\link{plot.fitEBW}}
#' 
#' Maximum-likelihood fitting for the CTP distribution: \code{\link{fitctp}}.
#' 
#' Maximum-likelihood fitting for the CBP distribution: \code{\link{fitcbp}}.
#'
#' @examples
#' set.seed(123)
#' x <- rebw(500, 2, rho = 5)
#' fitebw(x)
#' summary(fitebw(x, alphastart = 1, rhostart = 5))

fitebw <- function(x, alphastart = NULL, rhostart = NULL, gammastart = NULL, 
                      method = "L-BFGS-B", control = list(), ...)
{
  
  #Checking
  if ( mode(x) != "numeric")
    stop( paste ("non-numeric argument to mathematical function", sep = ""))
  
  if( is.null(control$maxit) )
    control$maxit = 10000
  
  if( is.null(control$trace) )
    control$trace = 0
  if( is.null(control$warn) )
    control$warn <- -2
  #set warning level
  defWarn <- getOption("warn")
  options(warn=control$warn)
  
  if ( !is.numeric(control$maxit) || control$maxit <= 0 )
    stop( "maximum number of iterations must be > 0" )
  
  if (!any(method == c("nlm", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"))){
    stop("method must be in c('nlm', 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN')")
  }
  moments<-TRUE
  if ( is.numeric(alphastart) && !is.nan(alphastart) && (
        (alphastart>0 && is.numeric(rhostart) && !is.nan(rhostart) && rhostart > 0) ||
        (alphastart<0 && is.numeric(gammastart) && !is.nan(gammastart) && gammastart > 0)
      )){
    moments <-FALSE
  }

if (method != "L-BFGS-B") {
  #log-likelihood as a function of alpha and rho
  logL1<- function(p){
    alpha <-   exp( p[1] )
    rho   <-    exp( p[2] ) 
    respuesta <- - sum(2 * lgamma(rho + alpha) + 2 * lgamma(alpha + x) - 2 * lgamma(alpha) - lgamma(rho) - lgamma(rho + 2 * alpha + x))
    return(respuesta)
  }
  #log-likelihood as a function of alpha and gamma
  logL2 <- function(p){
    alpha <-  - exp( p[1] )
    gamma <-    exp( p[2] ) 
    respuesta <- - sum(2 * lgamma(gamma - alpha) + 2 * lgamma(alpha + x) - 2 * lgamma(alpha) - lgamma(gamma - 2 * alpha) - lgamma(gamma + x))
    return(respuesta)
  }
} 
else{
  #log-likelihood as a function of alpha and rho
  logL1<- function(p){
    alpha <- p[1]
    rho <- p[2] 
    respuesta <- - sum(2 * lgamma(rho + alpha) + 2 * lgamma(alpha + x) - 2 * lgamma(alpha) - lgamma(rho) - lgamma(rho + 2 * alpha + x))
    return(respuesta)
  }
  #log-likelihood as a function of alpha and gamma
  logL2 <- function(p){
    alpha <- p[1]
    gamma <- p[2]
    respuesta <- - sum(2 * lgamma(gamma - alpha) + 2 * lgamma(alpha + x) - 2 * lgamma(alpha) - lgamma(gamma - 2 * alpha) - lgamma(gamma + x))
    return(respuesta)
  }
} 
  
  if (moments){
  	#Estimates by the method of moments
  	m1 <- mean(x)
  	m2 <- var(x)  
    alphastart1 <- (m1 ^ 2 + sqrt(m1 * m2 * ((m1 - 1) * m1 + m2))) / (m2 - m1)
    alphastart2 <- 2 * m1 ^ 2 / (m2 - m1) - alphastart1
    alphastart  <- c(alphastart1, alphastart2)
    
    if (m1 > 1){
      indexes <- which((alphastart <= - m1 + sqrt(m1 * (m1 - 1))) & (alphastart >= - m1 - sqrt(m1 * (m1 - 1))))
      alphastart <- alphastart[-indexes]    																																																				   
    }
    
    if (length(alphastart) < 1)
      stop("The method of moments does not provide any estimate. Enter initial values for the parameters.")
    
    param2 <- c()
    logLVect <- c()
    for (i in 1:length(alphastart)){
      if (alphastart[i] >0){
        param2 <- c(param2, alphastart[i] ^ 2 / m1 + 1)
        logLVect<-c(logLVect, logL1)
      } else {
        param2 <- c(param2, alphastart[i] ^ 2 / m1 + 2 * alphastart[i] + 1)
        logLVect<-c(logLVect, logL2)
      }
    }
  }
  else{
    if (alphastart < 0) {
      if (alphastart %% 1 == 0)
        stop("'alpha' cannot be a negative integer")
      if (missing(gammastart))
        stop("Introduce 'gammastart'")
      if (gammastart < 0)
        stop("'gammastart' must be positive")
      param2 <- c(gammastart)
      logLVect <- c(logL2)
    }else{ 
      if (missing(rhostart))
        stop("Introduce 'rhostart'")
      if (rhostart <= 0)
        stop("'rhostart' must be positive")
      param2 <- c(rhostart)
      logLVect <- c(logL1)
    }
    alphastart <- c(alphastart)
  	if ((mean(x) > var(x)) && alphastart > 0) warning("With underdispersed data alpha must be negative")
  }


  #Log-likelihood
  results <- NULL
  for (ind in 1:length(alphastart)){
    if (method != "L-BFGS-B") {
      if ( alphastart[[ind]] < 0 )
        pstart <- c(log(-alphastart[[ind]]),log(param2[[ind]]))
      else
        pstart<-c(log(alphastart[[ind]]),log(param2[[ind]]))
    }else{
      pstart<-c(alphastart[[ind]],param2[[ind]])
    }
	  converged <- FALSE
	  if (method == "nlm"){
	    methodText <- "nlm"
  		fit <- try (nlm(logLVect[[ind]],p=pstart, hessian = TRUE, 
  						iterlim = control$maxit, print.level = control$trace),silent = TRUE)
  		if ('try-error' %in% class(fit)){
  		  if (control$trace > 0) {
  			  cat(paste("Crashed '","nlm","' initial fit",sep=","),"\n")
  		  }	
  		  converged <- FALSE
  		}
  		else {
  		  fit$value <- fit$minimum
  		  fit$par <- fit$estimate
  		  fit$convergence <- fit$code
  		  if (fit$convergence < 3)
  			  converged <- TRUE
  		  else
  		    converged <- FALSE
  		}
	  }else if (any(method == c("Nelder-Mead", "BFGS", "CG", "SANN"))) {
	    methodText <- method
  		fit <- try(optim(pstart, logLVect[[ind]], method = method, 
  					hessian = TRUE, control = list(maxit = control$maxit, trace = control$trace)), silent=TRUE)
  		if ('try-error' %in% class(fit)){
  		  if (control$trace > 0) {
  			  cat(paste("Crashed '",method,"' initial fit",sep=","),"\n")
  		  }
  		  converged <- FALSE
  		}
  		else {
  		  if (fit$convergence == 0)
  			  converged <- TRUE
  		  else
  		    converged <- FALSE
  		}
	  }
	  else if (any(method == c("L-BFGS-B"))) {
	    methodText <- method
  		fit <- try(optim(pstart, logLVect[[ind]], 
  						lower = c(-Inf, 1e-10), upper = c(Inf, Inf), method = method, hessian = TRUE,
  					    control = list(maxit = control$maxit, trace = control$trace)), silent=TRUE)
    		if ('try-error' %in% class(fit)){
    		  if (control$trace > 0) {
    			  cat(paste("Crashed '",method,"' initial fit",sep=","),"\n")
    		  }
    		  converged <- FALSE
    		}
    		else {
    		  if (fit$convergence == 0)
    			  converged <- TRUE
    		  else
    		    converged <- FALSE
    		}
	  }else{
		  stop("Incorrect method")
	  }
    #Hessian exists
	  existHessian<-try(solve(fit$hessian))
	  if ('try-error' %in% class(existHessian)){
	     existHessian <- FALSE
	  }else{
	     existHessian <- TRUE
	  }
	  
    if (is.null(results) || results$converged == FALSE || (converged && (results$aic > (2 * (fit$value + sum(lfactorial(x))) + 4) )) ){
      if (converged && existHessian){
    		if (any(method == c("nlm", "Nelder-Mead", "BFGS", "CG","SANN"))){
    		  if (alphastart[ind]>0){
      			coef.table<-rbind(exp(fit$par),deparse.level=0)
      			dimnames(coef.table)<-list("",c("alpha","rho"))
        		se<-rbind(sqrt(diag(solve(fit$hessian))),deparse.level=0)
        		dimnames(se)<-list("",c("std error log(alpha)","std error log(rho)"))
    		  }
    		  else{
      			coef.table<-rbind(c(-exp(fit$par[1]),exp(fit$par[2])),deparse.level=0)
      			dimnames(coef.table)<-list("",c("alpha","gamma"))
        		se <-rbind(sqrt(diag(solve(fit$hessian))),deparse.level=0)
      			dimnames(se)<-list("",c("std error log(-alpha)","std error log(gamma)"))
      		}
    		}
    		else{
    		  if (alphastart[ind]>0){
      			coef.table<-rbind(c(fit$par,fit$par[2]+2*fit$par[1]),deparse.level=0)
      			dimnames(coef.table)<-list("",c("alpha","rho","gamma"))
        		se<-solve(fit$hessian)
        		se<-rbind(sqrt(c(diag(se),se[2,2]+4*se[1,1]+4*se[1,2])),deparse.level=0)
        		dimnames(se)<-list("",c("std error alpha","std error rho","std error gamma"))
    		  } else {
      			coef.table<-rbind(c(fit$par[1],fit$par[2]),deparse.level=0)
      			dimnames(coef.table)<-list("",c("alpha","gamma"))
        		se = rbind(sqrt(diag(solve(fit$hessian))))
        		dimnames(se)<-list("",c("std error alpha","std error gamma"))
    		  }
    		}
    		results<-list(
    		  x = x,
    		  n=length(x),
    		  call = call,
    		  loglik = - (fit$value + sum(lfactorial(x))),
    		  aic = 2 * (fit$value + sum(lfactorial(x))) + 4,
    		  bic = 2 * (fit$value + sum(lfactorial(x))) + 2 * log(length(x)),
    		  coefficients = coef.table,
    		  #expected.frequencies=length(data)*dctp(0:max(data),fit$par[1],fit$par[2],fit$par[3]),
    		  hessian = fit$hessian,
    		  cov = solve(fit$hessian),
    		  se = se, 
    		  corr = solve(fit$hessian)/(sqrt(diag(solve(fit$hessian))) %o% 
    		                                      sqrt(diag(solve(fit$hessian)))),
    		  code = fit$convergence, 
    		  converged = converged,
    		  method = methodText
    		)
  	  }
  	  else{
    		warning("fitebw have not converged")
    		results<-list(
    		  x = x,
    		  n=length(x),
    		  call = call,
    		  code = fit$convergence, 
    		  converged = fit$converged,
    		  method = methodText
    		)
  	  }
    }
  }
  
  #Restore warning level
  options(warn=defWarn)
  class(results)<-"fitEBW"
  return(results)
}

#' @export
logLik.fitCTP <- function (object, ...){
  val <- object$loglik
  attr(val, "nobs") <- object$n
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  val
}

#' @method print fitCTP
#' @export
print.fitCTP<-function (x, digits = getOption("digits"), ...) {
  if (length(coef(x))) {
    ans=format(rbind(dimnames(x$coefficients)[[2L]],
                     format(x$coefficients,digits = digits),
              sapply(format(x$se,digits = digits),function(v){paste("(",v,")",sep="")})),justify = "centre")
    colnames(ans)<-ans[1,]
    ans<-ans[2:3,]
    print.default(ans, print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")

  if(!x$converged){
    cat("Error of convergence")
  }
  invisible(x)
}

#' @export
logLik.fitCBP <- function (object, ...){
  val <- object$loglik
  attr(val, "nobs") <- object$n
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  val
}

#' @method print fitCBP
#' @export
print.fitCBP<-function (x, digits = getOption("digits"), ...) {
  if (length(coef(x))) {
    ans=format(rbind(dimnames(x$coefficients)[[2L]],
                     format(x$coefficients,digits = digits),
                     sapply(format(x$se,digits = digits),function(v){paste("(",v,")",sep="")})),justify = "centre")
    colnames(ans)<-ans[1,]
    ans<-ans[2:3,]
    print.default(ans, print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
  
  if(!x$converged){
    cat("Error of convergence")
  }
  invisible(x)
}

#' @method logLik fitEBW
#' @export
logLik.fitEBW <- function (object, ...){
  val <- object$loglik
  attr(val, "nobs") <- object$n
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  val
}

#' @method print fitEBW
#' @export
print.fitEBW<-function (x, digits = getOption("digits"), ...) {
  if (length(coef(x))) {
    ans=format(rbind(dimnames(x$coefficients)[[2L]],
                     format(x$coefficients,digits = digits),
              sapply(format(x$se,digits = digits),function(v){paste("(",v,")",sep="")})),justify = "centre")
    colnames(ans)<-ans[1,]
    ans<-ans[2:3,]
    print.default(ans, print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")

  if(!x$converged){
    cat("Error of convergence")
  }
  invisible(x)
}



#' @method summary fitCTP
#' @importFrom stats pchisq pnorm stepfun
#' @importFrom dgof ks.test
#' @export
summary.fitCTP<-function (object, ...) {
  if (!("fitCTP" %in% class(object))){ 
    stop("This method only works with fitCTP objects")
  }

  myfun <- stepfun(0:max(object$x), c(0, pctp(0:max(object$x), a = object$coefficients[1], b = object$coefficients[2], gamma = object$coefficients[3])))
  object$kstest <- ks.test(object$x,myfun , simulate.p.value=TRUE, B=1000)
  object$zvalue<-object$coefficients/object$se
  object$pvalue<- 2 * pnorm(abs(object$zvalue),lower.tail = FALSE)
  xmax<-max(object$x)
  p.esp<-dctp(0:(xmax-1),object$coefficients[1],object$coefficients[2],object$coefficients[3])
  p.esp[xmax+1]<-1-sum(p.esp[1:(xmax)])
  obs<-rep(0,xmax+1)
  for (i in 1:length(object$x))
    obs[object$x[i]+1]<-obs[object$x[i]+1]+1
  object$chi2test=chisq.test2(obs, p.esp, npar = 3, ...)
  class(object) <- "summary.fitCTP"
  object
}


#' @method print summary.fitCTP
#' @export
print.summary.fitCTP <- function(x, digits = getOption("digits"), ...){
  if (!("summary.fitCTP" %in% class(x))){ 
    stop("This method only works with fitCTP objects")
  }
  coefName<-rbind("",cbind(dimnames(x$coefficients)[[2L]]))
  coefValu<-rbind("Estimate",cbind(as.vector(format(x$coefficients,digits = digits))))
  se<-rbind("Std. Error",cbind(as.vector(format(x$se,digits = digits))))
  zvalue<-rbind("z-value",cbind(as.vector(format(x$zvalue,digits = digits))))
  prz<-rbind("Pr(>|z|)",cbind(as.vector(format(x$pvalue,digits = digits))))
  ans<-cbind(format(cbind(coefValu,se,zvalue,prz),justify = "centre"))
  cat("Parameters:")
  prmatrix(ans, rowlab=coefName, collab=rep("",4),quote=FALSE)
  cat(paste("\nLoglikelihood: ",round(x$loglik,2),"   AIC: ",round(x$aic,2),"   BIC: ",round(x$bic,2),"\n",sep=""))
  cat("\nGoodness-of-fit tests:\n")
  cat(paste("Chi-2: ", format(x$chi2test$statistic,digits=digits)," (p-value: ",x$chi2test$p.value,")   Kolmogorov-Smirnov: ",format(x$kstest$statistic,digits=digits)," (p-value: ",x$kstest$p.value,")\n",sep=""))
  cat("\nCorrelation Matrix:\n")
  prmatrix(x$corr, rowlab=dimnames(x$coefficients)[[2L]], collab=dimnames(x$coefficients)[[2L]],quote=FALSE)
  invisible(x)
}

#' @method summary fitCBP
#' @importFrom stats pchisq pnorm stepfun
#' @importFrom dgof ks.test
#' @export
summary.fitCBP<-function (object, ...) {
  if (!("fitCBP" %in% class(object))){ 
    stop("This method only works with fitCBP objects")
  }
  myfun <- stepfun(0:max(object$x), c(0, pcbp(0:max(object$x), b = object$coefficients[1], gamma = object$coefficients[2])))
  object$kstest <- ks.test(object$x,myfun , simulate.p.value=TRUE, B=1000)
  object$zvalue<-object$coefficients/object$se
  object$pvalue<- 2 * pnorm(abs(object$zvalue),lower.tail = FALSE)
  xmax<-max(object$x)
  p.esp<-dcbp(0:(xmax-1),object$coefficients[1],object$coefficients[2])
  p.esp[xmax+1]<-1-sum(p.esp[1:(xmax)])
  obs<-rep(0,xmax+1)
  for (i in 1:length(object$x))
    obs[object$x[i]+1]<-obs[object$x[i]+1]+1
  object$chi2test=chisq.test2(obs, p.esp, npar = 2, ...)
  class(object) <- "summary.fitCBP"
  object
}


#' @method print summary.fitCBP
#' @export
print.summary.fitCBP <- function(x, digits = getOption("digits"), ...){
  if (!("summary.fitCBP" %in% class(x))){ 
    stop("This method only works with fitCBP objects")
  }
  coefName<-rbind("",cbind(dimnames(x$coefficients)[[2L]]))
  coefValu<-rbind("Estimate",cbind(as.vector(format(x$coefficients,digits = digits))))
  se<-rbind("Std. Error",cbind(as.vector(format(x$se,digits = digits))))
  zvalue<-rbind("z-value",cbind(as.vector(format(x$zvalue,digits = digits))))
  prz<-rbind("Pr(>|z|)",cbind(as.vector(format(x$pvalue,digits = digits))))
  ans<-cbind(format(cbind(coefValu,se,zvalue,prz),justify = "centre"))
  cat("Parameters:")
  prmatrix(ans, rowlab=coefName, collab=rep("",4),quote=FALSE)
  cat(paste("\nLoglikelihood: ",round(x$loglik,2),"   AIC: ",round(x$aic,2),"   BIC: ",round(x$bic,2),"\n",sep=""))
  cat("\nGoodness-of-fit tests:\n")
  cat(paste("Chi-2: ", format(x$chi2test$statistic,digits=digits)," (p-value: ",x$chi2test$p.value,")   Kolmogorov-Smirnov: ",format(x$kstest$statistic,digits=digits)," (p-value: ",x$kstest$p.value,")\n",sep=""))
  cat("\nCorrelation Matrix:\n")
  prmatrix(x$corr, rowlab=dimnames(x$coefficients)[[2L]], collab=dimnames(x$coefficients)[[2L]],quote=FALSE)
  invisible(x)
}


#' @method summary fitEBW
#' @importFrom stats pchisq pnorm stepfun
#' @importFrom dgof ks.test
#' @export
summary.fitEBW<-function (object, ...) {
  if (!("fitEBW" %in% class(object))){ 
    stop("This method only works with fitEBW objects")
  }
  if (object$coefficients[1]<0)
     myfun <- stepfun(0:max(object$x), c(0, pebw(0:max(object$x), alpha = object$coefficients[1], gamma = object$coefficients[2])))
  else
     myfun <- stepfun(0:max(object$x), c(0, pebw(0:max(object$x), alpha = object$coefficients[1], rho = object$coefficients[2])))
	  
  object$kstest <- ks.test(object$x,myfun , simulate.p.value=TRUE, B=1000)
  object$zvalue<-object$coefficients/object$se
  object$pvalue<- 2 * pnorm(abs(object$zvalue),lower.tail = FALSE)
  xmax<-max(object$x)						   
  if (object$coefficients[1]<0)
    p.esp<-debw(0:(xmax-1), alpha = object$coefficients[1], gamma = object$coefficients[2])
  else
    p.esp<-debw(0:(xmax-1), alpha = object$coefficients[1], rho = object$coefficients[2])
  p.esp[xmax+1]<-1-sum(p.esp[1:(xmax)])
  obs<-rep(0,xmax+1)
  for (i in 1:length(object$x))
    obs[object$x[i]+1]<-obs[object$x[i]+1]+1
  object$chi2test=chisq.test2(obs, p.esp, npar = 2, ...)
  class(object) <- "summary.fitEBW"
  object
}


#' @method print summary.fitEBW
#' @export
print.summary.fitEBW <- function(x, digits = getOption("digits"), ...){
  if (!("summary.fitEBW" %in% class(x))){ 
    stop("This method only works with fitEBW objects")
  }
  coefName<-rbind("",cbind(dimnames(x$coefficients)[[2L]]))
  coefValu<-rbind("Estimate",cbind(as.vector(format(x$coefficients,digits = digits))))
  se<-rbind("Std. Error",cbind(as.vector(format(x$se,digits = digits))))
  zvalue<-rbind("z-value",cbind(as.vector(format(x$zvalue,digits = digits))))
  prz<-rbind("Pr(>|z|)",cbind(as.vector(format(x$pvalue,digits = digits))))
  ans<-cbind(format(cbind(coefValu,se,zvalue,prz),justify = "centre"))
  cat("Parameters:")
  prmatrix(ans, rowlab=coefName, collab=rep("",4),quote=FALSE)
  cat(paste("\nLoglikelihood: ",round(x$loglik,2),"   AIC: ",round(x$aic,2),"   BIC: ",round(x$bic,2),"\n",sep=""))
  cat("\nGoodness-of-fit tests:\n")
  cat(paste("Chi-2: ", format(x$chi2test$statistic,digits=digits)," (p-value: ",x$chi2test$p.value,")   Kolmogorov-Smirnov: ",format(x$kstest$statistic,digits=digits)," (p-value: ",x$kstest$p.value,")\n",sep=""))
  cat("\nCorrelation Matrix:\n")
  prmatrix(x$corr, rowlab=dimnames(x$coefficients)[[2L]], collab=dimnames(x$coefficients)[[2L]],quote=FALSE)
  invisible(x)
}


#' Plot of observed and theoretical frequencies for a CBP fit
#' 
#' @param x An object of class \code{'fitCBP'}
#' @param plty Plot type to be shown. Default is \code{"FREQ"} which shows the observed and theoretical frequencies for each value of the variable; \code{"CDF"} and \code{"PP"} are also available for plotting the empirical and theoretical cumulative distribution functions or the theoretical cumulative probabilities against the empirical cumulative probabilities, respectively.
#' @param maxValue maxValue you want to appear in the plot
#' @param ...  Additional parameters.
#' @importFrom graphics abline legend plot points segments
#' @method plot fitCBP
#' 
#' @examples 
#' set.seed(123)
#' x <- rcbp(500, 1.75, 3.5)
#' fit <- fitcbp(x)
#' plot(fit)
#' plot(fit, plty = "CDF")
#' plot(fit, plty = "PP")
#'  
#' @export
plot.fitCBP <- function(x,plty="FREQ",maxValue=NULL,...){
  if (!("fitCBP" %in% class(x))){ 
    stop("This method only works with fitCBP objects")
  }
  hLimit<-max(x$x)
  if (is.numeric(maxValue) && maxValue%%1==0){
    if (maxValue>hLimit || maxValue < 0){
      warning("maxValue can't be greater than maximum value neither negative. maxValue will be ignored")
      maxValue=hLimit
    }
  }else{
    maxValue=hLimit
  }
  values<-0:hLimit;
  
  p<-pcbp(values,x$coefficients[1],x$coefficients[2] )
  freq<-rep(0,hLimit+1)
  for (i in 1:x$n)
    freq[x$x[i]+1]<-freq[x$x[i]+1]+1
  cumFreq=cumsum(freq)
  cumFreq=cumFreq/x$n
  if (plty=="PP"){
    plot(range(0, 1), range(0, 1), type = "n", xlab = "Empirical Cumulative Probabilities", 
         ylab = "Theoretical Cumulative Pobabilities",main="PP plot",...)
    points(cumFreq,p,col="black",pch=19)
    abline(0,1, col = "blue", lty = 2)
  } else if(plty=="CDF") {
    plot(range(0, maxValue), range(0, 1), type = "n", xlab = "Values", 
         ylab = "Cumulative probability",main="Empirical and Theoretical CDFs")
    cFreq<-cumFreq[values+1]
    points(values,p,col="blue",pch=19)
    points(values,cFreq,col="red",pch=19)
    segments(values[-length(values)], cumFreq[-length(cumFreq)], values[-1], cumFreq[-length(cumFreq)], col= 'red')
    segments(values[-length(values)], p[-length(p)], values[-1], p[-length(p)], col= 'blue')
    abline(h = c(0, 1), col = "gray", lty = 2)
    legend(hLimit*4/5,0.25, legend=c("Empirical", "Theoretical"),
           col=c("red", "blue"), lty=1, cex=0.8)
  } else{
    n.esp<-dcbp(values,x$coefficients[1],x$coefficients[2] )*x$n
    fLimit <- max(n.esp,freq)
    plot(range(0, maxValue), range(0, fLimit), type = "n", xlab = "Values", 
         ylab = "Frequencies",main="Observed & Theoretical Frequencies")
    points(values,n.esp,col="blue",pch=19)
    points(values,freq,col="red",pch=19)
    segments(values[-length(values)], freq[-length(freq)], values[-1], freq[-1], col= 'red',lty=2)
    segments(values[-length(values)], n.esp[-length(n.esp)], values[-1], n.esp[-1], col= 'blue',lty=2)
    abline(h = 0, col = "gray", lty = 2)
    legend(hLimit*4/5,fLimit*4/5, legend=c("Observed", "Theoretical"),
           col=c("red", "blue"), lty=2, cex=0.8)
  }
}


#' Plot of observed and theoretical frequencies for a CTP fit
#' 
#' @param x An object of class \code{'fitCTP'}
#' @param plty Plot type to be shown. Default is \code{"FREQ"} which shows the observed and theoretical frequencies for each value of the variable; \code{"CDF"} and \code{"PP"} are also available for plotting the empirical and theoretical cumulative distribution functions or the theoretical cumulative probabilities against the empirical cumulative probabilities, respectively.
#' @param maxValue maxValue you want to appear in the plot
#' @param ...  Additional parameters.
#' @importFrom graphics abline legend plot points segments
#' @method plot fitCTP
#' 
#' @examples 
#' set.seed(123)
#' x <- rctp(500, -0.5, 1, 2)
#' fit <- fitctp(x)
#' plot(fit)
#' plot(fit, plty = "CDF")
#' plot(fit, plty = "PP")
#' 
#' @export
plot.fitCTP <- function(x,plty="FREQ",maxValue=NULL,...){
  if (!("fitCTP" %in% class(x))){ 
    stop("This method only works with fitCTP objects")
  }
  hLimit<-max(x$x)
  if (is.numeric(maxValue) && maxValue%%1==0){
    if (maxValue>hLimit || maxValue < 0){
      warning("maxValue can't be greater than maximum value neither negative. maxValue will be ignored")
      maxValue=hLimit
    }
  }else{
    maxValue=hLimit
  }
  values<-0:hLimit;
  
  p<-pctp(values,x$coefficients[1],x$coefficients[2],x$coefficients[3] )
  freq<-rep(0,hLimit+1)
  for (i in 1:x$n)
    freq[x$x[i]+1]<-freq[x$x[i]+1]+1
  cumFreq=cumsum(freq)
  cumFreq=cumFreq/x$n
  if (plty=="PP"){
    plot(range(0, 1), range(0, 1), type = "n", xlab = "Empirical Cumulative Probabilities", 
         ylab = "Theoretical Cumulative Probalilities",main="PP plot",...)
    points(cumFreq,p,col="black",pch=19)
    abline(0,1, col = "blue", lty = 2)
  } else if(plty=="CDF") {
    plot(range(0, maxValue), range(0, 1), type = "n", xlab = "Values", 
         ylab = "Cumulative probability",main="Empirical and Theoretical CDF")
    cFreq<-cumFreq[values+1]
    points(values,p,col="blue",pch=19)
    points(values,cFreq,col="red",pch=19)
    segments(values[-length(values)], cumFreq[-length(cumFreq)], values[-1], cumFreq[-length(cumFreq)], col= 'red')
    segments(values[-length(values)], p[-length(p)], values[-1], p[-length(p)], col= 'blue')
    abline(h = c(0, 1), col = "gray", lty = 2)
    legend(maxValue*4/5,0.25, legend=c("Empirical", "Theoretical"),
           col=c("red", "blue"), lty=1, cex=0.8)
  } else{
    n.esp<-dctp(values,x$coefficients[1],x$coefficients[2],x$coefficients[3] )*x$n
    fLimit <- max(n.esp,freq)
    plot(range(0, maxValue), range(0, fLimit), type = "n", xlab = "Values", 
         ylab = "Frequencies",main="Observed & Theoretical Frequencies")
    points(values,n.esp,col="blue",pch=19)
    points(values,freq,col="red",pch=19)
    segments(values[-length(values)], freq[-length(freq)], values[-1], freq[-1], col= 'red',lty=2)
    segments(values[-length(values)], n.esp[-length(n.esp)], values[-1], n.esp[-1], col= 'blue',lty=2)
    abline(h = 0, col = "gray", lty = 2)
    legend(maxValue*4/5,fLimit*4/5, legend=c("Observed", "Theoretical"),
           col=c("red", "blue"), lty=2, cex=0.8)
  }
}

#' Plot of observed and theoretical frequencies for a EBW fit
#' 
#' @param x An object of class \code{'fitEBW'}
#' @param plty Plot type to be shown. Default is \code{"FREQ"} which shows the observed and theoretical frequencies for each value of the variable; \code{"CDF"} and \code{"PP"} are also available for plotting the empirical and theoretical cumulative distribution functions or the theoretical cumulative probabilities against the empirical cumulative probabilities, respectively.
#' @param maxValue maxValue you want to appear in the plot
#' @param ...  Additional parameters.
#' @importFrom graphics abline legend plot points segments
#' @method plot fitEBW
#' 
#' @examples 
#' set.seed(123)
#' x <- rebw(500, -0.25, 1)
#' fit <- fitebw(x)
#' plot(fit)
#' plot(fit, plty = "CDF")
#' plot(fit, plty = "PP")
#' 
#' @export
plot.fitEBW <- function(x,plty="FREQ",maxValue=NULL,...){
  if (!("fitEBW" %in% class(x))){ 
    stop("This method only works with fitEBW objects")
  }
  hLimit<-max(x$x)
  if (is.numeric(maxValue) && maxValue%%1==0){
    if (maxValue>hLimit || maxValue < 0){
      warning("maxValue can't be greater than maximum value neither negative. maxValue will be ignored")
      maxValue=hLimit
    }
  }else{
    maxValue=hLimit
  }
  values<-0:hLimit;
  if (x$coefficients[1]>0) #rho parameterization
    p<-pebw(values,alpha=x$coefficients[1],rho=x$coefficients[2])
  else
    p<-pebw(values,alpha=x$coefficients[1],gamma=x$coefficients[2])
  freq<-rep(0,hLimit+1)
  for (i in 1:x$n)
    freq[x$x[i]+1]<-freq[x$x[i]+1]+1
  cumFreq=cumsum(freq)
  cumFreq=cumFreq/x$n
  if (plty=="PP"){
    plot(range(0, 1), range(0, 1), type = "n", xlab = "Empirical Cumulative Probabilities", 
         ylab = "Theoretical Cumulative Probalilities",main="PP plot",...)
    points(cumFreq,p,col="black",pch=19)
    abline(0,1, col = "blue", lty = 2)
  } else if(plty=="CDF") {
    plot(range(0, maxValue), range(0, 1), type = "n", xlab = "Values", 
         ylab = "Cumulative probability",main="Empirical and Theoretical CDF")
    cFreq<-cumFreq[values+1]
    points(values,p,col="blue",pch=19)
    points(values,cFreq,col="red",pch=19)
    segments(values[-length(values)], cumFreq[-length(cumFreq)], values[-1], cumFreq[-length(cumFreq)], col= 'red')
    segments(values[-length(values)], p[-length(p)], values[-1], p[-length(p)], col= 'blue')
    abline(h = c(0, 1), col = "gray", lty = 2)
    legend(maxValue*4/5,0.25, legend=c("Empirical", "Theoretical"),
           col=c("red", "blue"), lty=1, cex=0.8)
  } else{
    if (x$coefficients[1]>0) #rho parameterization
      n.esp<-debw(values,alpha=x$coefficients[1],rho=x$coefficients[2])*x$n
    else
      n.esp<-debw(values,alpha=x$coefficients[1],gamma=x$coefficients[2])*x$n
    fLimit <- max(n.esp,freq)
    plot(range(0, maxValue), range(0, fLimit), type = "n", xlab = "Values", 
         ylab = "Frequencies",main="Observed & Theoretical Frequencies")
    points(values,n.esp,col="blue",pch=19)
    points(values,freq,col="red",pch=19)
    segments(values[-length(values)], freq[-length(freq)], values[-1], freq[-1], col= 'red',lty=2)
    segments(values[-length(values)], n.esp[-length(n.esp)], values[-1], n.esp[-1], col= 'blue',lty=2)
    abline(h = 0, col = "gray", lty = 2)
    legend(maxValue*4/5,fLimit*4/5, legend=c("Observed", "Theoretical"),
           col=c("red", "blue"), lty=2, cex=0.8)
  }
}

#' Variance decomposition for a EBW fit
#' 
#' One of the main drawbacks of the Univariate Generalized Waring (UGW) distribution with parameters \eqn{a}, 
#' \eqn{k} and \eqn{\rho} is that the first two parameters are interchangeable, so it is not possible to distinguish 
#' the variance components 'liability' and 'proneness' without additional information. To solve this problem, 
#' an EBW distribution (where these components are uniquely identifiable) can be used since, 
#' given a UGW distribution, there always exists an EBW close to it.
#' 
#' @param object An object of class \code{'fitEBW'}
#' @param ...  Additional parameters.
#' @return A data frame with the variance components (randomness, liability and proneness) in absolute and relative terms.
#'
#' @examples
#'
#' set.seed(123)
#' x <- rebw(500, 2, rho = 5)
#' fit <- fitebw(x, alphastart = 1, rhostart = 5)
#' varcomp(fit)
#' 
#' @export
varcomp <- function(object,  ...){
  UseMethod("varcomp")
}

#' @export
varcomp.fitEBW <- function(object ,...){
  if (!("fitEBW" %in% class(object))){ 
    stop("This method only works with fitEBW objects")
  }
  alpha<-object$coefficients[1]
  rho <- object$coefficients[2]
  if (alpha >0 && rho >2){
    #proneness
    prone <- (alpha ^ 3 * (alpha + rho - 1)) / ((rho - 1) ^ 2 * (rho - 2))
    #randomness
    rand <- alpha ^ 2 / (rho - 1)
    #liability
    liabi <- (alpha ^ 2 * (alpha + 1)) / ((rho - 1) * (rho - 2))
    var <- rand + liabi + prone
    Value <- c(rand, liabi, prone)
    Ratio <- c(rand / var, liabi / var, prone / var)
    ans <- data.frame(cbind(Value, Ratio), row.names = c("Randomness","Liability","Proneness"))
  }else{
    warning("alpha must be greater than 0 and rho greater than 2")
    ans<-NULL
  }
  ans
}






#' The Complex Triparametric Pearson (CTP) Distribution
#'
#' @description
#' Probability mass function, distribution function, quantile function and random generation for the Complex Triparametric Pearson (CTP) distribution with parameters \eqn{a}, \eqn{b} and \eqn{\gamma}. 
#'
#' @usage
#' dctp(x, a, b, gamma)
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param a parameter a (real)
#' @param b parameter b (real)
#' @param gamma parameter \eqn{\gamma} (positive)
#' @param lower.tail if TRUE (default), probabilities are \eqn{P(X<=x)}, otherwise, \eqn{P(X>x)}.
#'
#' @details
#' The CTP distribution with parameters \eqn{a}, \eqn{b} and \eqn{\gamma} has pmf
#' \deqn{f(x|a,b,\gamma) = C \Gamma(a+ib+x) \Gamma(a-ib+x) / (\Gamma(\gamma+x) x!), x=0,1,2,...} 
#' where \eqn{i} is the imaginary unit, \eqn{\Gamma(·)} the gamma function and 
#' \deqn{C = \Gamma(\gamma-a-ib) \Gamma(\gamma-a+ib) / (\Gamma(\gamma-2a) \Gamma(a+ib) \Gamma(a-ib))}
#' the normalizing constant.
#' 
#' If \eqn{a=0} the CTP is a Complex Biparametric Pearson (CBP) distribution, so the pmf of the CBP distribution is obtained.
#' If \eqn{b=0} the CTP is an Extended Biparametric Waring (EBW) distribution, so the pmf of the EBW distribution is obtained.																															  
#'
#' The mean and the variance of the CTP distribution are
#' \eqn{E(X)=\mu=(a^2+b^2)/(\gamma-2a-1)} and \eqn{Var(X)=\mu(\mu+\gamma-1)/(\gamma-2a-2)}
#' so \eqn{\gamma > 2a + 2}.
#'
#' It is underdispersed if \eqn{a < - (\mu + 1) / 2}, equidispersed if \eqn{a = - (\mu + 1) / 2} or overdispersed
#' if \eqn{a > - (\mu + 1) / 2}. In particular, if \eqn{a >= 0} the CTP is always overdispersed.
#'
#' @return 
#' \code{dctp} gives the pmf, \code{pctp} gives the distribution function, \code{qctp} gives the quantile function and \code{rctp} generates random values.
#'
#' If \eqn{a = 0} the probability mass function, distribution function, quantile function and random generation function for the CBP distribution arise.
#' If \eqn{b = 0} the probability mass function, distribution function, quantile function and random generation function for the EBW distribution arise.
#' 
#' @references 
#' \insertRef{RCS2003}{cpd}
#' 
#' \insertRef{RCSO2004}{cpd}
#' 
#' \insertRef{ROC2018}{cpd}
#' 
#' \insertRef{COR2021}{cpd}
#'
#' @seealso
#' Functions for maximum-likelihood fitting of the CTP, CBP and EBW distributions: \code{\link{fitctp}}, \code{\link{fitcbp}} and \code{\link{fitebw}}.
#'
#' @examples
#' # Examples for the function dctp
#' dctp(3,1,2,5)
#' dctp(c(3,4),1,2,5)
#'
#' @name ctp
#'
#' @rdname ctp
#' @importFrom hypergeo complex_gamma
#' @export
#'
dctp <- function(x, a, b, gamma) {

  if ( mode(c(x,a,b,gamma)) != "numeric")
    stop( "non-numeric argument to mathematical function" )

  if( (gamma <= 2 * a) || (gamma <= 0) )
    stop("gamma must be greater than max(0,2a)")

  
  x <- as.vector(x)
  errors<-c()


  for ( i in 1:length(x) ){
      if ( round(x[[i]]) != x[[i]] || x[[i]]<0){
        warning( paste ("non-integer or negative x[[", i, "]]", sep = ""))
        errors <- c(errors,i)
      }

  }

  if (b==0){
     lpmf<-2*(lgamma(gamma-a)+lgamma(x+a)-lgamma(a))-lgamma(gamma-2*a)-lgamma(gamma+x)-lgamma(x+1)
  } else {
     i <- sqrt(as.complex(-1))
     lpmf<-Re(complex_gamma(gamma-a-b*i,log=TRUE)+complex_gamma(gamma-a+b*i,log=TRUE)+complex_gamma(x+a+b*i,log=TRUE)+complex_gamma(x+a-b*i,log=TRUE)) -lgamma(gamma-2*a)-Re(complex_gamma(a+b*i,log=TRUE)+complex_gamma(a-b*i,log=TRUE))-lgamma(gamma+x)-lgamma(x+1)
  }
  result <- exp(lpmf)
  for (i in 1:length(errors)){
    result[ errors[[i]] ] = 0
  }
  return(result)
}


#' The Complex Biparametric Pearson (CBP) Distribution
#'
#' @description
#' Probability mass function, distribution function, quantile function and random generation for the Complex Biparametric Pearson (CBP) distribution with parameters \eqn{b} and \eqn{\gamma}. 
#'
#' @usage
#' dcbp(x, b, gamma)
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param b parameter b (real)
#' @param gamma parameter gamma (positive)
#' @param lower.tail if TRUE (default), probabilities are \eqn{P(X<=x)}, otherwise, \eqn{P(X>x)}.
#'
#' @details
#' The CBP distribution with parameters \eqn{b} and \eqn{\gamma} has pmf
#' \deqn{f(x|b,\gamma) = C \Gamma(ib+x) \Gamma(-ib+x) / (\Gamma(\gamma+x) x!), x=0,1,2,...} 
#' where \eqn{i} is the imaginary unit, \eqn{\Gamma(·)} the gamma function and 
#' \deqn{C = \Gamma(\gamma-ib) \Gamma(\gamma+ib) / (\Gamma(\gamma) \Gamma(ib) \Gamma(-ib))}
#' the normalizing constant.
#' 
#' The CBP is a particular case of the CTP when \eqn{a=0}.
#'
#' The mean and the variance of the CBP distribution are
#' \eqn{E(X)=\mu=b^2/(\gamma-1)} and \eqn{Var(X)=\mu(\mu+\gamma-1)/(\gamma-2)}
#' so \eqn{\gamma > 2}.
#'
#' It is always overdispersed.
#'
#' @return 
#' \code{dcbp} gives the pmf, \code{pcbp} gives the distribution function, \code{qcbp} gives the quantile function and \code{rcbp} generates random values.
#'
#' @references 
#' 
#' \insertRef{RCS2003}{cpd}
#' 
#' @seealso
#' Probability mass function, distribution function, quantile function and random generation for the CTP distribution: \code{\link{dctp}}, \code{\link{pctp}}, \code{\link{qctp}} and \code{\link{rctp}}.
#' Functions for maximum-likelihood fitting of the CBP distribution: \code{\link{fitcbp}}.
#'
#' @examples
#' # Examples for the function dcbp
#' dcbp(3,2,5)
#' dcbp(c(3,4),2,5)
#' 
#' @name cbp
#'
#' @rdname cbp
#' @importFrom hypergeo complex_gamma
#' @export
#'
#'

dcbp <- function(x, b, gamma) {

  if ( mode(c(x,b,gamma)) != "numeric")
    stop( "non-numeric argument to mathematical function" )

  if( gamma <= 0)
    stop("gamma must be greater than 0")

  x <- as.vector(x)
  errors<-c()


  for ( i in 1:length(x) ){
    if ( round(x[[i]]) != x[[i]] || x[[i]]<0){
      warning( paste ("non-integer or negative x[[", i, "]]", sep = ""))
      errors <- c(errors,i)
    }

  }

  if (b==0){
    lpmf<-2*(lgamma(gamma)+lgamma(x))-lgamma(gamma)-lgamma(gamma+x)-lgamma(x+1)
  } else {
    i <- sqrt(as.complex(-1))
    lpmf<-Re(complex_gamma(gamma-b*i,log=TRUE)+complex_gamma(gamma+b*i,log=TRUE)+complex_gamma(x+b*i,log=TRUE)+complex_gamma(x-b*i,log=TRUE)) -lgamma(gamma)-Re(complex_gamma(b*i,log=TRUE)+complex_gamma(b*i,log=TRUE))-lgamma(gamma+x)-lgamma(x+1)
  }
  result <- exp(lpmf)
  for (i in 1:length(errors)){
    result[ errors[[i]] ] = 0
  }
  return(result)
}


#' The Extended Biparametric Waring (EBW) Distribution
#'
#' @description
#' Probability mass function, distribution function, quantile function and random generation for the Extended Biparametric Waring (EBW) distribution with parameters \eqn{\alpha} and \eqn{\gamma} (or \eqn{\rho}). 
#'
#' @usage
#' debw(x, alpha, gamma, rho)
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param alpha parameter alpha (real)
#' @param rho parameter rho (positive)
#' @param gamma parameter \eqn{\gamma} (positive)
#' @param lower.tail if TRUE (default), probabilities are \eqn{P(X<=x)}, otherwise, \eqn{P(X>x)}.
#'
#' @details
#' The EBW distribution with parameters \eqn{\alpha} and \eqn{\gamma} has pmf
#' \deqn{f(x|a,\alpha,\gamma) = C \Gamma(\alpha+x)^2 / (\Gamma(\gamma+x) x!), x=0,1,2,...} 
#' where \eqn{\Gamma(·)} is the gamma function and 
#' \deqn{C = \Gamma(\gamma-\alpha^2 / (\Gamma(\alpha)^2 \Gamma(\gamma-2a))}
#' the normalizing constant.
#' 
#' There is an alternative parametrization in terms of \eqn{\alpha} and \eqn{\rho=\gamma-2\alpha>0} 
#' when \eqn{\alpha>0}. So, introduce only \eqn{\alpha} and \eqn{\gamma} or \eqn{\alpha} and \eqn{\rho},
#' depending on the parametrization you wish to use.
#' 
#' 
#' The mean and the variance of the EBW distribution are
#' \eqn{E(X)=\mu=\alpha^2/(\gamma-2\alpha-1)} and \eqn{Var(X)=\mu(\mu+\gamma-1)/(\gamma-2\alpha-2)}
#' so \eqn{\gamma > 2a + 2}.
#'
#' It is underdispersed if \eqn{\alpha < - (\mu + 1) / 2}, equidispersed if \eqn{\alpha = - (\mu + 1) / 2} or overdispersed
#' if \eqn{\alpha > - (\mu + 1) / 2}. In particular, if \eqn{\alpha >= -0.5} the EBW is overdispersed, whereas if 
#' \eqn{\alpha < -1} the EBW is underdispersed. In the case \eqn{-1 < \alpha <= -0.5}, the EBW may be under-, equi- or 
#' overdispersed depending on the value of \eqn{\gamma}.
#' 
#'
#' @return 
#' \code{debw} gives the pmf, \code{pebw} gives the distribution function, \code{qebw} gives the quantile function and \code{rebw} generates random values.																																	
#'
#' If \eqn{\alpha > 0} the probability mass function, distribution function, quantile function and random generation function for the UGW\eqn{(\alpha,\alpha,\rho)} distribution arise.
#' 
#' If \eqn{\alpha < 0} the probability mass function, distribution function, quantile function and random generation function for the CTP\eqn{(\alpha,0,\gamma)} distribution arise.
#' 
#' @references 
#' \insertRef{RCS2003}{cpd}
#' 
#' \insertRef{RCSO2004}{cpd}
#' 
#' \insertRef{ROC2018}{cpd}
#'
#' @seealso
#' Functions for maximum-likelihood fitting of the CTP and CBP distributions: \code{\link{fitctp}} and \code{\link{fitcbp}}.
#'
#' @examples
#' # Examples for the function dctp
#' debw(3,1,rho=5)
#' debw(c(3,4),2,rho=5)
#'
#' @name ebw
#'
#' @rdname ebw
#' @importFrom hypergeo complex_gamma
#' @export
#'

debw<-function(x,alpha,gamma,rho){
  
  if (!missing(gamma) & !missing(rho))
    stop("Specify only 'gamma' or 'rho'")
  
  if (missing(gamma) & missing(rho))
    stop("Specify 'gamma' or 'rho'")
  
  if ( !((missing(rho) && (mode(c(x,alpha,gamma)) == "numeric")) | 
         (missing(gamma) && (mode(c(x,alpha,rho)) == "numeric"))))
    stop( "non-numeric argument to mathematical function" )
  
  if (!missing(gamma)){

    if (alpha>0){
      warning("The usual parametrization when alpha>0 is ('alpha','rho')")
      if (gamma <= 2*alpha)
        stop("gamma must be greater than 2a")
    } else{
      if (gamma <= 0)
        stop("gamma must be positive")  
    }
  }
  
  if (!missing(rho)){
    if (alpha<0)
      stop("In the parametrization ('alpha','rho'), 'alpha' must be positive")
    if (rho<=0)
      stop (("rho must be positive") )
  }
  
  x <- as.vector(x)
  errors<-c()
  
  for ( i in 1:length(x) ){
    if ( round(x[[i]]) != x[[i]] || x[[i]]<0){
      warning( paste ("non-integer or negative x[[", i, "]]", sep = ""))
      errors <- c(errors,i)
    }
    
  }
  
  if (alpha>0 && !missing(rho)){
    lpmf<-2*lgamma(rho+alpha)+2*lgamma(x+alpha)-lgamma(rho)-2*lgamma(alpha)-lgamma(rho+2*alpha+x)-lgamma(x+1)
  }else{
    lpmf<-2*lgamma(gamma-alpha)+2*lgamma(x+alpha)-lgamma(gamma-2*alpha)-2*lgamma(alpha)-lgamma(gamma+x)-lgamma(x+1)
  }
  
  result <- exp(lpmf)
  for (i in 1:length(errors)){
    result[ errors[[i]] ] = 0
  }
  return(result)
}

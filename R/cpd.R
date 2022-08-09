#' @title
#' The Complex Triparametric Pearson (CTP) Distribution
#' @name
#' cpd
#' @details
#' The Complex Triparametric Pearson (CTP) distribution with parameters \eqn{a}, \eqn{b} and \eqn{\gamma} has pmf
#' \deqn{f(x|a,b,\gamma) = C \Gamma(a+ib+x) \Gamma(a-ib+x) / (\Gamma(\gamma+x) x!), x=0,1,2,...} 
#' where \eqn{i} is the imaginary unit, \eqn{\Gamma(·)} the gamma function and 
#' \deqn{C = \Gamma(\gamma-a-ib) \Gamma(\gamma-a+ib) / (\Gamma(\gamma-2a) \Gamma(a+ib) \Gamma(a-ib))}
#' the normalizing constant.
#' 
#' If \eqn{a=0} the CTP is a Complex Biparametric Pearson (CBP) distribution, so the pmf of the CBP distribution is obtained.
#' 
#' If \eqn{b=0} the CTP is an Extended Biparametric Waring (EBW) distribution, so the pmf of the EBW distribution is obtained.
#' In this case, \eqn{a} is call \eqn{\alpha}.
#'
#' The mean and the variance of the CTP distribution are
#' \eqn{E(X)=\mu=(a^2+b^2)/(\gamma-2a-1)} and \eqn{Var(X)=E(X)·(E(X)+\gamma-1)/(\gamma-2a-2)}
#' so \eqn{\gamma>2a+2}.
#'
#' It is underdispersed if \eqn{a<-(\mu+1)/2}, equidispersed if \eqn{a=-(\mu+1)/2} or overdispersed
#' if \eqn{a>-(\mu+1)/2}. In particular, if \eqn{a>=0} the CTP is always overdispersed.
#' 
#' @importFrom Rdpack reprompt
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
#' @keywords internal
"_PACKAGE"
#> [1] "_PACKAGE"

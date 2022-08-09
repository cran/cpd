#' @rdname ctp
#' @importFrom hypergeo complex_gamma
#' @export
#'
#' @examples
#' # Examples for the function pctp
#' pctp(3,1,2,3)
#' pctp(c(3,4),1,2,3)
#' 

pctp <- function(q, a, b, gamma, lower.tail = TRUE) {
  if ( mode(c(q,a,b,gamma)) != "numeric")
    stop( "non-numeric argument to mathematical function" )

  if( gamma <= 2 * a )
    stop( "gamma must be greater than 2a" )

  q<-as.vector(q)
  maxX=q[1]
  n<-length(q)
  for ( i in 1:n ){
    q[i] = floor( q[i] )
    if (q[i] > maxX)
      maxX=q[i]
  }
  icomplex <- sqrt(as.complex(-1))
  lf0 <- 2 * Re(complex_gamma(gamma - a + b * icomplex, log = TRUE)) - lgamma(gamma) - lgamma(gamma - 2 * a)
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  
  #Generating distribution function
  digits=options()$digit
  while( i <= maxX+1 && Fd[i]<(1-10^-digits)){
    pmfAux <- exp(log(pmfAux)+log(((a+i-1)^2+b^2))-log((gamma+i-1))-log(i))
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    i <- i + 1
  }

  result<-vector(mode="numeric",length=n)
  for ( i in 1:n ){
    if ( q[i] < 0 )
      result[i]=0
    else
      result[i]=Fd[q[i]+1]

    if (! lower.tail){
      result[i]<-1-result[i]
    }
  }

  return (result)
}

#' @rdname cbp
#' @importFrom hypergeo complex_gamma
#' @export
#'
#'
#' @examples
#' # Examples for the function pcbp
#' pcbp(3,2,3)
#' pcbp(c(3,4),2,3)
#' 

pcbp <- function(q, b, gamma, lower.tail = TRUE ) {
  if ( mode(c(q,b,gamma)) != "numeric")
    stop( "non-numeric argument to mathematical function" )

  if( gamma <= 0 )
    stop( "gamma must be greater than 0" )

  q<-as.vector(q)
  maxX=q[1]
  n<-length(q)
  for ( i in 1:n ){
    q[i] = floor( q[i] )
    if (q[i] > maxX)
      maxX=q[i]
  }
  icomplex <- sqrt(as.complex(-1))
  lf0 <- 2 * ( Re(complex_gamma(gamma  + b * icomplex, log = TRUE)) - lgamma(gamma))
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  
  #Generating distribution function
  while( i <= maxX+1 ){
    pmfAux <- exp(log(pmfAux)+log(((i-1)^2+b^2))-log((gamma+i-1))-log(i))
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    i <- i + 1
  }

  result<-vector(mode="numeric",length=n)
  for ( i in 1:n ){
    if ( q[i] < 0 )
      result[i]=0
    else
      result[i]=Fd[q[i]+1]

    if (! lower.tail){
      result[i]<-1-result[i]
    }
  }

  return (result)
}

#' @rdname ebw
#' @importFrom hypergeo complex_gamma
#' @export
#'
#'
#' @examples
#' # Examples for the function pebw
#' pebw(3,2,rho=5)
#' pebw(c(3,4),2,rho=5)
#' 

pebw <- function(q,alpha,gamma,rho, lower.tail = TRUE ) {
  if (!missing(gamma) & !missing(rho))
    stop("Specify only 'gamma' or 'rho'")
  
  if (missing(gamma) & missing(rho))
    stop("Specify 'gamma' or 'rho'")
  
  if ( !((missing(rho) && (mode(c(q,alpha,gamma)) == "numeric")) | 
         (missing(gamma) && (mode(c(q,alpha,rho)) == "numeric"))))
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
  
  if (alpha>0 && !missing(rho)){
    auxgamma=2*alpha+rho
  }else{
    auxgamma=gamma
  }

  q<-as.vector(q)
  maxX=q[1]
  n<-length(q)
  for ( i in 1:n ){
    q[i] = floor( q[i] )
    if (q[i] > maxX)
      maxX=q[i]
  }

  lf0 <- 2 * lgamma(auxgamma-alpha)-lgamma(auxgamma-2*alpha)-lgamma(auxgamma)
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  
  #Generating distribution function
  while( i <= maxX+1 ){
    pmfAux <- exp(log(pmfAux)+log((alpha+i-1)^2)-log(auxgamma+i-1)-log(i))
    #pmfAux <- pmfAux * (alpha+i-1)^2 / ((auxgamma+i-1) *i)
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    i <- i + 1
  }
  
  result<-vector(mode="numeric",length=n)
  for ( i in 1:n ){
    if ( q[i] < 0 )
      result[i]=0
    else
      result[i]=Fd[q[i]+1]
    
    if (! lower.tail){
      result[i]<-1-result[i]
    }
  }
  
  return (result)
  
}
#' @rdname ctp
#' @importFrom hypergeo complex_gamma
#' @export
#'
#' @examples
#' # Examples for the function rctp
#' rctp(10,1,1,3)
#' 

rctp<-function(n, a, b, gamma){
  if ( mode(c(n,a,b,gamma)) != "numeric")
    stop( "non-numeric argument to mathematical function")

  if( (gamma <= 2 * a) || (gamma <= 0) )
    stop("gamma must be greater than max(0,2a)")

  icomplex <- sqrt(as.complex(-1))

  randoms <- runif(n, 0, 1)
  maxRandoms <- max(randoms)
  result<-vector(mode="numeric",length=n)

  lf0 <- 2 * Re(complex_gamma(gamma - a + b * icomplex, log = TRUE)) - lgamma(gamma) - lgamma(gamma - 2 * a)
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  
  #Generating distribution function
  while( Fd[[i]] < maxRandoms ){
    pmfAux <- exp(log(pmfAux)+log(((a+i-1)^2+b^2))-log((gamma+i-1))-log(i))
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    i <- i + 1
  }
  #Searching values
  for (i in 1:length(randoms)){
    pMin=1
    pMax=length(Fd)
    while (Fd[pMin] < randoms[i]){
      mitad = floor ((pMin + pMax)/2)
      if (Fd[mitad] <= randoms[i] && Fd[mitad+1] >= randoms[i])
        pMin = mitad + 1
      else if (Fd[mitad] <= randoms[i])
          pMin = mitad
      else
          pMax = mitad
    }

    result[[i]]=pMin-1
  }
  return (result)
}


#' @rdname cbp
#' @importFrom hypergeo complex_gamma
#' @export
#'
#' @examples
#' # Examples for the function rcbp
#' rcbp(10,1,3)
#' 

rcbp<-function(n, b, gamma){
  if ( mode(c(n,b,gamma)) != "numeric")
    stop( "non-numeric argument to mathematical function")

  if( gamma <= 0)
    stop("gamma must be greater than 0")

  icomplex <- sqrt(as.complex(-1))

  randoms <- runif(n, 0, 1)
  maxRandoms <- max(randoms)
  result<-vector(mode="numeric",length=n)

  lf0 <- 2 * Re(complex_gamma(gamma + b * icomplex, log = TRUE)) - lgamma(gamma) - lgamma(gamma)
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  
  #Generating distribution function
  while( Fd[[i]] < maxRandoms ){
    pmfAux <- exp(log(pmfAux)+log(((i-1)^2+b^2))-log((gamma+i-1))-log(i))
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    i <- i + 1
  }
  #Searching values
  for (i in 1:length(randoms)){
    pMin=1
    pMax=length(Fd)

    while (Fd[pMin] < randoms[i]){
      mitad = floor ((pMin + pMax)/2)
      if (Fd[mitad] <= randoms[i] && Fd[mitad+1] >= randoms[i])
        pMin = mitad + 1
      else if (Fd[mitad] <= randoms[i])
        pMin = mitad
      else
        pMax = mitad
    }

    result[[i]]=pMin-1
  }
  return (result)
}



#' @rdname ebw
#' @importFrom hypergeo complex_gamma
#' @export
#'
#'
#' @examples
#' # Examples for the function rebw
#' rebw(10,2,rho=5)
#' rebw(10,-2.1,gamma=5)
#' 

rebw <- function(n,alpha,gamma,rho, lower.tail = TRUE ) {
  if ( !((missing(rho) && (mode(c(n,alpha,gamma)) == "numeric")) | 
         (missing(gamma) && (mode(c(n,alpha,rho)) == "numeric"))))
    stop( "non-numeric argument to mathematical function" )
  if (n<0)
    stop( "n must be postive integer")
  else if (trunc(n)==0)
    return (integer(0))
  
  if (!missing(gamma) & !missing(rho))
    stop("Specify only 'gamma' or 'rho'")
  
  if (missing(gamma) & missing(rho))
    stop("Specify 'gamma' or 'rho'")
  
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
  
  randoms <- runif(n, 0, 1)
  maxRandoms <- max(randoms)
  result<-vector(mode="numeric",length=n)
  
  lf0 = 2 * lgamma(auxgamma-alpha)-lgamma(auxgamma-2*alpha)-lgamma(auxgamma)
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  
  #Generating distribution function
  while( Fd[[i]] < maxRandoms ){
    pmfAux <- exp(log(pmfAux)+log((alpha+i-1)^2)-log(auxgamma+i-1)-log(i))
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    i <- i + 1
  }
  #Searching values
  for (i in 1:length(randoms)){
    pMin=1
    pMax=length(Fd)
    
    while (Fd[pMin] < randoms[i]){
      mitad = floor ((pMin + pMax)/2)
      if (Fd[mitad] <= randoms[i] && Fd[mitad+1] >= randoms[i])
        pMin = mitad + 1
      else if (Fd[mitad] <= randoms[i])
        pMin = mitad
      else
        pMax = mitad
    }
    
    result[[i]]=pMin-1
  }
  return (result)
}
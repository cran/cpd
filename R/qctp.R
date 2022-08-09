#' @rdname ctp
#' @importFrom hypergeo complex_gamma
#' @export
#'
#' @examples
#' # Examples for the function qctp
#' qctp(0.5,1,2,3)
#' qctp(c(.8,.9),1,2,3)
#' 

qctp <- function(p, a, b, gamma, lower.tail = TRUE ){
  if ( mode(c(p,a,b,gamma)) != "numeric")
    stop( "non-numeric argument to mathematical function")

  if( gamma <= 2 * a )
    stop("gamma must be greater than 2a")

  icomplex <- sqrt(as.complex(-1))
  auxP=p[p<1 & p>0]
  if (length(auxP)>0)
     if (lower.tail){
       maxP <- max(auxP)
     }else{
       maxP <- 1-min(auxP)
     }
  else
     maxP=0
  
  n=length(p)
  result<-vector(mode="numeric",length=n)

  lf0 <- 2 * Re(complex_gamma(gamma - a + b * icomplex, log = TRUE)) - lgamma(gamma) - lgamma(gamma - 2 * a)
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  
  #Generating distribution function
  while( Fd[[i]] < maxP ){
    pmfAux <- exp(log(pmfAux)+log(((a+i-1)^2+b^2))-log((gamma+i-1))-log(i))
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    i <- i + 1
  }
  #Searching values
  for (i in 1:n){
    pMin=1
    pMax=length(Fd)
    
    if (! lower.tail){
      p[i]<-1-p[i]
    }
    
    
    if (p[i]>1 || p[i]<0){
      warning( paste ("p[[", i, "]] must be a probability", sep = ""))
      result[[i]]=NaN
    }else if (p[i]==1){
      result[[i]]=Inf
    }else{
      while (Fd[pMin] < p[i]){
        mitad = floor ((pMin + pMax)/2)
        if (Fd[mitad] <= p[i] && Fd[mitad+1] >= p[i])
          pMin = mitad + 1
        else if (Fd[mitad] <= p[i])
          pMin = mitad
        else
          pMax = mitad
      }

      result[[i]]=pMin-1
    }
  }
  return (result)
}


#' @rdname cbp
#' @importFrom hypergeo complex_gamma
#' @export
#'
#' @examples
#' # Examples for the function qcbp
#' qcbp(0.5,2,3)
#' qcbp(c(.8,.9),2,3)
#' 
qcbp <- function(p, b, gamma, lower.tail = TRUE)  {

  if ( mode(c(p,b,gamma)) != "numeric")
    stop( "non-numeric argument to mathematical function")

  if( gamma <= 0 )
    stop("gamma must be greater than 0")

  icomplex <- sqrt(as.complex(-1))

  auxP=p[p<1 & p>0]
  if (length(auxP)>0)
    if (lower.tail){
      maxP <- max(auxP)
    }else{
      maxP <- 1-min(auxP)
    }
  else
    maxP=0
  n=length(p)
  result<-vector(mode="numeric",length=n)

  lf0 <- 2 * Re(complex_gamma(gamma  + b * icomplex, log = TRUE)) - lgamma(gamma) - lgamma(gamma )
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  
  #Generating distribution function
  while( Fd[[i]] < maxP ){
    pmfAux <- exp(log(pmfAux)+log(((i-1)^2+b^2))-log((gamma+i-1))-log(i))
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    i <- i + 1
  }
  #Searching values
  for (i in 1:n){
    pMin=1
    pMax=length(Fd)

    if (! lower.tail){
      p[i]<-1-p[i]
    }
    
    if (p[i]>1 || p[i]<0){
      warning( paste ("p[[", i, "]] must be a probability", sep = ""))
      result[[i]]=NaN
    }else if (p[i]==1){
      result[[i]]=Inf
    }else{
      while (Fd[pMin] < p[i]){
        mitad = floor ((pMin + pMax)/2)
        if (Fd[mitad] <= p[i] && Fd[mitad+1] >= p[i])
          pMin = mitad + 1
        else if (Fd[mitad] <= p[i])
          pMin = mitad
        else
          pMax = mitad
      }

      result[[i]]=pMin-1
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
#' # Examples for the function qebw
#' qebw(0.5,-2.1,gamma=0.1)
#' qebw(c(.8,.9),-2.1,gamma=0.1)
#' qebw(0.5,2,rho=5)
#' qebw(c(.8,.9),2,rho=5)
#' 

qebw <- function(p,alpha,gamma,rho, lower.tail = TRUE ) {
  if (!missing(gamma) & !missing(rho))
    stop("Specify only 'gamma' or 'rho'")
  
  if (missing(gamma) & missing(rho))
    stop("Specify 'gamma' or 'rho'")
  
  if ( !((missing(rho) && (mode(c(p,alpha,gamma)) == "numeric")) | 
         (missing(gamma) && (mode(c(p,alpha,rho)) == "numeric"))))
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

  auxP=p[p<1 & p>0]
  if (length(auxP)>0)
    if (lower.tail){
      maxP <- max(auxP)
    }else{
      maxP <- 1-min(auxP)
    }
  else
    maxP=0
  n=length(p)
  result<-vector(mode="numeric",length=n)
  
  lf0  <- 2 * lgamma(auxgamma-alpha)-lgamma(auxgamma-2*alpha)-lgamma(auxgamma)
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  
  #Generating distribution function
  while( Fd[[i]] < maxP ){
    pmfAux <- exp(log(pmfAux)+log((alpha+i-1)^2)-log(auxgamma+i-1)-log(i))
    #pmfAux <- pmfAux * (alpha+i-1)^2 / ((auxgamma+i-1) *i)
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    i <- i + 1
  }
  #Searching values
  for (i in 1:n){
    pMin=1
    pMax=length(Fd)
    
    if (! lower.tail){
      p[i]<-1-p[i]
    }
    
    if (p[i]>1 || p[i]<0){
      warning( paste ("p[[", i, "]] must be a probability", sep = ""))
      result[[i]]=NaN
    }else if (p[i]==1){
      result[[i]]=Inf
    }else{
      while (Fd[pMin] < p[i]){
        mitad = floor ((pMin + pMax)/2)
        if (Fd[mitad] <= p[i] && Fd[mitad+1] >= p[i])
          pMin = mitad + 1
        else if (Fd[mitad] <= p[i])
          pMin = mitad
        else
          pMax = mitad
      }
      
      result[[i]]=pMin-1
    }
  }
  return (result)
}

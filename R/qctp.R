#' @rdname ctp
#' @importFrom fAsianOptions cgamma
#' @export
#'
#' @examples
#' # Examples for the function qctp
#' qctp(0.5,1,2,3)
#' qctp(c(.8,.9),1,2,3)

qctp <- function(p, a, b, gamma, lower.tail = TRUE ){
  if ( mode(c(p,a,b,gamma)) != "numeric")
    stop( "non-numeric argument to mathematical function")

  if( gamma <= 2 * a )
    stop("gamma must be greater than 2a")

  icomplex <- sqrt(as.complex(-1))
  auxP=p[p!=1]
  if (length(auxP)>0)
      maxP <- max(auxP)
  else
     maxP=0
  n=length(p)
  result<-vector(mode="numeric",length=n)

  lf0 <- 2 * Re(cgamma(gamma - a + b * icomplex, log = TRUE)) - lgamma(gamma) - lgamma(gamma - 2 * a)
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  #Generating Density Distribution
  while( Fd[[i]] < maxP ){
    pmfAux <- exp(log(pmfAux)+log(((a+i-1)^2+b^2))-log((gamma+i-1))-log(i))
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    i <- i + 1
  }
  #Searching values
  for (i in 1:n){
    pMin=1
    pMax=length(Fd)

    if (p[i]==1){
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
#' @importFrom fAsianOptions cgamma
#' @export
#'
#' @export
#'
#' @examples
#' # Examples for the function qcbp
#' qcbp(0.5,2,3)
#' qcbp(c(.8,.9),2,3)
qcbp <- function(p, b, gamma, lower.tail = TRUE)  {

  if ( mode(c(p,b,gamma)) != "numeric")
    stop( "non-numeric argument to mathematical function")

  if( gamma <= 0 )
    stop("gamma must be greater than 0")

  icomplex <- sqrt(as.complex(-1))

  auxP=p[p!=1]
  if (length(auxP)>0)
    maxP <- max(auxP)
  else
    maxP=0
  n=length(p)
  result<-vector(mode="numeric",length=n)

  lf0 <- 2 * Re(cgamma(gamma  + b * icomplex, log = TRUE)) - lgamma(gamma) - lgamma(gamma )
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  #Generating Density Distribution
  while( Fd[[i]] < maxP ){
    pmfAux <- exp(log(pmfAux)+log(((i-1)^2+b^2))-log((gamma+i-1))-log(i))
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    i <- i + 1
  }
  #Searching values
  for (i in 1:n){
    pMin=1
    pMax=length(Fd)

    if (p[i]==1){
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

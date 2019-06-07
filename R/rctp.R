#' @rdname ctp
#' @importFrom fAsianOptions cgamma
#' @export
#'
#' @examples
#' # Examples for the function rctp
#' rctp(4,1,1,3)

rctp<-function(n, a, b, gamma){
  if ( mode(c(n,a,b,gamma)) != "numeric")
    stop( "non-numeric argument to mathematical function")

  if( gamma <= 2 * a )
    stop("gamma must be greater than 2a")

  icomplex <- sqrt(as.complex(-1))

  randoms <- runif(n, 0, 1)
  maxRandoms <- max(randoms)
  result<-vector(mode="numeric",length=n)

  lf0 <- 2 * Re(cgamma(gamma - a + b * icomplex, log = TRUE)) - lgamma(gamma) - lgamma(gamma - 2 * a)
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  #Generating Density Distribution
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
#' @importFrom fAsianOptions cgamma
#' @export
#'
#' @export
#'
#' @examples
#' # Examples for the function rcbp(4,1,3)

rcbp<-function(n, b, gamma){
  if ( mode(c(n,b,gamma)) != "numeric")
    stop( "non-numeric argument to mathematical function")

  if( gamma <= 0)
    stop("gamma must be greater than 0")

  icomplex <- sqrt(as.complex(-1))

  randoms <- runif(n, 0, 1)
  maxRandoms <- max(randoms)
  result<-vector(mode="numeric",length=n)

  lf0 <- 2 * Re(cgamma(gamma + b * icomplex, log = TRUE)) - lgamma(gamma) - lgamma(gamma)
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  #Generating Density Distribution
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


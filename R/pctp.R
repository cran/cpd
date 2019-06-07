#' @rdname ctp
#' @importFrom fAsianOptions cgamma
#' @export
#'
#' @examples
#' # Examples for the function pctp
#' pctp(3,1,2,3)
#' pctp(c(3,4),1,2,3)

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
  lf0 <- 2 * Re(cgamma(gamma - a + b * icomplex, log = TRUE)) - lgamma(gamma) - lgamma(gamma - 2 * a)
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  #Generating distribution function
  digits=options()$digit
  while( i <= maxX+1 && Fd[i]<(1-10^-digits)){
    pmfAux <- exp(log(pmfAux)+log(((a+i-1)^2+b^2))-log((gamma+i-1))-log(i))
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    print(paste(i,Fd[i],pmfAux,sep=":"))
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
#' @importFrom fAsianOptions cgamma
#' @export
#'
#'
#' @examples
#' # Examples for the function pcbp
#' pcbp(3,2,3)
#' pcbp(c(3,4),2,3)

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
  lf0 <- 2 * Re(cgamma(gamma  + b * icomplex, log = TRUE)) - lgamma(gamma) - lgamma(gamma )
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

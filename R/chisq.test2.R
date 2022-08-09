#' @name chisq.test2
#' @title Pearson's Chi-squared Test for Count Data
#' @description \code{chisq.test2} performs Pearson chi-squared goodness-of-fit test for count data
#' @param obs a numeric vector with the counts
#' @param p.esp a numeric vector with the expected probabilities of the same length of \code{obs}. They must sum 1.
#' @param npar an integer specifying the number of parameters of the model. By default \code{npar} is \code{NULL}, so the degrees of freedom of
#' the chi-squared statistics are the number of classes minus 1.
#' @param grouping a logical indicating whether to group in classes with expected frequency greater than or equal to 5. By default \code{grouping} is \code{FALSE}.
#' @return A list with class \code{"htest"} containing the following components:
#' \itemize{
#' \item \code{statistic}: the value of the chi-squared test statistic.
#' \item \code{df}: the degrees of freedom of the approximate chi-squared distribution.
#' \item \code{p.value}: the p-value for the test.
#' \item \code{observed}: the observed counts.
#' \item \code{observed.grouped}: the observed counts grouped in classes with expected frequency greather or equal to 5.
#' \item \code{expected}: the expected counts under the null hypothesis.
#' \item \code{expected.grouped}: the expected counts under the null hypothesis grouped in classes with expected frequency greather or equal to 5.
#' \item \code{residuals}: the Pearson residuals computed as \code{(observed - expected) / sqrt(expected)}.
#' }
#' @importFrom stats pchisq
#' @export
#' 
#' @examples
#' set.seed(123)
#' x <- rctp(500, -1.5, 1, 2)
#' table(x)
#' obs <- c(172, 264, 57, 6, 0, 1)
#' fit <- fitctp(x)
#' p.esp <- c(dctp(0:(length(obs)-1),fit$coefficients[1],fit$coefficients[2],
#'            fit$coefficients[3])[1:(length(obs)-1)],1-sum(dctp(0:(length(obs)-1),
#'            fit$coefficients[1],fit$coefficients[2],fit$coefficients[3])[1:(length(obs)-1)]))
#' chisq.test2(obs, p.esp)
#' chisq.test2(obs, p.esp, grouping = TRUE)
#' chisq.test2(obs, npar= 3, p.esp)

chisq.test2 <- function (obs, p.esp, npar = NULL, grouping = FALSE) {
  
  if (length(obs) != length(p.esp)) 
    stop("'obs' and 'p.esp' must have the same length")
  if (any(obs < 0) || anyNA(obs)) 
    stop("all entries of 'obs' must be nonnegative and finite")
  if ((n <- sum(obs)) == 0) 
    stop("at least one entry of 'obs' must be positive")
  if (sum(p.esp) < 0.9999999999999999) 
    stop("expected probabilities must sum 1")
  
  n <- sum(obs)
  m <- length(obs)
  
  esp <- n * p.esp
  
  nesp <- c()
  nobs <- c()
  i <- 1
  j <- 1
  suma.e <- 0
  suma.o <- 0
  
  if(grouping == TRUE){
    METHOD <- "Pearson's Chi-squared test grouping in classes with expected frequency >= 5"
    while(i <= length(esp)){
      suma.e <- suma.e + esp[i]
      suma.o <- suma.o + obs[i]
      if(suma.e >= 5){
        nesp[j] <- suma.e
        nobs[j] <- suma.o
        j <- j+1
        suma.e <- 0
        suma.o <- 0
      }
      i <- i+1
    }
    if (suma.e > 0){
      nesp[j-1] <- nesp[j-1] + suma.e
      nobs[j-1] <- nobs[j-1] + suma.o
    }
  }
  else{
    METHOD <- "Pearson's Chi-squared test"
    nesp <- esp
    nobs <- obs
  }
  
  if(is.null(npar)){
    DF <- length(nesp) - 1
  }
  else DF <- length(nesp) - npar - 1
  if(DF <= 0) DF <- NaN
  
  RESIDUALS <- (nobs - nesp) / sqrt(nesp)
  
  STATISTIC <- sum(RESIDUALS ^ 2)
  
  PVALUE <- pchisq(STATISTIC, DF, lower.tail = FALSE)
  
  names(STATISTIC) <- "X-squared"
  
  names(DF) <- "df"
  
  if (any(esp < 5) && is.finite(DF) && isFALSE(grouping)) 
    warning("Chi-squared approximation may be incorrect (expected frequencies less than 5)")
  
  ans <- list(statistic = STATISTIC, df = DF, p.value = PVALUE, observed = obs, observed.grouped = nobs,
                 expected = esp, expected.grouped = nesp, residuals = RESIDUALS, method = METHOD)
  
  class(ans) <- "htest"
  
  ans
}
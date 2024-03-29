% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitctp.R
\name{fitebw}
\alias{fitebw}
\title{Maximum-likelihood fitting of the EBW distribution}
\usage{
fitebw(x, alphastart = NULL, rhostart = NULL, gammastart = NULL, 
          method = "L-BFGS-B", control = list(),...)
}
\arguments{
\item{x}{A numeric vector of length at least one containing only finite values.}

\item{alphastart}{A starting value for the parameter \eqn{\alpha}; by default \code{NULL}.}

\item{rhostart}{A starting value for the parameter \eqn{\rho>0}; by default \code{NULL}.}

\item{gammastart}{A starting value for the parameter \eqn{\gamma>max(0,2\alpha)}; by default \code{NULL}.}

\item{method}{The method to be used in fitting the model. The default method is "L-BFGS-B" (optim).}

\item{control}{A list of parameters for controlling the fitting process.}

\item{...}{Additional parameters.}
}
\value{
An object of class \code{'fitEBW'} is a list containing the following components:

\itemize{
\item \code{n}, the number of observations,
\item \code{initialValues}, a vector with the starting values used,
\item \code{coefficients}, the parameter ML estimates of the CTP distribution,
\item \code{se}, a vector of the standard error estimates,
\item \code{hessian}, a symmetric matrix giving an estimate of the Hessian at the solution found in the optimization of the log-likelihood function,
\item \code{cov}, an estimate of the covariance matrix of the model coefficients,
\item \code{corr}, an estimate of the correlation matrix of the model estimates,
\item \code{loglik}, the maximized log-likelihood,
\item \code{aic}, Akaike Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters,
\item \code{bic}, Bayesian Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters,
\item \code{code}, a code that indicates successful convergence of the fitter function used (see nlm and optim helps),
\item \code{converged},  logical value that indicates if the optimization algorithms succesfull.
\item \code{method}, the name of the fitter function used.
}

Generic functions:

\itemize{
\item \code{print}: The print of a \code{'fitEBW'} object shows the ML parameter estimates and their standard errors in parenthesis.
\item \code{summary}: The summary provides the ML parameter estimates, their standard errors and the statistic and p-value of the Wald test to check if the parameters are significant.
This summary also shows the loglikelihood, AIC and BIC values, as well as the results for the chi-squared goodness-of-fit test and the Kolmogorov-Smirnov test for discrete variables. Finally, the correlation matrix between parameter estimates appears.
\item \code{coef}: It extracts the fitted coefficients from a \code{'fitEBW'} object.
\item \code{logLik}: It extracts the estimated log-likelihood from a \code{'fitEBW'} object.
\item \code{AIC}: It extracts the value of the Akaike Information Criterion from a \code{'fitEBW'} object.
\item \code{BIC}: It extracts the value of the Bayesian Information Criterion from a \code{'fitEBW'} object.
\item \code{plot}: It shows the plot of a \code{'fitEBW'} object. Observed and theoretical probabilities, empirical and theoretical cumulative distribution functions or empirical cumulative probabilities against theoretical cumulative probabilities are the three plot types.
}
}
\description{
Maximum-likelihood fitting of the Extended Biparametric Waring (EBW) distribution with parameters \eqn{\alpha}, \eqn{\rho} and \eqn{\gamma}. Generic
methods are \code{print}, \code{summary}, \code{coef}, \code{logLik}, \code{AIC}, \code{BIC} and \code{plot}. The method to be used in fitting the 
model is "L-BFGS-B" which allows constraints for each variable (see details in \code{\link{optim}} funtion).
}
\details{
If the starting value for \eqn{\alpha} is positive, the parameterization \eqn{(\alpha,\rho)} is used; 
otherwise, the parameterization \eqn{(\alpha,\gamma)} is used.

If the starting values of the parameters \eqn{\alpha}, \eqn{\gamma} or \eqn{\rho} are omitted (default option), 
they are computing by the method of moments (if possible; otherwise they must be entered).

The default method is \code{"L-BFGS-B"} (see details in \code{\link{optim}} function),
but non-linear minimization (\code{\link{nlm}}) or those included in the \code{optim} function 
(\code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"} and \code{"SANN"}) may be used.

Standard error (SE) estimates for \eqn{\alpha}, \eqn{\gamma} or \eqn{\rho} are provided by the default method; 
otherwise, SE for \eqn{\alpha_0} and \eqn{\gamma_0} where \eqn{\alpha=-exp(\alpha_0)} and \eqn{\gamma=exp(\gamma_0)}
(or for \eqn{\alpha_0} and \eqn{\rho_0} where \eqn{\alpha=exp(\alpha_0)} and \eqn{\rho=exp(\rho_0)}) are computed.
}
\examples{
set.seed(123)
x <- rebw(500, 2, rho = 5)
fitebw(x)
summary(fitebw(x, alphastart = 1, rhostart = 5))
}
\references{
\insertRef{COR2021}{cpd}
}
\seealso{
Plot of observed and theoretical frequencies for a EBW fit: \code{\link{plot.fitEBW}}

Maximum-likelihood fitting for the CTP distribution: \code{\link{fitctp}}.

Maximum-likelihood fitting for the CBP distribution: \code{\link{fitcbp}}.
}

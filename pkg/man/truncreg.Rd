\name{truncreg}
\alias{truncreg}
\alias{print.truncreg}
\alias{fitted.truncreg}
\alias{residuals.truncreg}
\alias{summary.truncreg}
\alias{print.summary.truncreg}
\alias{vcov.truncreg}
\alias{logLik.truncreg}
\title{Models for truncated regressions}
\description{
  Estimation of models with truncated explanatory variables by maximum likelihood
}
\usage{
truncreg(formula, data, subset, weights, na.action,
       point = 0, direction = "left", ...)
\method{print}{truncreg}(x, digits = max(3, getOption("digits") - 2),
    width = getOption("width"), ...)
\method{summary}{truncreg}(object, ...)
\method{print}{summary.truncreg}(x, digits = max(3, getOption("digits") - 2),
    width = getOption("width"), ...)
\method{logLik}{truncreg}(object, ...)
\method{vcov}{truncreg}(object, ...)
\method{residuals}{truncreg}(object, ...)
\method{fitted}{truncreg}(object, ...)
}

\arguments{
  \item{x, object}{an object of class \code{truncreg}}
  \item{formula}{a symbolic description of the model to be estimated,}
  \item{data}{the data,}
  \item{subset}{an optional vector specifying a subset of observations,}
  \item{weights}{an optional vector of weights,}
  \item{na.action}{a function which indicates what should happen when
    the data contains '\code{NA}'s,}
  \item{point}{the value of truncation (the default is 0),}
  \item{direction}{the direction of the truncation, either \code{"left"}
    (the default) or \code{"right"},}
  \item{digits}{the number of digits,}
  \item{width}{the width of the printing,}

  \item{...}{further arguments.}
  }

\value{
  An object of class \code{"truncreg"}, a list with elements:

  \item{coefficients}{the named vector of coefficients,}
  \item{vcov}{the variance matrix of the coefficients,}
  \item{fitted.values}{the fitted values,}
  \item{logLik}{the value of the log-likelihood,}
  \item{gradient}{the gradient of the log-likelihood at convergence,}
  \item{model}{the model frame used,}
  \item{call}{the matched call,}
  \item{est.stat}{some information about the estimation (time used,
    optimisation method),}
}
\details{
  
  The model is estimated with the \code{maxLik} package and the
  Newton-Raphson method, using analytic gradient and hessian.
}

\examples{

## Simulate a data.frame
n <- 10000
sigma <- 4
alpha <- 2
beta <- 1
x <- rnorm(n,0,2)
eps <- rnorm(n)
y <- alpha+beta*x+eps*sigma
d <- data.frame(y = y, x = x)

## Use a truncated subsample
dl1 <- subset(d, y>1)

## Use truncreg to estimate consistently the model

truncreg(y~x, dl1, point = 1, direction = "left")
}

\author{Yves Croissant}

\references{

  Hausman, J.A. and D.A. Wise (1976) ``The evaluation of results from
  truncated samples: the New-Jersey negative invome tax experiment'',
  \emph{Annals of Economic ans Social Measurment}, 5, pp.421--45.

  Hausman, J.A. and D.A. Wise (1976)  ``Social experimentation,
  truncated distributions and efficient estimation'',
  \emph{Econometrica}, 45, pp.421--5.
}
\keyword{regression}

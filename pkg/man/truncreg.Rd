\name{truncreg}
\alias{truncreg}
\alias{coef.truncreg}
\alias{fitted.truncreg}
\alias{logLik.truncreg}
\alias{model.frame.truncreg}
\alias{model.matrix.truncreg}
\alias{predict.truncreg}
\alias{print.truncreg}
\alias{residuals.truncreg}
\alias{summary.truncreg}
\alias{print.summary.truncreg}
\alias{vcov.truncreg}

\title{Truncated Gaussian Response Models}

\description{
  Estimation of models for truncated Gaussian variables by maximum likelihood.
}

\usage{
truncreg(formula, data, subset, weights, na.action,
  point = 0, direction = "left",
  model = TRUE, y = FALSE, x = FALSE, ...)
}

\arguments{
  \item{formula}{a symbolic description of the model to be estimated,}
  \item{data}{the data,}
  \item{subset}{an optional vector specifying a subset of observations,}
  \item{weights}{an optional vector of weights,}
  \item{na.action}{a function which indicates what should happen when
    the data contains '\code{NA}'s,}
  \item{point}{the value of truncation (the default is 0),}
  \item{direction}{the direction of the truncation, either \code{"left"}
    (the default) or \code{"right"},}
  \item{model, y, x}{logicals. If \code{TRUE} the corresponding components
    of the fit (model frame, response, model matrix) are returned.}
  \item{...}{further arguments.}
}

\value{
  An object of class \code{"truncreg"}, a list with elements:

  \item{coefficients}{the named vector of coefficients,}
  \item{vcov}{the variance matrix of the coefficients,}
  \item{fitted.values}{the fitted values,}
  \item{logLik}{the value of the log-likelihood,}
  \item{gradient}{the gradient of the log-likelihood at convergence,}
  \item{nobs}{the number of observations,}
  \item{call}{the matched call,}
  \item{terms}{the model terms,}
  \item{model}{the model frame used (if \code{model = TRUE}),}
  \item{y}{the response vector (if \code{y = TRUE}),}
  \item{x}{the model matrix (if \code{x = TRUE}),}
  \item{point}{the truncation point used,}
  \item{direction}{the truncation direction used,}
  \item{est.stat}{some information about the estimation (time used,
    optimization method),}
}

\details{
  The model is estimated with the \code{\link[maxLik]{maxLik}} package and the
  Newton-Raphson method, using analytic gradient and Hessian.

  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"truncreg"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[stats]{coef}}, 
  \code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}}, 
  \code{\link[stats]{predict}}, \code{\link[stats]{fitted}}, \code{\link[stats]{model.frame}}, 
  and \code{\link[stats]{model.matrix}}.
}

\seealso{\code{\link[maxLik]{maxLik}}, \code{\link[mhurdle]{mhurdle}}}

\references{
  Cragg JG (1971). Some Statistical Models for Limited Dependent Variables
    with Application to the Demand for Durable Goods.
    \emph{Econometrica}, \bold{39}, 829--844.

  Hausman JA, Wise DA (1976). The Evaluation of Results from Truncated Samples:
    The New-Jersey Negative Income Tax Experiment.
    \emph{Annals of Economic ans Social Measurment}, \bold{5}, 421--445.

  Hausman JA, Wise DA (1976). Social Experimentation,
    Truncated Distributions and Efficient Estimation.
    \emph{Econometrica}, \bold{45}, 421--425.

  Tobin J (1958). Estimation of Relationships for Limited Dependent Variables.
    \emph{Econometrica}, \bold{26}, 24--36.
}

\examples{
########################
## Artificial example ##
########################

## simulate a data.frame
set.seed(1071)
n <- 10000
sigma <- 4
alpha <- 2
beta <- 1
x <- rnorm(n, mean = 0, sd = 2)
eps <- rnorm(n, sd = sigma)
y <- alpha + beta * x + eps
d <- data.frame(y = y, x = x)

## truncated response
d$yt <- ifelse(d$y > 1, d$y, NA)

## binary threshold response
d$yb <- factor(d$y > 0)

## censored response
d$yc <- pmax(1, d$y)

## compare estimates for full/truncated/censored/threshold response
fm_full <- lm(y ~ x, data = d)
fm_trunc <- truncreg(yt ~ x, data = d, point = 1, direction = "left")
fm_thresh <- glm(yb ~ x, data = d, family = binomial(link = "probit"))
library("survival")
fm_cens <- survreg(Surv(yc, yc > 1, type = "left") ~ x, data = d, dist = "gaussian")

## compare scaled regression coefficients
cbind(
  "True"      = c(alpha, beta) / sigma,
  "Full"      = coef(fm_full) / summary(fm_full)$sigma,
  "Truncated" = coef(fm_trunc)[1:2] / coef(fm_trunc)[3],
  "Censored"  = coef(fm_cens) / fm_cens$scale,
  "Threshold" = coef(fm_thresh)
)


################################
## Tobin's durable goods data ##
################################

## Tobit model (Tobin 1958)
data("tobin", package = "survival")
tobit <- survreg(Surv(durable, durable > 0, type = "left") ~ age + quant,
  data = tobin, dist = "gaussian")

## Two-part model (Cragg 1971)
## (see "mhurdle" package for a combined solution)
cragg_probit <- glm(factor(durable > 0) ~ age + quant,
  data = tobin, family = binomial(link = "logit"))
cragg_trunc <- truncreg(durable ~ age + quant, data = tobin, subset = durable > 0)

## Scaled coefficients
cbind(
  "Tobit"     = coef(tobit) / tobit$scale,
  "Binary"    = coef(cragg_probit),
  "Truncated" = coef(cragg_trunc)[1:3] / coef(cragg_trunc)[4])

## likelihood ratio test and BIC
ll <- c("Tobit" = tobit$loglik[1],
        "Two-Part" = as.vector(logLik(cragg_probit) + logLik(cragg_trunc)))
df <- c(4, 3 + 4)
pchisq(2 * diff(ll), diff(df), lower.tail = FALSE)
-2 * ll + log(nrow(tobin)) * df
}

\keyword{regression}

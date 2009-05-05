ml.truncreg <- function(param, X, y, gradient = FALSE, hessian = FALSE, fit = FALSE, point, direction){
  beta <- param[1:ncol(X)]
  sigma <- param[length(param)]
  bX <- as.numeric(crossprod(beta,t(X)))
  resid <- y-bX
  if (direction == "left"){
    trunc <- bX-point
    sgn <- 1
  }
  else{
    trunc <- point-bX
    sgn <- -1
  }
  mills <- dnorm(trunc/sigma)/pnorm(trunc/sigma)

  lnL <- sum(log(dnorm(resid/sigma))-log(sigma)
               -log(pnorm(trunc/sigma)))
  if (gradient){
    gbX <- resid/sigma^2 -sgn/sigma*mills
    gsigma <- resid^2/sigma^3-1/sigma+trunc/sigma^2*mills
    gradi <- cbind(gbX*X,as.numeric(gsigma))
    attr(lnL,"gradient") <- gradi
  }
  if (fit){
    attr(lnL,"fit") <- bX
  }
  if (hessian){
    bb <- mills*(trunc/sigma+mills)/sigma^2-1/sigma^2
    ss <- -3*resid^2/sigma^4+1/sigma^2+trunc^2/sigma^4*mills*(mills+trunc/sigma)-
      2*trunc*mills/sigma^3
    bs <- -2*resid/sigma^3+sgn*(mills/sigma^2-trunc/sigma^3*mills*(mills+trunc/sigma))
    bb <- crossprod(bb*X,X)
    bs <- apply(bs*X,2,sum)
    ss <- sum(ss)
    h <- rbind(cbind(bb,bs),c(bs,ss))
    attr(lnL,"hessian") <- h
  }
  lnL
}

truncreg <- function(formula, data, subset, weights, na.action, point = 0, direction = "left", ...){
  formula.type <- TRUE
  if (class(formula[[3]]) == "name"){
    X <- try(eval(formula[[3]],sys.frame(which=-3)),silent = TRUE)
    if (class(X) == "matrix") formula.type <- FALSE else formula.type <- TRUE
  }
  
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  
  if (formula.type){
    m <- match(c("formula", "data", "subset", "na.action","weights"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    X <- model.matrix(formula, data = mf)
    y <- model.response(mf)
  }
  else{
    y <- eval(formula[[2]],sys.frame(which=-3))
  }
  result <- truncreg.fit(X, y, point, direction, ...)
  result$model <- mf
  result$call <- cl
  result
}

truncreg.fit <- function(X, y, point, direction, ...){
  dots <- list(...)
  if (is.null(dots$method)) method <- "nr" else method <- dots$method
  if (is.null(dots$iterlim)) iterlim <- 50 else iterlim <- dots$iterlim
  if (is.null(dots$print.level)) print.level <- 0 else print.level <- dots$print.level
  
  oldoptions <- options(warn=-1)
  on.exit(options(oldoptions))
  start.time <- proc.time()

  f <- function(param)  ml.truncreg(param,  X = X, y = y,
                                    gradient = FALSE, hessian = FALSE,
                                    fit = FALSE, point = point,
                                    direction = direction)
  g <- function(param){
    attr(ml.truncreg(param, X = X, y = y,
                     gradient = TRUE, hessian = FALSE,
                     fit = FALSE, point = point,
                     direction = direction),"gradient")
  } 
  h <- function(param){
    attr(ml.truncreg(param, X = X, y = y,
                     gradient = FALSE, hessian = TRUE,
                     fit = FALSE, point = point,
                     direction = direction),"hessian")
  } 

  linmod <- lm(y~X-1)
  start <- c(coef(linmod),summary(linmod)$sigma)
  maxl <- maxLik(f, g, h, start = start, method = method,
                 iterlim = iterlim, print.level = print.level)
  grad.conv <- g(maxl$estimate)
  coefficients <- maxl$estimate
  vcov <- -solve(maxl$hessian)
  fit <- attr(ml.truncreg(coefficients, X = X, y = y,
                          gradient = FALSE, hessian = FALSE,
                          fit = TRUE, point = point,
                          direction = direction),"fit")
  logLik <- maxl$maximum
  attr(logLik,"df") <- length(coefficients)
  hessian <- maxl$hessian
  convergence.OK <- maxl$code<=2
  elaps.time <- proc.time() - start.time
  nb.iter <- maxl$iterations
  eps <- maxl$gradient%*%solve(-maxl$hessian)%*%maxl$gradient

  est.stat <- list(elaps.time = elaps.time,
                   nb.iter = nb.iter,
                   eps = eps,
                   method = maxl$type,
                   message = maxl$message
                   )
  class(est.stat) <- "est.stat"
  coef.names <- c(colnames(X),"sigma")
  names(coefficients) <- rownames(vcov) <- colnames(vcov) <- coef.names

  result <- list(coefficients = coefficients,
                 vcov = vcov,
                 fitted.values = fit,
                 logLik = logLik,
                 gradient = grad.conv,
                 model = NULL,
                 call = NULL,
                 est.stat = est.stat
                 )
  class(result) <- c("truncreg","maxLik")
  result
}

fitted.truncreg <- function(object, ...){
  object$fitted.values
}

residuals.truncreg <- function(object, ...){
  model.frame(object)[[1]]-fitted(object)
}

coef.truncreg <- function(object, ...){
  object$coefficients
}

vcov.truncreg <- function(object, ...){
  object$vcov
}

logLik.truncreg <- function(object, ...){
  x <- object$logLik
  attr(x,"df") <- NULL
  x
}

print.truncreg <- function (x, digits = max(3, getOption("digits") - 2), width = getOption("width"), ...){
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

summary.truncreg <- function (object,...){
  b <- coef(object)
  std.err <- sqrt(diag(vcov(object)))
  z <- b/std.err
  p <- 2*(1-pnorm(abs(z)))
  CoefTable <- cbind(b,std.err,z,p)
  colnames(CoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
  object$CoefTable <- CoefTable
  class(object) <- c("summary.truncreg","truncreg")
  return(object)
}

print.summary.truncreg <- function(x,digits= max(3, getOption("digits") - 2),width=getOption("width"),...){
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  cat("\nCoefficients :\n")
  printCoefmat(x$CoefTable,digits=digits)
  cat("\n")
  df <- attr(x$logLik,"df")
  cat(paste("Log-Likelihood: ",
            signif(x$logLik,digits),
            " on ",df," Df\n",sep=""))
  invisible(x)
}

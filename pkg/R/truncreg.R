ml.truncreg <- function(param, X, y, gradient = FALSE, hessian = FALSE, fit = FALSE, point, direction, scaled){
    beta <- param[1:ncol(X)]
    sigma <- param[length(param)]
    bX <- as.numeric(crossprod(beta, t(X)))
    sgn <- ifelse(direction == "left", +1, -1)
    if (!scaled){
        resid <- (y - bX)
        trunc <- (bX - point)
        # Update of the mills function, use logs to avoid Inf
        #mills <- dnorm(sgn * trunc / sigma) / pnorm(sgn * trunc / sigma)
        mills <- exp(dnorm(sgn * trunc / sigma, log = TRUE) - pnorm(sgn * trunc / sigma, log.p = TRUE))
        lnL <-  - log(sigma) + dnorm(resid / sigma, log = TRUE) - pnorm(sgn * trunc / sigma, log.p = TRUE)
        if (gradient){
            gbX <- resid / sigma ^ 2 - sgn * mills / sigma
            gsigma <- - 1 / sigma + resid ^ 2 / sigma ^ 3  + sgn * mills * trunc / sigma ^ 2
            gradi <- cbind(gbX * X, as.numeric(gsigma))
            attr(lnL, "gradient") <- gradi
        }
        if (fit) attr(lnL, "fit") <- bX
        if (hessian){
            bb <- -1 / sigma ^ 2 + mills * (sgn * trunc / sigma + mills) / sigma ^ 2
            ss <- 1 / sigma ^ 2 - 3 * resid ^ 2 / sigma ^ 4 - 2 * mills * sgn * trunc / sigma ^ 3 +
                mills * (sgn * trunc / sigma + mills) * trunc / sigma ^ 3
            ss <- 1 / sigma ^ 2 - 3 * resid ^ 2 / sigma ^ 4 - 2 * mills * sgn * trunc / sigma ^ 3 +
                mills * (sgn * trunc / sigma + mills) * trunc^2 / sigma ^ 4
            
            bs <- - 2 * resid / sigma ^ 3 + sgn * mills / sigma ^ 2 -
                mills * (mills + sgn * trunc / sigma) * trunc / sigma ^ 3
            bb <- crossprod(bb * X, X)
            bs <- apply(bs * X, 2, sum)
            ss <- sum(ss)
            h <- rbind(cbind(bb, bs), c(bs, ss))
            attr(lnL,"hessian") <- h
        }
    }
    else{
        lnL <- - log(sigma) + dnorm(y / sigma - bX, log = TRUE) - pnorm(sgn * (bX - point / sigma), log.p = TRUE)
        # Update of the mills function, use logs to avoid Inf (YC 2015/12/11)
        # mills <- dnorm(sgn * (bX - point / sigma)) / pnorm(sgn * (bX - point / sigma))
        mills <- exp(dnorm(sgn * (bX - point / sigma), log = TRUE) - pnorm(sgn * (bX - point / sigma), log.p = TRUE))
        if (gradient){
            gbX <- (y / sigma - bX) - mills * sgn
            gsigma <- - 1 / sigma + (y / sigma - bX) * y / sigma ^ 2 - sgn * mills * point / sigma ^ 2
            gradi <- cbind(gbX * X, as.numeric(gsigma))
            attr(lnL, "gradient") <- gradi
        }
        if (fit) attr(lnL, "fit") <- bX * sigma
        if(hessian){
            bb <- -1 + mills * (mills + sgn * (bX - point / sigma))
            bs <- - y / sigma ^ 2 + (mills + sgn * (bX - point / sigma)) * mills * point / sigma ^ 2
            ss <- 1 / sigma ^ 2 - 3 * y ^ 2 / sigma ^ 4 + 2 * bX * y / sigma ^ 3 +
                mills * (mills + sgn * (bX - point / sigma)) * point ^ 2 / sigma ^ 4 +
                    2 * sgn * mills * point / sigma ^ 3
            bb <- crossprod(bb * X, X)
            bs <- apply(bs * X, 2, sum)
            ss <- sum(ss)
            h <- rbind(cbind(bb, bs), c(bs, ss))
            attr(lnL,"hessian") <- h
        }
    }
    lnL
}

truncreg <- function(formula, data, subset, weights, na.action, point = 0, direction = "left",
  model = TRUE, y = FALSE, x = FALSE, scaled = FALSE, ...){
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
        Y <- model.response(mf)
        mt <- attr(mf, "terms")
    }
    else{
        Y <- eval(formula[[2]], sys.frame(which=-3))
        mt <- terms(formula)
    }
    
     ## process options
    direction <- match.arg(direction, c("left", "right"))
    point <- rep(point, length.out = length(Y))
    
    result <- truncreg.fit(X, Y, point, direction, scaled, ...)
    result$call <- cl
    result$terms <- mt
    if(model) result$model <- mf
    if(y) result$y <- Y
    if(x) result$x <- X
    result
}

truncreg.fit <- function(X, y, point, direction, scaled, ...){
    ## check input
    if(direction == "left" & any(y < point)) stop("response not truncated, contains observations < 'point'")
    if(direction == "right" & any(y > point)) stop("response not truncated, contains observations > 'point'")

    dots <- list(...)
    if (is.null(dots$method)) method <- "bfgs" else method <- dots$method
    if (is.null(dots$iterlim)) iterlim <- 50 else iterlim <- dots$iterlim
    if (is.null(dots$print.level)) print.level <- 0 else print.level <- dots$print.level
    
    oldoptions <- options(warn = - 1)
    on.exit(options(oldoptions))
    start.time <- proc.time()

    f <- function(param)  ml.truncreg(param,  X = X, y = y,
                                      gradient = TRUE, hessian = TRUE,
                                      fit = FALSE, point = point,
                                      direction = direction, scaled = scaled)
    linmod <- lm(y ~ X - 1)
    start <- c(coef(linmod), summary(linmod)$sigma)
    if (scaled) start[1:ncol(X)] <- start[1:ncol(X)] / start[ncol(X) + 1]
    
    if (FALSE){
        f0 <-  ml.truncreg(start,  X = X, y = y,
                           gradient = TRUE, hessian = FALSE,
                           fit = FALSE, point = point, direction = direction, scaled = scaled)
        agrad <- apply(attr(f0, "gradient"), 2, sum)
        ostart <- start
        of <- sum(f(start))
        ngrad <- c()
        eps <- 1E-5
        for (i in 1:length(start)){
            start <- ostart
            start[i] <- start[i] + eps
            ngrad <- c(ngrad, (sum(f(start)) - of) / eps)
        }
        start <- ostart
        print(cbind(agrad, ngrad))
        stop()
    }
    
    maxl <- maxLik(f, start = start, method = method,
                   iterlim = iterlim, print.level = print.level)
    coefficients <- maxl$estimate
    
    if (scaled){
        coefficients[1:ncol(X)] <- coefficients[1:ncol(X)] * coefficients[ncol(X) + 1]
        f <- function(param)  ml.truncreg(param,  X = X, y = y,
                                          gradient = TRUE, hessian = TRUE,
                                          fit = FALSE, point = point,
                                          direction = direction, scaled = FALSE)
        maxl <- maxLik(f, start = coefficients, method = method,
                       iterlim = 0, print.level = print.level)
    }
    coefficients <- maxl$estimate
    grad.conv <- maxl$gradient
    vcov <- - solve(maxl$hessian)
    fit <- attr(ml.truncreg(coefficients, X = X, y = y,
                            gradient = TRUE, hessian = TRUE,
                            fit = TRUE, point = point,
                            direction = direction, scaled = scaled), "fit")
    names(fit) <- rownames(X)
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
                   gradientObs = maxl$gradientObs,
                   nobs = length(y),
                   call = NULL,
                   terms = NULL,
                   model = NULL,
                   y = NULL,
                   x = NULL,
                   point = if(isTRUE(all.equal(rep(point[1], length(point)), point))) point[1] else point,
                   direction = direction,
                   est.stat = est.stat
                   )
    class(result) <- c("truncreg", "maxLik")
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

logLik.truncreg <- function(object, ...)
    structure(object$logLik, df = length(object$coefficients), nobs = object$nobs, class = "logLik")

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
    object$coefficients <- CoefTable
    class(object) <- c("summary.truncreg","truncreg")
    return(object)
}

print.summary.truncreg <- function(x,digits= max(3, getOption("digits") - 2),width=getOption("width"),...){
    cat("\nCall:\n")
    print(x$call)
    if (!is.null(x$est.stat)){
        cat("\n")
        print(x$est.stat)
    }
    
    cat("\n")
    cat("\nCoefficients :\n")
    printCoefmat(x$coefficients,digits=digits)
    cat("\n")
    df <- attr(x$logLik,"df")
    cat(paste("Log-Likelihood: ",
              signif(x$logLik,digits),
              " on ",df," Df\n",sep=""))
    invisible(x)
}

model.frame.truncreg <- function(formula, ...) {
    if(!is.null(formula$model)) return(formula$model)
    NextMethod()
}

model.matrix.truncreg <- function(object, ...)
    if(!is.null(object$x)) object$x else model.matrix(object$terms, model.frame(object), ...)


predict.truncreg <- function(object, newdata = NULL, na.action = na.pass, ...) 
{
    if(missing(newdata)) {
        rval <- object$fitted.values
    } else {
        mt <- delete.response(object$terms)
        X <- model.matrix(mt, model.frame(mt, newdata, na.action = na.action))
        rval <- drop(X %*% head(object$coefficients, -1))
    }
    return(rval)
}

print.est.stat <- function(x, ...){
    et <- x$elaps.time[3]
    i <- x$nb.iter[1]
    halton <- x$halton
    method <- x$method
    if (!is.null(x$type) && x$type != "simple"){
        R <- x$nb.draws
        cat(paste("Simulated maximum likelihood with", R, "draws\n"))
    }
    s <- round(et,0)
    h <- s%/%3600
    s <- s-3600*h
    m <- s%/%60
    s <- s-60*m
    cat(paste(method, "method\n"))
    tstr <- paste(h, "h:", m, "m:", s, "s", sep="")
    cat(paste(i,"iterations,",tstr,"\n"))
    if (!is.null(halton)) cat("Halton's sequences used\n")
    if (!is.null(x$eps)) cat(paste("g'(-H)^-1g =", sprintf("%5.3G", as.numeric(x$eps)),"\n"))
    if (is.numeric(x$code)){
        msg <- switch(x$code,
                      "1" = "gradient close to zero",
                      "2" = "successive function values within tolerance limits",
                      "3" = "last step couldn't find higher value",
                      "4" = "iteration limit exceeded"
                      )
        cat(paste(msg, "\n"))
    }
    else cat(paste(x$code, "\n"))
}

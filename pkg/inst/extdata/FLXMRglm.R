# Patch provided by Bettina Grün (flexmix author).
# https://cran.r-project.org/web/packages/flexmix/index.html

FLXMRglm <- function(formula=.~., family=gaussian, offset=NULL)
{
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  glmrefit <- function(x, y, w) {
    fit <- c(glm.fit(x, y, weights=w, offset=offset, family=family),
      list(call = sys.call(), offset = offset,
        control = eval(formals(glm.fit)$control),
        method = "weighted.glm.fit"))
    fit$df.null <- sum(w) + fit$df.null - fit$df.residual - fit$rank
    fit$df.residual <- sum(w) - fit$rank
    fit$x <- x
    fit
  }

  z <- new("FLXMRglm", weighted=TRUE, formula=formula,
    name=paste("FLXMRglm", family$family, sep=":"), offset = offset,
    family=family$family, refit=glmrefit)

  z@preproc.y <- function(x) {
    if (ncol(x) > 1)
      stop(paste("for the", family$family, "family y must be univariate"))
    x
  }

  if (family$family=="gaussian") {
    z@defineComponent <- function(para) {
      predict <- function(x, ...) {
        dotarg = list(...)
        if("offset" %in% names(dotarg)) offset <- dotarg$offset
        p <- x %*% para$coef
        if (!is.null(offset)) p <-  p + offset
        family$linkinv(p)
      }

      logLik <- function(x, y, ...)
        dnorm(y, mean=predict(x, ...), sd=para$sigma, log=TRUE)

      new("FLXcomponent",
        parameters=list(coef=para$coef, sigma=para$sigma),
        logLik=logLik, predict=predict,
        df=para$df)
    }

    z@fit <- function(x, y, w, component){
      fit <- glm.fit(x, y, w=w, offset=offset, family = family)
      z@defineComponent(para = list(coef = coef(fit), df = ncol(x)+1,
        sigma =  sqrt(sum(fit$weights * fit$residuals^2 /
        mean(fit$weights))/ (nrow(x)-fit$rank))))
    }
  }

  else if (family$family=="binomial") {
    z@preproc.y <- function(x) {
      if (ncol(x) != 2)
      {
        stop("for the binomial family, y must be a 2 column matrix\n",
          "where col 1 is no. successes and col 2 is no. failures")
      }
      if (any(x < 0))
        stop("negative values are not allowed for the binomial family")
      x
    }

    z@defineComponent <- function(para) {
      predict <- function(x, ...) {
        dotarg = list(...)
        if("offset" %in% names(dotarg))
          offset <- dotarg$offset
        p <- x %*% para$coef
        if (!is.null(offset))
          p <- p + offset
        family$linkinv(p)
      }
      logLik <- function(x, y, ...)
      dbinom(y[,1], size=rowSums(y), prob=predict(x, ...), log=TRUE)

      new("FLXcomponent",
        parameters=list(coef=para$coef),
        logLik=logLik, predict=predict,
        df=para$df)
    }

    z@fit <- function(x, y, w, component) {
      fit <- glm.fit(x, y, weights=w, family=family, offset=offset, start=component$coef)
      z@defineComponent(para = list(coef = coef(fit), df = ncol(x)))
    }
  }

  else if (family$family=="poisson") {
    z@defineComponent <- function(para) {
      predict <- function(x, ...) {
        dotarg = list(...)
        if("offset" %in% names(dotarg)) offset <- dotarg$offset
        p <- x %*% para$coef
        if (!is.null(offset)) p <- p + offset
        family$linkinv(p)
      }
      logLik <- function(x, y, ...)
        dpois(y, lambda=predict(x, ...), log=TRUE)

      new("FLXcomponent",
        parameters=list(coef=para$coef),
        logLik=logLik, predict=predict,
        df=para$df)
    }

    z@fit <- function(x, y, w, component) {
      fit <- glm.fit(x, y, weights=w, family=family, offset=offset, start=component$coef)
      z@defineComponent(para = list(coef = coef(fit), df = ncol(x)))
    }
  }

  else if (family$family=="Gamma") {
    z@defineComponent <- function(para) {
      predict <- function(x, ...) {
        dotarg = list(...)
        if("offset" %in% names(dotarg)) offset <- dotarg$offset
        p <- x %*% para$coef
        if (!is.null(offset)) p <- p + offset
        family$linkinv(p)
      }
      logLik <- function(x, y, ...)
        dgamma(y, shape = para$shape, scale=predict(x, ...)/para$shape, log=TRUE)

      new("FLXcomponent",
        parameters = list(coef = para$coef, shape = para$shape),
        predict = predict, logLik = logLik,
        df = para$df)
    }

    z@fit <- function(x, y, w, component) {
      fit <- glm.fit(x, y, weights=w, family=family, offset=offset, start=component$coef)
      z@defineComponent(para = list(coef = coef(fit), df = ncol(x)+1,
        shape = sum(fit$prior.weights)/fit$deviance))
    }
  }

  else stop(paste("Unknown family", family))
  z
}

# extractParam
#
# Extract successive values of a projection of the parameter(s).
# The method works both on a list of lists of results,
# or on a single list of parameters matrices.
#
# @inheritParams plotHist
#
.extractParam <- function(mr, x=1, y=1)
{
  if (is.list(mr[[1]])) {
    # Obtain L vectors where L = number of res lists in mr
    return ( lapply( mr, function(mr_list) {
      sapply(mr_list, function(m) m[x,y])
    } ) )
  }
  sapply(mr, function(m) m[x,y])
}

#' plotHist
#'
#' Plot compared histograms of a single parameter (scalar)
#'
#' @name plotHist
#'
#' @param mr Output of multiRun(), list of lists of functions results
#' @param x Row index of the element inside the aggregated parameter
#' @param y Column index of the element inside the aggregated parameter
#' @param ... Additional graphical parameters (xlab, ylab, ...)
#'
#' @examples
#' # mr[[i]] is a list of estimated parameters matrices (here random matrices).
#' # Should be mr <- multiRun(...) --> see bootstrap example in ?multiRun.
#' simmr_path <- system.file("extdata", "simulateMr.R", package = "morpheus")
#' source(simmr_path)
#' mr <- simulateMr(c(2,2), 10)$mr
#' plotHist(mr, 2, 1) #second row, first column
#'
#' @return No return value, called for side effects.
#'
#' @export
plotHist <- function(mr, x, y, ...)
{
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  params <- .extractParam(mr, x, y)
  L <- length(params)
  # Plot histograms side by side
  par(mfrow=c(1,L), cex.axis=1.5, cex.lab=1.5, mar=c(4.7,5,1,1))
  args <- list(...)
  for (i in 1:L) {
    hist(params[[i]], breaks=40, freq=FALSE,
      xlab=ifelse("xlab" %in% names(args), args$xlab, "Parameter value"),
      ylab=ifelse("ylab" %in% names(args), args$ylab, "Density"))
  }
  NULL
}

# NOTE: roxygen2 bug, "@inheritParams plotHist" fails in next header:

#' plotBox
#'
#' Draw compared boxplots of a single parameter (scalar)
#'
#' @name plotBox
#'
#' @param mr Output of multiRun(), list of lists of functions results
#' @param x Row index of the element inside the aggregated parameter
#' @param y Column index of the element inside the aggregated parameter
#' @param ... Additional graphical parameters (xlab, ylab, ...)
#'
#' @examples
#' # mr[[i]] is a list of estimated parameters matrices (here random matrices).
#' # Should be mr <- multiRun(...) --> see bootstrap example in ?multiRun.
#' simmr_path <- system.file("extdata", "simulateMr.R", package = "morpheus")
#' source(simmr_path)
#' source(simmr_path)
#' mr <- simulateMr(c(2,2), 10)$mr
#' plotBox(mr, 2, 1) #second row, first column
#'
#' @return No return value, called for side effects.
#'
#' @export
plotBox <- function(mr, x, y, ...)
{
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  params <- .extractParam(mr, x, y)
  L <- length(params)
  # Plot boxplots side by side
  par(mfrow=c(1,L), cex.axis=1.5, cex.lab=1.5, mar=c(4.7,5,1,1))
  for (i in 1:L)
    boxplot(params[[i]], ...)
  NULL
}

#' plotCoefs
#'
#' Draw a graph of (averaged) coefficients estimations with their standard,
#' deviations ordered by mean values.
#' Note that the drawing does not correspond to a function; it is just a
#' convenient way to visualize the estimated parameters.
#'
#' @name plotCoefs
#'
#' @param mr List of parameters matrices
#' @param params True value of the parameters matrix
#' @param ... Additional graphical parameters
#'
#' @examples
#' # mr[[i]] is a list of estimated parameters matrices (here random matrices).
#' # Should be mr <- multiRun(...) --> see bootstrap example in ?multiRun.
#' simmr_path <- system.file("extdata", "simulateMr.R", package = "morpheus")
#' source(simmr_path)
#' mr_θ <- simulateMr(c(3,2), 10)
#' mr <- mr_θ$mr ; θ <- mr_θ$θ
#' plotCoefs(mr[[1]], θ)
#'
#' @return No return value, called for side effects.
#'
#' @export
plotCoefs <- function(mr, params, ...)
{
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  d <- nrow(mr[[1]])
  K <- ncol(mr[[1]])

  params_hat <- matrix(nrow=d, ncol=K)
  stdev <- matrix(nrow=d, ncol=K)
  for (x in 1:d) {
    for (y in 1:K) {
      estims <- .extractParam(mr, x, y)
      params_hat[x,y] <- mean(estims)
      # Another way to compute stdev: using distances to true params
#      stdev[x,y] <- sqrt( mean( (estims - params[x,y])^2 ) )
      # HACK remove extreme quantile in estims[[i]] before computing sd()
      stdev[x,y] <- sd(estims) #[ estims < max(estims) & estims > min(estims) ] )
    }
  }

  par(cex.axis=1.5, cex.lab=1.5, mar=c(4.7,5,1,1))
  params <- as.double(params)
  o <- order(params)
  avg_param <- as.double(params_hat)
  std_param <- as.double(stdev)
  args <- list(...)
  matplot(
    cbind(params[o],avg_param[o],
      avg_param[o]+std_param[o],avg_param[o]-std_param[o]),
    col=1, lty=c(1,5,3,3), type="l", lwd=2,
    xlab=ifelse("xlab" %in% names(args), args$xlab, "Parameter index"),
    ylab=ifelse("ylab" %in% names(args), args$ylab, "") )
  NULL
}

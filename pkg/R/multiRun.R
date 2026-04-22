#' multiRun
#'
#' Estimate N times some parameters, outputs of some list of functions.
#' This method is thus very generic, allowing typically bootstrap or
#' Monte-Carlo estimations of matrices μ or β.
#' Passing a list of functions opens the possibility to compare them on a fair
#' basis (exact same inputs). It's even possible to compare methods on some
#' deterministic design of experiments.
#'
#' @name multiRun
#'
#' @param fargs List of arguments for the estimation functions
#' @param estimParams List of nf function(s) to apply on fargs
#' @param packages Vector of packages to load on each node (default: morpheus)
#' @param prepareArgs Prepare arguments for the functions inside estimParams
#' @param N Number of runs
#' @param ncores Number of cores for parallel runs (<=1: sequential)
#' @param agg Aggregation method (default: lapply)
#' @param verbose TRUE to indicate runs + methods numbers
#'
#' @return A list of nf aggregates of N results (matrices).
#'
#' @examples
#' \donttest{
#' β <- matrix(c(1,-2,3,1),ncol=2)
#'
#' # This example requires a patch for the flexmix package written by Bettina Grün:
#' patch_path <- system.file("extdata", "FLXMRglm.R", package = "morpheus")
#' source(patch_path)
#'
#' # Bootstrap + computeMu, morpheus VS flexmix
#' io <- generateSampleIO(n=1000, p=1/2, β=β, b=c(0,0), "logit")
#' μ <- normalize(β)
#' res <- multiRun(list(X=io$X,Y=io$Y,K=2), list(
#'   # morpheus
#'   function(fargs) {
#'     ind <- fargs$ind
#'     computeMu(fargs$X[ind,], fargs$Y[ind], list(K=fargs$K))
#'   },
#'   # flexmix
#'   function(fargs) {
#'     ind <- fargs$ind
#'     K <- fargs$K
#'     dat <- as.data.frame( cbind(fargs$Y[ind],fargs$X[ind,]) )
#'     out <- refit( flexmix( cbind(V1, 1 - V1) ~ 0+., data=dat, k=K,
#'       model=FLXMRglm(family="binomial") ) )
#'     normalize( matrix(out@@coef[1:(ncol(fargs$X)*K)], ncol=K) )
#'   } ),
#'   packages = c("morpheus", "flexmix"),
#'   prepareArgs = function(fargs,index) {
#'     if (index == 1)
#'       fargs$ind <- 1:nrow(fargs$X)
#'     else
#'       fargs$ind <- sample(1:nrow(fargs$X),replace=TRUE)
#'     fargs
#'   }, N=10, ncores=1) #ncores=3
#' for (i in 1:2)
#'   res[[i]] <- alignMatrices(res[[i]], ref=μ, ls_mode="exact")
#'
#' # Monte-Carlo + optimParams from X,Y, morpheus VS flexmix
#' res <- multiRun(list(n=1000,p=1/2,β=β,b=c(0,0),link="logit"), list(
#'   # morpheus
#'   function(fargs) {
#'     K <- fargs$K
#'     μ <- computeMu(fargs$X, fargs$Y, list(K=fargs$K))
#'     o <- optimParams(fargs$X, fargs$Y, fargs$K, fargs$link, fargs$M)
#'     o$run(list(β=μ))$β
#'   },
#'   # flexmix
#'   function(fargs) {
#'     K <- fargs$K
#'     dat <- as.data.frame( cbind(fargs$Y,fargs$X) )
#'     out <- refit( flexmix( cbind(V1, 1 - V1) ~ ., data=dat, k=K,
#'       model=FLXMRglm(family="binomial") ) )
#'     sapply( seq_len(K), function(i)
#'       as.double( out@@components[[1]][[i]][2:(1+ncol(fargs$X)),1] ) )
#'   } ),
#'   packages = c("morpheus", "flexmix"),
#'   prepareArgs = function(fargs,index) {
#'     library(morpheus)
#'     io <- generateSampleIO(fargs$n, fargs$p, fargs$β, fargs$b, fargs$link)
#'     fargs$X <- io$X
#'     fargs$Y <- io$Y
#'     fargs$K <- ncol(fargs$β)
#'     fargs$link <- fargs$link
#'     fargs$M <- computeMoments(io$X,io$Y)
#'     fargs
#'   }, N=10, ncores=1) #ncores=3
#' for (i in 1:2)
#'   res[[i]] <- alignMatrices(res[[i]], ref=β, ls_mode="exact")}
#'
#' @export
multiRun <- function(fargs, estimParams, packages = c("morpheus"),
  prepareArgs = function(x,i) x, N=10, ncores=3, agg=lapply, verbose=FALSE)
{
  if (!is.list(fargs))
    stop("fargs: list")
  # No checks on fargs: supposedly done in estimParams[[i]]()
  if (!is.list(estimParams))
    estimParams = list(estimParams)
  # Verify that the provided parameters estimations are indeed functions
  lapply(seq_along(estimParams), function(i) {
    if (!is.function(estimParams[[i]]))
      stop("estimParams: list of function(fargs)")
  })
  if (!is.numeric(N) || N < 1)
    stop("N: positive integer")

  estimParamAtIndex <- function(index)
  {
    fargs <- prepareArgs(fargs, index)
    if (verbose)
      cat("Run ",index,"\n")
    lapply(seq_along(estimParams), function(i) {
      if (verbose)
        cat("   Method ",i,"\n")
      out <- estimParams[[i]](fargs)
      if (is.list(out))
        do.call(rbind, out)
      else
        out
    })
  }

  loadPackages <- function() {
    for (p in packages)
      library(p, character.only = TRUE)
  }

  if (ncores > 1)
  {
    cl <- parallel::makeCluster(ncores, outfile="")
    parallel::clusterExport(cl, c("fargs", "verbose", "loadPackages"), environment())
    parallel::clusterEvalQ(cl, loadPackages() )
    list_res <- parallel::clusterApplyLB(cl, 1:N, estimParamAtIndex)
    parallel::stopCluster(cl)
  }
  else {
    loadPackages()
    list_res <- lapply(1:N, estimParamAtIndex)
  }

  # De-interlace results: output one list per function
  nf <- length(estimParams)
  lapply( seq_len(nf), function(i) lapply(seq_len(N), function(j) list_res[[j]][[i]]) )
}

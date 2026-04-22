# simulateMr
#
# Simulate an output of multiRun(...), because this function is a little
# complicated to call in an example (see ?multiRun).
#
#' β <- matrix(c(1,-2,3,1),ncol=2)


# theta, not beta here (more general)
.simulateMr(
#' # mr[[i]] is a list of estimated parameters matrices (here random matrices).
#' # Should be mr <- multiRun(...) --> see bootstrap example in ?multiRun.
#' mr <- list()
#' μ <- normalize(β)
#' for (i in 1:2) {
#'   mr[[i]] <- list()
#'   for (j in 1:3)
#'     mr[[i]][[j]] <- β + matrix(rnorm(4,sd=0.25),ncol=2)
#'   mr[[i]] <- alignMatrices(mr[[i]], ref=β, ls_mode="exact") --> TODO: align in same function
#' }

            list(theta = ..., mr = ...)

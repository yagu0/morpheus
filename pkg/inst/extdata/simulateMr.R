# simulateMr
#
# Simulate an output of multiRun(...), because this function is a little
# complicated to call in an example (see ?multiRun).
#
# Usage: simulateMr(c(nrow, ncol), N)
# where nrow,ncol = dimension of the parameter and N = number of "runs"
#
# Always length(mr) == 2 and random parameter + noise, to not over-parameterize.
#
# Return a list with mr = simulated multiRun() output, and θ = parameter matrix.
simulateMr <- function(dims, N)
{
  library(morpheus)
  mr <- list()
  θ <- matrix(runif(min=-5,max=5,n=dims[1]*dims[2]), ncol=dims[1])
  for (i in 1:2) {
    mr[[i]] <- list()
    for (j in 1:N)
      mr[[i]][[j]] <- θ + matrix(rnorm(dims[1]*dims[2],sd=0.25),ncol=dims[1])
    mr[[i]] <- alignMatrices(mr[[i]], ref=θ, ls_mode="exact")
  }
  list(mr = mr, θ = θ)
}

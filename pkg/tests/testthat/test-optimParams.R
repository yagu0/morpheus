naive_f <- function(link, M1,M2,M3, p,β,b)
{
  d <- length(M1)
  K <- length(p)
  λ <- sqrt(colSums(β^2))

  # Compute β x2,3 (self) tensorial products
  β2 <- array(0, dim=c(d,d,K))
  β3 <- array(0, dim=c(d,d,d,K))
  for (k in 1:K)
  {
    for (i in 1:d)
    {
      for (j in 1:d)
      {
        β2[i,j,k] = β[i,k]*β[j,k]
        for (l in 1:d)
          β3[i,j,l,k] = β[i,k]*β[j,k]*β[l,k]
      }
    }
  }

  res <- 0
  for (i in 1:d)
  {
    term <- 0
    for (k in 1:K)
      term <- term + p[k]*.G(link,1,λ[k],b[k])*β[i,k]
    res <- res + (term - M1[i])^2
    for (j in 1:d)
    {
      term <- 0
      for (k in 1:K)
        term <- term + p[k]*.G(link,2,λ[k],b[k])*β2[i,j,k]
      res <- res + (term - M2[i,j])^2
      for (l in 1:d)
      {
        term <- 0
        for (k in 1:K)
          term <- term + p[k]*.G(link,3,λ[k],b[k])*β3[i,j,l,k]
        res <- res + (term - M3[i,j,l])^2
      }
    }
  }
  res
}

# TODO: understand why delta is so large (should be 10^-6 10^-7 ...)
test_that("naive computation provides the same result as vectorized computations",
{
  set.seed(7)
  h <- 1e-7 #for finite-difference tests
  n <- 10
  for (dK in list( c(2,2), c(5,3)))
  {
    d <- dK[1]
    K <- dK[2]

    M1 <- runif(d, -1, 1)
    M2 <- matrix(runif(d^2, -1, 1), ncol=d)
    M3 <- array(runif(d^3, -1, 1), dim=c(d,d,d))

    for (link in c("logit","probit"))
    {
      # X and Y are unused here (W not re-computed)
      op <- optimParams(X=matrix(runif(n*d),ncol=d), Y=rbinom(n,1,.5),
        K, link, M=list(M1,M2,M3))
      op$W <- diag(d + d^2 + d^3)

      for (var in seq_len((2+d)*K-1))
      {
        p <- runif(K, 0, 1)
        p <- p / sum(p)
        β <- matrix(runif(d*K,-5,5),ncol=K)
        b <- runif(K, -5, 5)
        x <- c(p[1:(K-1)],as.double(β),b)

        # Test functions values (TODO: 1 is way too high)
        expect_equal( op$f(x)[1], naive_f(link,M1,M2,M3, p,β,b), tolerance=1 )

        # Test finite differences ~= gradient values
        dir_h <- rep(0, (2+d)*K-1)
        dir_h[var] = h
        expect_equal( op$grad_f(x)[var], ((op$f(x+dir_h) - op$f(x)) / h)[1], tolerance=0.5 )
      }
    }
  }
})

test_that("W computed in C and in R are the same",
{
  set.seed(7)
  tol <- 1e-8
  n <- 10
  for (dK in list( c(2,2))) #, c(5,3)))
  {
    d <- dK[1]
    K <- dK[2]
    link <- ifelse(d==2, "logit", "probit")
    θ <- list(
      p=rep(1/K,K),
      β=matrix(runif(d*K),ncol=K),
      b=rep(0,K))
    io <- generateSampleIO(n, θ$p, θ$β, θ$b, link)
    X <- io$X
    Y <- io$Y
    dd <- d + d^2 + d^3
    p <- θ$p
    β <- θ$β
    λ <- sqrt(colSums(β^2))
    b <- θ$b
    β2 <- apply(β, 2, function(col) col %o% col)
    β3 <- apply(β, 2, function(col) col %o% col %o% col)
    M <- c(
      β  %*% (p * .G(link,1,λ,b)),
      β2 %*% (p * .G(link,2,λ,b)),
      β3 %*% (p * .G(link,3,λ,b)))
    Id <- as.double(diag(d))
    E <- diag(d)
    v1 <- Y * X
    v2 <- Y * t( apply(X, 1, function(Xi) Xi %o% Xi - Id) )
    v3 <- Y * t( apply(X, 1, function(Xi) { return (Xi %o% Xi %o% Xi
      - Reduce('+', lapply(1:d, function(j)
        as.double(Xi %o% E[j,] %o% E[j,])), rep(0, d*d*d))
      - Reduce('+', lapply(1:d, function(j)
        as.double(E[j,] %o% Xi %o% E[j,])), rep(0, d*d*d))
      - Reduce('+', lapply(1:d, function(j)
        as.double(E[j,] %o% E[j,] %o% Xi)), rep(0, d*d*d))) } ) )
    Omega1 <- matrix(0, nrow=dd, ncol=dd)
    for (i in 1:n)
    {
      gi <- t(as.matrix(c(v1[i,], v2[i,], v3[i,]) - M))
      Omega1 <- Omega1 + t(gi) %*% gi / n
    }
    W <- matrix(0, nrow=dd, ncol=dd)
    Omega2 <- matrix( .C("Compute_Omega",
      X=as.double(X), Y=as.integer(Y), M=as.double(M),
      pnc=as.integer(1), pn=as.integer(n), pd=as.integer(d),
      W=as.double(W), PACKAGE="morpheus")$W, nrow=dd, ncol=dd )
    rg <- range(Omega1 - Omega2)
    expect_equal(rg[1], rg[2], tolerance=tol)
  }
})

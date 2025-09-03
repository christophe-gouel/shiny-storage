library(dplyr)
library(ggplot2)
library(markdown)
library(moments)
library(statmod)
library(tidyr)

SolveStorage <- function(model) {
  # Initialization
  pbar <- model$params$pbar
  dbar <- model$params$dbar
  k <- model$params$k
  delta <- model$params$delta
  r <- model$params$r
  elastD <- model$params$elastD
  SDe <- model$params$SDe

  e <- model$shocks$e
  w <- model$shocks$w

  beta <- (1 - delta) / (1 + r)
  P <- model$P

  demand <- function(p) dbar * (1 + elastD * (p - pbar) / pbar)
  invdemand <- function(d) pbar * (1 + (d - dbar) / (elastD * dbar))

  A <- model$s

  PA <- invdemand(A) # Precalculation for speed

  dis <- Inf
  Display <- FALSE
  Iter <- 0
  MaxIter <- 1E3
  TolX <- 1E-5

  n <- length(A)
  kshocks <- length(e)

  PriceInterp <- splinefun(A, P, method = "monoH.FC")

  while (dis > TolX & Iter < MaxIter) {
    Iter <- Iter + 1
    Pold <- P

    # Calculate next-period availability
    S <- A - demand(P)
    Anext <-
      apply((1 - delta) * S, MARGIN = 1, function(x) {
        x + dbar * (1 + SDe * e)
      }) |>
      t() |>
      as.vector()

    # Update the price and its approximation
    Pnext <- PriceInterp(Anext)
    P <- pmax(
      PA,
      beta * matrix(Pnext, nrow = n, ncol = kshocks) %*% w - k * pbar
    )
    PriceInterp <- splinefun(A, P, method = "monoH.FC")

    dis <- max(abs(P - Pold))
    if (Display) print(c(Iter, dis))
  }

  model$SolveStat <- list(dis = dis, Iter = Iter, exitflag = Iter < MaxIter)
  model$P <- P
  model$S <- S
  model$PriceInterp <- PriceInterp
  model$StorageInterp <- splinefun(A, S, method = "monoH.FC")

  return(model)
}

SimulateStorage <- function(model, A0, nper = 100, nrep = 100, nburn = 20) {
  StorageInterp <- model$StorageInterp
  PriceInterp <- model$PriceInterp
  dbar <- model$params$dbar
  delta <- model$params$delta
  SDe <- model$params$SDe
  ntot <- nper + nburn
  A <- matrix(data = A0, nrow = nrep, ncol = ntot)
  S <- matrix(nrow = nrep, ncol = ntot)
  P <- matrix(nrow = nrep, ncol = ntot)
  set.seed(1)
  e <- matrix(
    data = 1 + SDe * c(rep(NA, nrep), rnorm(nrep * (ntot - 1))),
    nrow = nrep,
    ncol = ntot
  )
  for (t in 1:ntot) {
    if (t > 1) {
      A[, t] <- (1 - delta) * S[, t - 1] + dbar * e[, t]
    }
    S[, t] <- StorageInterp(A[, t])
    P[, t] <- PriceInterp(A[, t])
  }
  return(list(
    A = A[, (nburn + 1):ntot],
    S = S[, (nburn + 1):ntot],
    P = P[, (nburn + 1):ntot],
    e = e[, (nburn + 1):ntot]
  ))
}

gherm <- gauss.quad(9, kind = "hermite")

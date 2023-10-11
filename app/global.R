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
  δ <- model$params$δ
  r <- model$params$r
  elastD <- model$params$elastD
  SDe <- model$params$SDe

  e <- model$shocks$e
  w <- model$shocks$w

  β <- (1 - δ)/(1 + r)
  P <- model$P

  demand <- function(p) dbar * (p / pbar)^elastD
  invdemand <- function(d) pbar * (d / dbar)^(1 / elastD)

  S <- model$A

  dis <- Inf
  Display <- FALSE
  Iter <- 0
  MaxIter <- 1E3
  TolX <- 1E-5

  n <- length(A)
  kshocks <- length(e)

  PriceInterp <- splinefun(A, P, method = "monoH.FC")

  while(dis > TolX & Iter < MaxIter) {
    Iter <- Iter + 1
    Pold <- P

    # Calculate next-period availability
    Anext <-
      apply((1 - δ) * S, MARGIN = 1,
            function(x) x + dbar * exp(SDe * e)) %>%
      t() %>%
      as.vector()

    # Update the price and its approximation
    Pnext <- invdemand(DemandFunction(Anext))
    P <- β * matrix(Pnext, nrow = n, ncol = kshocks) %*% w - k * pbar
    D <- invdemand(P)
    A <- S + D
    DemandInterp <- splinefun(A, D, method = "monoH.FC")
    DemandFunction <- function(a) {
      D <- a
      D[a > A[1]] <- DemandInterp(a[a > A[1]])
    }

    dis <- max(abs(P - Pold))
    if(Display) print(c(Iter, dis))
  }

  model$SolveStat <- list(dis = dis, Iter = Iter, exitflag = Iter < MaxIter)
  model$A <- A
  model$P <- P
  model$PriceInterp <- PriceInterp
  model$StorageInterp <- splinefun(A, S, method = "monoH.FC")

  return(model)
}

SimulateStorage <- function(model, A0, nper = 100, nrep = 100, nburn = 20) {
  StorageInterp <- model$StorageInterp
  PriceInterp <- model$PriceInterp
  dbar  <- model$params$dbar
  δ <- model$params$δ
  SDe <- model$params$SDe
  ntot <- nper + nburn
  A <- matrix(data = A0, nrow = nrep, ncol = ntot)
  S <- matrix(nrow = nrep, ncol = ntot)
  P <- matrix(nrow = nrep, ncol = ntot)
  set.seed(1)
  e <- c(rep(NA, nrep), rnorm(nrep * (ntot - 1))) %>%
    matrix(nrow = nrep, ncol = ntot)
  for (t in 1:ntot) {
    if (t > 1) A[,t] <- (1 - δ) * S[,t-1] + dbar * (1 + SDe * e[,t])
    S[,t] <- StorageInterp(A[,t])
    P[,t] <- PriceInterp(A[,t])
  }
  return(list(A = A[,(nburn+1):ntot],
              S = S[,(nburn+1):ntot],
              P = P[,(nburn+1):ntot],
              e = e[,(nburn+1):ntot]))
}

gherm <- gauss.quad(9, kind = "hermite")


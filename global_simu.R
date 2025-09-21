rm(list = ls())
library(BPST)
library(TPST)
library(dplyr)
library(Triangulation)
library(pracma)
library(MASS)
library(mgcv)
library(splines2)
library(Matrix)
require(MGLM)
source('funcs/simul2.data.gen.R')
source('funcs/basis.tensor.R')
source('funcs/energy.tensor.R')
source('funcs/gc.fit.R')
source('funcs/gc.fit.gcv.R')

simulation <- function(iter){
  set.seed(iter)
  # Set boundary 
  xm <- seq(-1, 3.5, length = 200)
  yn <- seq(-1, 1, length = 200)
  xy_grid <- pracma::meshgrid(xm, yn)
  uu <- c(xy_grid$X)
  vv <- c(xy_grid$Y)
  V <- as.matrix(cbind(uu,vv))
  
  # Load horseshoe and find triangulation 
  bb <- Triangulation::hs
  VT <- TriMesh(bb, n = nT)
  Tr <- as.matrix(VT$Tr) 
  Ver <- as.matrix(VT$V) 
  
  # Generate simulation data
  t0 <- proc.time()
  dat <- simul.data.gen(n, Tr, Ver, V)
  time.bound<- c(min(dat$tij), max(dat$tij))
  probs <- seq(0, 1, length.out = L + 2)
  time.knots <- quantile(dat$tij, probs = probs)[-c(1, L + 2)]
  Basis <- basis.tensor(ss = dat$V, tt = dat$tij, V = Ver, Tri = Tr, 
                        d = d, r = r, time.knots = time.knots, time.bound = time.bound, rho = rho)
  X <- dat$X
  Y <- dat$Y
  V <- dat$V
  M <- dat$M.vec
  tij <- dat$tij
  U0 <- Basis$U0
  B0 <- Basis$B0
  Q2 <- Basis$Q2
  K <- Basis$K
  P11 <- Basis$D1
  P12 <- Basis$D2
  P2 <- kronecker(diag(ncol(X)), as.matrix(crossprod(Q2, K) %*% Q2))
  
  # find optimal tuning parameters 
  gcv.fit <- gc.fit.gcv(X, Y, V, M, U0, B0, Q2, K, P11, P12, P2, 
                        lambda11 = exp(seq(log(0.001), log(1000), length.out = 5)),
                        lambda12 = exp(seq(log(0.001), log(1000), length.out = 5)),
                        lambda2 = exp(seq(log(0.001), log(1000), length.out = 5)))
  lambda11.optimal <- gcv.fit$lambda11
  lambda12.optimal <- gcv.fit$lambda12
  lambda2.optimal <- gcv.fit$lambda2
  
  # fit the model with optimal lambdas
  res <- gc.fit(X, Y, V, M, U0, B0, Q2, K, P11, P12, P2, L, rho,
                 lambda11=lambda11.optimal, lambda12=lambda12.optimal, lambda2=lambda2.optimal)
  t1 <- proc.time() - t0
  bivariate.true <- cbind(dat$bivar.alpha, dat$bivar.beta)
  bivar.est <- res$bivar.est
  trivariate.true <- dat$tri.alpha
  trivariate.est <- res$trivar.est
  trivar.all <- c()
  for(i in 1:nrow(trivariate.est)){
    trivar.all[i] <- mean((trivariate.est[i, ] - trivariate.true[i, ])^2, na.rm = TRUE)
  }
  final.res <- c(round(c(mean(trivar.all), apply((bivariate.true - bivar.est)^2, 2, mean, na.rm = TRUE)), 4), round(t1['elapsed']))
  names(final.res)[1] <- 'alpha0'
  return(final.res)
}

# Set parameters
L <- 3
rho <- 3
nT <- 6
n <- 100; d <- 2; r <- 1
result <- lapply(c(1:100), simulation)
round(colMeans(do.call(rbind, result)), 4)

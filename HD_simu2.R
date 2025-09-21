rm(list = ls())
library(devtools)
library(BPST)
library(TPST)
library(dplyr)
library(Triangulation)
library(pracma)
library(MASS)
library(mgcv)
library(splines2)
require(MGLM)
library(Matrix)
library(gghilbertstrings)
library(parallel)
library(prodlim)

source('funcs/simul2.data.gen.R')
source('funcs/basis.tensor.R')
source('funcs/energy.tensor.R')
source('funcs/gc.fit.R')
source('funcs/gc.fit.gcv.R')
source('funcs/gc.fit.dc.R')
source('funcs/sampling.HC.R')
source('funcs/ring.dc.R')
source('funcs/basis.tensor.local.R')
source('funcs/local.fit.R')

ddc.fit.simul <- function(iter){
  set.seed(iter)
  bb <- Triangulation::hs
  VT <- TriMesh(bb, n = nT)
  Tr <- as.matrix(VT$Tr) 
  Ver <- as.matrix(VT$V) 
  
  # Set boundary 
  xm <- seq(-1, 3.5, length = 200)
  yn <- seq(-1, 1, length = 200)
  xy_grid <- pracma::meshgrid(xm, yn)
  uu <- c(xy_grid$X)
  vv <- c(xy_grid$Y)
  V <- as.matrix(cbind(uu,vv))
  
  # Generate dataset
  t0 <- proc.time()
  dat <- simul.data.gen(n, Tr, Ver, V)
  L <- 3
  rho <- 3
  time.bound<- c(min(dat$tij), max(dat$tij))
  probs <- seq(0, 1, length.out = L + 2)
  time.knots <- quantile(dat$tij, probs = probs)[-c(1, L + 2)]
  
  sampling.res <- sampling.HC(n.samp, Ver, Tr, n.layer)
  count.tri <- rep(0, nrow(Tr))
  TV <- tdata(Ver, Tr)$TV
  idx.sample.tri <- sampling.res$sample.tri
  
  # fits estimators using HD method
  lambda11.r = exp(seq(log(0.001), log(1000), length.out = 5))
  lambda12.r = exp(seq(log(0.001), log(1000), length.out = 5))
  lambda2.r = exp(seq(log(0.001), log(1000), length.out = 5))
  fit.all <- mclapply(idx.sample.tri, function(iter) {
    local.fit(iter, Ver0 = Ver, Tr0 = Tr, TV0 = TV, n.layer = n.layer, 
              X = dat$X, Y = dat$Y, V = dat$V, M = dat$M.vec, tij = dat$tij,
              d = d, r = r, L = L, rho = rho, 
              time.knots = time.knots, time.bound = time.bound,
              lambda11 = exp(seq(log(0.001), log(1000), length.out = 5)), 
              lambda12 = exp(seq(log(0.001), log(1000), length.out = 5)), 
              lambda2 = exp(seq(log(0.001), log(1000), length.out = 5)))
  }, mc.cores = n.core, mc.preschedule = TRUE)
  
  gamma.all <- matrix(0, ncol = L+rho, nrow = nrow(Tr)*(d+1)*(d+2)/2)
  theta.all <- matrix(0, ncol = ncol(dat$X), nrow = nrow(Tr)*(d+1)*(d+2)/2)
  
  for(iter in 1:length(fit.all)) {
    idx.tr <- fit.all[[iter]]$idx.tr
    count.tri[idx.tr] <- count.tri[idx.tr] + 1
    gamma.all[fit.all[[iter]]$idx.psi,] = gamma.all[fit.all[[iter]]$idx.psi,] + fit.all[[iter]]$gamma.hat
    theta.all[fit.all[[iter]]$idx.psi,] = theta.all[fit.all[[iter]]$idx.psi,] + fit.all[[iter]]$theta.hat
  }
  
  nbasis.tri <- (d+1)*(d+2)/2
  count.psi <- rep(count.tri, each = nbasis.tri)
  mean.gamma <- 1/count.psi * gamma.all
  mean.theta <- 1/count.psi * theta.all
  
  p = dim(dat$X)[2]
  Basis <- basis.tensor(ss = dat$V, tt = dat$tij, V = Ver, Tri = Tr,
                        d = d, r = r, time.knots = time.knots, time.bound = time.bound, rho = rho)
  
  est.gamma = solve(crossprod(Basis$Q2), crossprod(Basis$Q2, mean.gamma))
  est.theta = solve(crossprod(Basis$Q2), crossprod(Basis$Q2, mean.theta))
  
  alpha0.hat <- tcrossprod(Basis$U0, (Basis$B0 %*% Basis$Q2 %*% est.gamma)) 
  bivar.est <- Basis$B0 %*% Basis$Q2 %*% est.theta 
  
  t1 <- proc.time() - t0
  bivariate.true <- cbind(dat$bivar.alpha, dat$bivar.beta)
  trivariate.true <- dat$tri.alpha
  trivariate.est <- alpha0.hat
  trivar.all <- c()
  for(i in 1:nrow(trivariate.est)){
    trivar.all[i] <- mean((trivariate.est[i, ] - trivariate.true[i, ])^2, na.rm = TRUE)
  }
  
  final.res <- c(round(c(mean(trivar.all), apply((bivariate.true - bivar.est)^2, 2, mean, na.rm = TRUE)), 4), round(t1['elapsed']))
  names(final.res)[1] <- 'alpha0'
  return(final.res)
}

# Set parameters
n <- 100
d <- 2; r <- 1
nT <- 6
n.samp <- 20; n.layer <- 3; n.core <- 8                  
res <- lapply(c(1:100), ddc.fit.simul)
round(colMeans(do.call(rbind, res)), 4)

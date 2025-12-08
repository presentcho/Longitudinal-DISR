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

simulation <- function(iter){
  set.seed(iter)
  # Generate simulation data
  dat <- simul.data.gen(n, Tr, Ver, V.all)
  X <- dat$X
  Y <- dat$Y
  V <- dat$V
  M <- dat$M.vec
  tij <- dat$tij
  bivar.true <- cbind(dat$bivar.alpha, dat$bivar.beta)
  trivar.true <- dat$tri.alpha

  cv.fit <- gc.fit.cv.k(nfold=5,X=X, Y=Y, V=V, M=M, tij=tij)
  
  lambda11.optimal <- cv.fit$lambda11
  lambda12.optimal <- cv.fit$lambda12
  lambda2.optimal <- cv.fit$lambda2
  
  sampling.res <- sampling.HC(n.samp, Ver, Tr, n.layer)
  count.tri <- rep(0, nrow(Tr))
  TV <- tdata(Ver, Tr)$TV
  idx.sample.tri <- sampling.res$sample.tri
  
  fit.dc <- gc.fit.dc.k(idx.sample.tri = idx.sample.tri, Tr = Tr, Ver = Ver, TV = TV, 
                                  n.layer = n.layer, X = X, Y = Y, V = V, M= M, tij = tij,
                                  d=d, r=r, L=L, rho=rho, 
                      lambda11=lambda11.optimal, lambda12 = lambda12.optimal, lambda2= lambda2.optimal)
  bivar.est <- as.matrix(fit.dc$bivar.est)
  trivar.est <- as.matrix(fit.dc$trivar.est)
  
  boot_trivar <- list()
  boot_beta <- array(NA_real_, dim = c(nrow(bivar.est), ncol(bivar.est), nboot))
  for(b in 1:nboot){
    sampling <- sample(nrow(X), replace = TRUE)
    X.s <- X[sampling,]
    Y.s <- Y[sampling,]
    M.s <- M[sampling]
    t.s <- tij[sampling]
    
    fit.b <- gc.fit.dc.k(idx.sample.tri = idx.sample.tri, Tr = Tr, Ver = Ver, TV = TV, 
                          n.layer = n.layer, X = X.s, Y = Y.s, V = V, M= M.s, tij = t.s,
                          d=d, r=r, L=L, rho=rho,
                          lambda11=lambda11.optimal, lambda12 = lambda12.optimal, lambda2= lambda2.optimal)
    boot_beta[,,b] <- as.matrix(fit.b$bivar.est)
    boot_trivar[[b]] <- as.matrix(fit.b$trivar.est)
  }
  # obtain bootstrap variance matrix
  bivar.se.mat <- apply(boot_beta, c(1,2), sd, na.rm = TRUE)
  trivar.se.mat <- apply(simplify2array(boot_trivar), c(1,2), sd)
  bivar.var.mat <- apply(boot_beta, c(1,2), var, na.rm = TRUE)
  trivar.var.mat <- apply(simplify2array(boot_trivar), c(1,2), var)
  var.res <- round(c(mean(trivar.var.mat), colMeans(bivar.var.mat)), 4)
  
  # obtain coverage probability
  # 95% coverage rate
  trivar.coverage.lower.95 <- trivar.est + qnorm(0.05/2)*trivar.se.mat
  trivar.coverage.upper.95 <- trivar.est + qnorm(1-0.05/2)*trivar.se.mat
  bivar.coverage.lower.95 <- bivar.est + qnorm(0.05/2)*bivar.se.mat
  bivar.coverage.upper.95 <- bivar.est + qnorm(1-0.05/2)*bivar.se.mat
  trivar.covered.95 <- (trivar.true >= trivar.coverage.lower.95) & (trivar.true <= trivar.coverage.upper.95)
  bivar.covered.95 <- (bivar.true >= bivar.coverage.lower.95) & (bivar.true <= bivar.coverage.upper.95)
  covered.bivar.95 <- colMeans(bivar.covered.95)
  covered.trivar.95 <- mean(rowMeans(trivar.covered.95))
  covered.95 <- c(covered.trivar.95, covered.bivar.95)
  
  # 90% coverage rate
  trivar.coverage.lower.90 <- trivar.est + qnorm(0.1/2)*trivar.se.mat
  trivar.coverage.upper.90 <- trivar.est + qnorm(1-0.1/2)*trivar.se.mat
  bivar.coverage.lower.90 <- bivar.est + qnorm(0.1/2)*bivar.se.mat
  bivar.coverage.upper.90 <- bivar.est + qnorm(1-0.1/2)*bivar.se.mat
  trivar.covered.90 <- (trivar.true >= trivar.coverage.lower.90) & (trivar.true <= trivar.coverage.upper.90)
  bivar.covered.90 <- (bivar.true >= bivar.coverage.lower.90) & (bivar.true <= bivar.coverage.upper.90)
  covered.bivar.90 <- colMeans(bivar.covered.90)
  covered.trivar.90 <- mean(rowMeans(trivar.covered.90))
  covered.90 <- c(covered.trivar.90, covered.bivar.90)
  
  # 99% coverage rate
  trivar.coverage.lower.99 <- trivar.est + qnorm(0.01/2)*trivar.se.mat
  trivar.coverage.upper.99 <- trivar.est + qnorm(1-0.01/2)*trivar.se.mat
  bivar.coverage.lower.99 <- bivar.est + qnorm(0.01/2)*bivar.se.mat
  bivar.coverage.upper.99 <- bivar.est + qnorm(1-0.01/2)*bivar.se.mat
  trivar.covered.99 <- (trivar.true >= trivar.coverage.lower.99) & (trivar.true <= trivar.coverage.upper.99)
  bivar.covered.99 <- (bivar.true >= bivar.coverage.lower.99) & (bivar.true <= bivar.coverage.upper.99)
  covered.bivar.99 <- colMeans(bivar.covered.99)
  covered.trivar.99 <- mean(rowMeans(trivar.covered.99))
  covered.99 <- c(covered.trivar.99, covered.bivar.99)
  
  # bias and MISE
  trivar.all <- c()
  trivar.bias <- c()
  for(i in 1:nrow(trivar.est)){
    trivar.all[i] <- mean((trivar.est[i, ] - trivar.true[i, ])^2, na.rm = TRUE)
    trivar.bias[i] <- mean((trivar.est[i, ] - trivar.true[i, ]), na.rm = TRUE)
  }
  
  MISE.res <- c(round(c(mean(trivar.all), apply((bivar.true - bivar.est)^2, 2, mean, na.rm = TRUE)), 4))
  names(MISE.res)[1] <- 'alpha0'
  bias.res <- c(round(c(mean(trivar.bias), apply((bivar.true - bivar.est), 2, mean, na.rm = TRUE)), 4))
  
  
  return(list(bias = bias.res, MISE = MISE.res,
              bivar.covered.90=bivar.covered.90,
              trivar.covered.90=trivar.covered.90,
              bivar.covered.95=bivar.covered.95,
              trivar.covered.95=trivar.covered.95,
              bivar.covered.99=bivar.covered.99,
              trivar.covered.99=trivar.covered.99,
              variance=var.res ))
}

L <- 3
rho <- 3
nT <- 6
n.samp <- 20; n.layer <- 2; n.core <- 2
n <- 100; d <- 2; r <- 1
nboot<- 100

# Set boundary 
xm <- seq(-1, 3.5, length = 101)
yn <- seq(-1, 1, length = 101)
xy_grid <- pracma::meshgrid(xm, yn)
uu <- c(xy_grid$X)
vv <- c(xy_grid$Y)
V.all <- as.matrix(cbind(uu,vv))

# Load horseshoe and find triangulation 
bb <- Triangulation::hs
VT <- TriMesh(bb, n = nT)
Tr <- as.matrix(VT$Tr) 
Ver <- as.matrix(VT$V) 

t0 <- proc.time()
res <- lapply(1:100, function(iter) simulation(iter))
proc.time() - t0
bias_mat <- do.call(rbind, lapply(res, `[[`, "bias"))
MISE_mat <- do.call(rbind, lapply(res, `[[`, "MISE"))
var_mat  <- do.call(rbind, lapply(res, `[[`, "variance"))

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
source('funcs/gc.fit.cv.R')
source('funcs/gc.fit.ho.R')

simulation <- function(iter){
  set.seed(iter)
  # Generate simulation data
  dat <- simul.data.gen(n, Tr, Ver, V.all)
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
  cv.fit <- gc.fit.cv(nfold = 5, X, Y, V, M, tij, Basis)
  lambda11.optimal <- cv.fit$lambda11
  lambda12.optimal <- cv.fit$lambda12
  lambda2.optimal <- cv.fit$lambda2
  res <- gc.fit(X, Y, V, M, U0, B0, Q2, K, P11, P12, P2, L, rho,
                 lambda11=lambda11.optimal, lambda12=lambda12.optimal, lambda2=lambda2.optimal)
  bivar.true <- cbind(dat$bivar.alpha, dat$bivar.beta)
  bivar.est <- res$bivar.est
  trivar.true <- dat$tri.alpha
  trivar.est <- res$trivar.est
  
  boot_trivar <- list()
  boot_beta <- array(NA_real_, dim = c(nrow(bivar.est), ncol(bivar.est), nboot))
  for(b in 1:nboot){
    sampling <- sample(nrow(X), replace = TRUE)
    X.s <- X[sampling,]
    Y.s <- Y[sampling,]
    M.s <- M[sampling]
    U0.s <- U0[sampling,]
    b.res <- gc.fit(X.s, Y.s, V.s, M.s, U0.s, B0, Q2, K, P11, P12, P2, L, rho,
                   lambda11=lambda11.optimal, lambda12=lambda12.optimal, lambda2=lambda2.optimal)
    boot_beta[,,b] <- b.res$bivar.est
    boot_trivar[[b]] <- b.res$trivar.est
  }
  # obtain bootstrap variance matrix
  bivar.se.mat <- apply(boot_beta, c(1,2), sd)
  trivar.se.mat <- apply(simplify2array(boot_trivar), c(1,2), sd)
  bivar.var.mat <- apply(boot_beta, c(1,2), var)
  trivar.var.mat <- apply(simplify2array(boot_trivar), c(1,2), var)
  var.res <- round(c(mean(trivar.var.mat), colMeans(bivar.var.mat)), 4)
  
  # obtain coverage probability
  # 95% coverage rate
  bivar.coverage.lower.95 <- bivar.est + qnorm(0.05/2)*bivar.se.mat
  bivar.coverage.upper.95 <- bivar.est + qnorm(1-0.05/2)*bivar.se.mat
  bivar.covered.95 <- (bivar.true >= bivar.coverage.lower.95) & (bivar.true <= bivar.coverage.upper.95)
  covered.bivar.95 <- colMeans(bivar.covered.95)
  
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
              bivar.covered.95=bivar.covered.95,
              variance=var.res))
}


L <- 3
rho <- 3
nT <- 6
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

res <- lapply(1:100, function(iter) simulation(iter))

bias_mat <- do.call(rbind, lapply(res, `[[`, "bias"))
MISE_mat <- do.call(rbind, lapply(res, `[[`, "MISE"))
var_mat  <- do.call(rbind, lapply(res, `[[`, "variance"))
bivar_cr_mat_95  <- do.call(rbind, lapply(res, `[[`, "bivar.covered.95"))
round(colMeans(bias_mat),4)
round(colMeans(MISE_mat),4)
round(colMeans(var_mat),4)
round(colMeans(bivar_cr_mat_95), 4)

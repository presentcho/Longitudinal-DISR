rm(list = ls())
library(BPST)
library(TPST)
library(dplyr)
library(Triangulation)
library(pracma)
library(MASS)
library(mgcv)
library(colorRamps)

library(splines2)
require(MGLM)
source('funcs/simul1.data.gen.R')
source('funcs/bivar.gc.fit.R')

################################################################################
simulation <- function(iter){
  set.seed(iter)
  # Generate simulation data
  dat <- simul.data.gen(n, Tr, Ver, V.all)
  Y <- dat$Y 
  X <- dat$X 
  V <- dat$V
  M <- dat$M
  ind.inside <- dat$ind.inside
  alpha.true <- dat$bivar.alpha
  beta.true <- dat$bivar.beta
  bivar.true <- cbind(alpha.true, beta.true)
  
  # bootstrap sampling
  lambda <- 10^(seq(-6, 0, by = 1))
  Y.est <- Y[,ind.inside]
  V <- V.all[ind.inside, ]
  
  # get basis
  if(d > 1){
    Ball <- basis(Ver, Tr, d, r, V)
    K <- Ball$K
    Q2 <- Ball$Q2
    B <- Ball$B 
  }
  
  X_tilde <- M * X
  Y_tilde <- M * Y.est
  mfit0 <- bivar.gc.fit(B, Q2, K, lambda, X_tilde, Y_tilde)
  bivar.est <- mfit0$beta
  
  boot_beta <- array(NA_real_, dim = c(nrow(bivar.est), ncol(bivar.est), nboot))
  for(b in 1:nboot){
    sampling <- sample(nrow(X), replace = TRUE)
    X.s <- X[sampling,]
    Y.s <- Y[sampling,]
    M.s <- M[sampling]
    Y.est.s <- Y.s[,ind.inside]
    X_tilde.s <- M.s * X.s
    Y_tilde.s <- M.s * Y.est.s
    boot_beta[,,b] <- bivar.gc.fit(B, Q2, K, lambda, X_tilde.s, Y_tilde.s)$beta
  }
  # obtain bootstrap variance matrix
  variance.mat <- apply(boot_beta, c(1,2), var, na.rm = TRUE)
  se.mat <- apply(boot_beta, c(1,2), sd, na.rm = TRUE)
  # obtain coverage probability
  # 95% coverage rate
  coverage.prob.lower_95 <- bivar.est + qnorm(0.05/2)*se.mat
  coverage.prob.upper_95 <- bivar.est + qnorm(1-0.05/2)*se.mat
  covered_95 <- (bivar.true[ind.inside,] >= coverage.prob.lower_95) & (bivar.true[ind.inside,] <= coverage.prob.upper_95)
  covered.mean_95 <- colMeans(covered_95)
  
  # obtain bias
  bias.res <- apply((bivar.est - bivar.true[ind.inside,]), 2, mean, na.rm = TRUE)
  # obtain mise
  mse.res <- apply((bivar.est - bivar.true[ind.inside,])^2, 2, mean, na.rm = TRUE)
  return(list(bias = bias.res, MISE = mse.res,
              coverage.rate.95=covered.mean_95, 
              covered_95_mat = covered_95,
              variance = colMeans(variance.mat)))
}

nT <- 6
n <- 100; d <- 2; r <- 1; alpha <- 0.05
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
cr_mat_95  <- do.call(rbind, lapply(res, `[[`, "coverage.rate.95"))
var_mat  <- do.call(rbind, lapply(res, `[[`, "variance"))
cr_95_map <- Reduce("+", lapply(res, `[[`, 'covered_95_mat')) / length(res)

round(colMeans(bias_mat),4)
round(colMeans(MISE_mat),4)
round(colMeans(var_mat),4)
round(colMeans(cr_95_map), 4)

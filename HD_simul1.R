rm(list = ls())
library(BPST)
library(TPST)
library(dplyr)
library(Triangulation)
library(pracma)
library(MASS)
library(mgcv)
library(colorRamps)
library(gghilbertstrings)
library(parallel)
library(splines2)
require(MGLM)
source('funcs/simul1.data.gen.R')
source('funcs/bivar.gc.fit.R')
source('funcs/ring.dc.R')
source('funcs/sampling.HC.R')
source('funcs/local.fit.bivar.R')

################################################################################
simulation <- function(iter){
  set.seed(iter)
  # Generate datasets and true coefficient parameters
  dat <- simul.data.gen(n, Tr, Ver, V.all)
  Y <- dat$Y 
  X <- dat$X 
  V <- dat$V
  M <- dat$M
  ind.inside <- dat$ind.inside
  alpha.true <- dat$bivar.alpha
  beta.true <- dat$bivar.beta
  bivar.true <- cbind(alpha.true, beta.true)
  
  lambda <- 10^(seq(-6, 0, by = 1))
  Y.est <- Y[,ind.inside]
  V <- V.all[ind.inside, ]
  
  # Construct bivariate spline basis functions over the triangulation
  if(d > 1){
    Ball <- basis(Ver, Tr, d, r, V)
    K <- Ball$K
    Q2 <- Ball$Q2
    B <- Ball$B 
  }
  
  X_tilde <- M * X
  Y_tilde <- M * Y.est
  result <- sampling.HC(n.samp, Ver, Tr, n.layer = n.layer)
  gamma.all <- matrix(0, ncol = ncol(X_tilde), nrow = nrow(Tr) * (d+1) * (d+2)/2)
  count.tri <- rep(0, nrow(Tr))
  nbasis.tri <- (d+1) * (d+2)/2
  TV <- tdata(Ver, Tr)$TV
  idx.sample.tri <- result$sample.tri
  
  # Execute HD method- local fitting across subregions and aggregate
  fit.all <- mclapply(idx.sample.tri, FUN = local.fit.bivar, mc.cores = core,
                      V0 = Ver, Tr0 = Tr, TV0 = TV, n.layer = n.layer, 
                      Y.all = Y_tilde, X.all = X_tilde,
                      Z.all = V, d = d, r = r, lambda = 10^(seq(-6, 0, by = 1)))
  gamma.all <- matrix(0, ncol = ncol(X_tilde), nrow = nrow(Tr) * (d+1) * (d+2)/2)
  count.tri <- rep(0, nrow(Tr))
  
  for(iter in 1:length(fit.all)) {
    idx.tr <- fit.all[[iter]]$idx.tr
    count.tri[idx.tr] <- count.tri[idx.tr] + 1
    gamma.all[fit.all[[iter]]$idx.gamma, ] = 
      gamma.all[fit.all[[iter]]$idx.gamma, ] + fit.all[[iter]]$gamma.local
  }
  
  count.gamma <- rep(count.tri, each = nbasis.tri)
  mean.gamma <- 1/count.gamma * gamma.all
  
  p = dim(X)[2]
  Ball = basis(Ver,Tr,d,r,V)
  theta.hat <- solve(crossprod(Ball$Q2), crossprod(Ball$Q2, mean.gamma))
  bivar.est = Ball$B %*% Ball$Q2 %*% theta.hat
  
  boot_beta <- array(NA_real_, dim = c(nrow(bivar.est), ncol(bivar.est), nboot))
  for(b in 1:nboot){
    sampling <- sample(nrow(X), replace = TRUE)
    X.s <- X[sampling,]
    Y.s <- Y[sampling,]
    M.s <- M[sampling]
    Y.est.s <- Y.s[,ind.inside]
    X_tilde.s <- M.s * X.s
    Y_tilde.s <- M.s * Y.est.s
    
    # Distributed learning
    fit.boot <- mclapply(idx.sample.tri, FUN = local.fit.bivar, mc.cores = core,
                         V0 = Ver, Tr0 = Tr, TV0 = TV, n.layer = n.layer, 
                         Y.all = Y_tilde.s, X.all = X_tilde.s,
                         Z.all = V, d = d, r = r, lambda = 10^(seq(-6, 0, by = 1)))
    gamma.all.b <- matrix(0, ncol = ncol(X_tilde.s), nrow = nrow(Tr) * (d+1) * (d+2)/2)
    count.tri.b <- rep(0, nrow(Tr))
    
    for(iter in 1:length(fit.boot)) {
      idx.tr.b <- fit.boot[[iter]]$idx.tr
      count.tri.b[idx.tr.b] <- count.tri.b[idx.tr.b] + 1
      gamma.all.b[fit.boot[[iter]]$idx.gamma, ] = 
        gamma.all.b[fit.boot[[iter]]$idx.gamma, ] + fit.boot[[iter]]$gamma.local
    }
    
    count.gamma.b <- rep(count.tri.b, each = nbasis.tri)
    mean.gamma.b <- 1/count.gamma.b * gamma.all.b
    
    #Ball.b = basis(Ver,Tr,d,r,V)
    theta.hat.b <- solve(crossprod(Ball$Q2), crossprod(Ball$Q2, mean.gamma.b))
    bivar.hat.b <- Ball$B %*% Ball$Q2 %*% theta.hat.b
    boot_beta[,,b] <- as.matrix(bivar.hat.b)
  }
  # Obtain bootstrap variance matrix
  variance.mat <- apply(boot_beta, c(1,2), var, na.rm = TRUE)
  se.mat <- apply(boot_beta, c(1,2), sd, na.rm = TRUE)
  # Obtain coverage probability (95%)
  coverage.prob.lower_95 <- bivar.est + qnorm(0.05/2)*se.mat
  coverage.prob.upper_95 <- bivar.est + qnorm(1-0.05/2)*se.mat
  covered_95 <- (bivar.true[ind.inside,] >= coverage.prob.lower_95) & (bivar.true[ind.inside,] <= coverage.prob.upper_95)
  covered.mean_95 <- colMeans(covered_95)
  
  # Calculate Bias
  bias.res <- apply((bivar.est - bivar.true[ind.inside,]), 2, mean, na.rm = TRUE)
  # Calculate MISE
  mse.res <- apply((bivar.est - bivar.true[ind.inside,])^2, 2, mean, na.rm = TRUE)
  return(list(bias = bias.res, MISE = mse.res,
              coverage.rate.95=covered.mean_95, 
              covered_95_mat = covered_95,
              variance = colMeans(variance.mat)))
}

nT <- 6
n.samp <- 20
n.layer <- 3
core <- 2
n <- 100; d <- 2; r <- 1
nboot<- 100

# Define 2D spatial grid and extract coordinate vectors
xm <- seq(-1, 3.5, length = 101)
yn <- seq(-1, 1, length = 101)
xy_grid <- pracma::meshgrid(xm, yn)
uu <- c(xy_grid$X)
vv <- c(xy_grid$Y)
V.all <- as.matrix(cbind(uu,vv))

# Load horseshoe domain and construct triangulation mesh
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
round(colMeans(var_mat), 4)
round(colMeans(cr_95_map), 4)

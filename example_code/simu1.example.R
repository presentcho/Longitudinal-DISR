rm(list = ls())
library(BPST)
library(TPST)
library(dplyr)
library(Triangulation)
library(pracma)
library(MASS)
library(mgcv)
library(colorRamps)
library(ggplot2)
library(splines2)
require(MGLM)
library(gghilbertstrings)
library(parallel)
source('../funcs/bivar.gc.fit.R')
source('../funcs/ring.dc.R')
source('../funcs/sampling.HC.R')
source('../funcs/local.fit.bivar.R')
source('../funcs/plot.fun.R')
source('../funcs/simul1.data.gen.R')

# Initialize model parameters (Sample size, degree, and smoothness)
nT <- 6
n <- 100; d <- 2; r <- 1;

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

# Generate the dataset and true coefficient parameters
dat <- simul.data.gen(n, Tr, Ver, V.all)
X <- dat$X
Y <- dat$Y
M <- dat$M.vec
tij <- dat$tij
ind.inside <- dat$ind.inside

# Visualize the true bivariate coefficient functions
alpha.true <- dat$bivar.alpha
beta.true <- dat$bivar.beta
plot.fun(uu, vv, alpha.true, beta.true, hs, ind.inside)

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

# Fit the global model to estimate bivariate coefficient functions
X_tilde <- as.matrix(M * X)
Y_tilde <- as.matrix(M * Y.est)
mfit0 <- bivar.gc.fit(B, Q2, K, lambda, X_tilde, Y_tilde)
global.bivar.est <- mfit0$beta

# Visualize the estimated global coefficient functions
global.alpha.est <- global.bivar.est[,1:3]
global.beta.est <- global.bivar.est[,4:6]
plot.fun(uu, vv, global.alpha.est, global.beta.est, hs, ind.inside)

# Calculate Mean Integrated Squared Error (MISE) for the global method
global.res <- round(c(
  apply((as.matrix(global.alpha.est) - as.matrix(alpha.true[ind.inside,]))^2, 2, mean), 
  apply((as.matrix(global.beta.est) - as.matrix(beta.true[ind.inside,]))^2, 2, mean)), 3)

print(global.res)

# Initialize parameters and sample subregions for DDC
n.samp <- 20; n.layer <- 3; core <- 2
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
HD.bivar.est = Ball$B %*% Ball$Q2 %*% theta.hat

# Visualize the final aggregated HD coefficient functions
HD.alpha.est <- HD.bivar.est[,1:3]
HD.beta.est <- HD.bivar.est[,4:6]
plot.fun(uu, vv, HD.alpha.est, HD.beta.est, hs, ind.inside)

# Calculate Mean Integrated Squared Error (MISE) for the HD method
HD.res <- round(c(
  apply((as.matrix(HD.alpha.est) - as.matrix(alpha.true[ind.inside,]))^2, 2, mean), 
  apply((as.matrix(HD.beta.est) - as.matrix(beta.true[ind.inside,]))^2, 2, mean)), 3)

print(HD.res)


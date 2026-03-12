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
library(ggplot2)
library(MGLM)
library(parallel)
library(colorRamps)
library(gghilbertstrings)
source('../funcs/simul2.data.gen.R')
source('../funcs/basis.tensor.R')
source('../funcs/energy.tensor.R')
source('../funcs/gc.fit.R')
source('../funcs/gc.fit.ho.R')
source('../funcs/plot.fun.R')
source('../funcs/gc.fit.dc.R')
source('../funcs/sampling.HC.R')
source('../funcs/ring.dc.R')
source('../funcs/basis.tensor.local.R')
source('../funcs/local.fit.R')
source('../funcs/gc.fit.cv.R')

# Initialize model parameters (Sample size, degree, and smoothness)
L <- 3
rho <- 3
nT <- 6
n <- 100; d <- 2; r <- 1

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

# Load simulated datasets and true coefficient parameters
Y <- as.matrix(read.csv('../data/simul2/Y.csv'))
X <- as.matrix(read.csv('../data/simul2/X.csv'))
V <- as.matrix(read.csv('../data/simul2/V.csv'))
M <- as.vector(read.csv('../data/simul2/M.csv')$x)
tij <- as.vector(read.csv('../data/simul2/tij.csv')$x)
ind.inside <- as.vector(read.csv('../data/simul2/ind.inside.csv')$x)
bivar.alpha.true <- as.matrix(read.csv('../data/simul2/bivar.alpha.csv'))
bivar.beta.true <- as.matrix(read.csv('../data/simul2/bivar.beta.csv'))
trivar.true <- as.matrix(read.csv('../data/simul2/tri.alpha.csv'))

# Visualize the true coefficient functions
alpha.true <- bivar.alpha.true
beta.true <- bivar.beta.true
plot.fun(uu, vv, alpha.true, beta.true, hs, ind.inside)

# Construct tensor product basis and penalty matrices for time and spatial domains
time.bound<- c(min(tij), max(tij))
probs <- seq(0, 1, length.out = L + 2)
time.knots <- quantile(tij, probs = probs)[-c(1, L + 2)]
Basis <- basis.tensor(ss = V, tt = tij, V = Ver, Tri = Tr, 
                      d = d, r = r, time.knots = time.knots, time.bound = time.bound, rho = rho)

U0 <- Basis$U0
B0 <- Basis$B0
Q2 <- Basis$Q2
K <- Basis$K
P11 <- Basis$D1
P12 <- Basis$D2
P2 <- kronecker(diag(ncol(X)), as.matrix(crossprod(Q2, K) %*% Q2))

# Find optimal tuning parameters via Cross-Validation 
cv.fit <- gc.fit.cv(nfold=5,X=X, Y=Y, V=V, M=M, tij=tij)
lambda11.optimal <- cv.fit$lambda11
lambda12.optimal <- cv.fit$lambda12
lambda2.optimal <- cv.fit$lambda2

# Fit the global model to estimate proposed coefficient functions
res <- gc.fit(X, Y, V, M, U0, B0, Q2, K, P11, P12, P2, L, rho,
              lambda11=lambda11.optimal, lambda12=lambda12.optimal, lambda2=lambda2.optimal)
global.bivar.est <- res$bivar.est
global.trivar.est <- res$trivar.est

# Visualize the estimated global coefficient functions
global.alpha.est <- global.bivar.est[,1:3]
global.beta.est <- global.bivar.est[,4:6]
plot.fun(uu, vv, global.alpha.est, global.beta.est, hs, ind.inside)

# Calculate Mean Integrated Squared Error (MISE) for the global method
global.trivar.all <- rowMeans((global.trivar.est - trivar.true)^2, na.rm = TRUE)
global.res <- round(c(mean(global.trivar.all), apply((global.alpha.est-alpha.true)^2, 2, mean), apply((global.beta.est-beta.true)^2, 2, mean)), 3)
names(global.res)[1] <- 'alpha0'
print(global.res)

# Initialize parameters and sample sub-domains for Divide-and-Conquer (DDC)
n.samp <- 20; n.layer <- 3; n.core <- 2

sampling.res <- sampling.HC(n.samp, Ver, Tr, n.layer)
count.tri <- rep(0, nrow(Tr))
TV <- tdata(Ver, Tr)$TV
idx.sample.tri <- sampling.res$sample.tri

# Fit the HD model to estimate proposed coefficient functions
fit.dc <- gc.fit.dc(idx.sample.tri = idx.sample.tri, Tr = Tr, Ver = Ver, TV = TV, 
                    n.layer = n.layer, X = X, Y = Y, V = V, M= M, tij = tij,
                    d=d, r=r, L=L, rho=rho, 
                    lambda11=lambda11.optimal, lambda12 = lambda12.optimal, lambda2= lambda2.optimal)
HD.bivar.est <- as.matrix(fit.dc$bivar.est)
HD.trivar.est <- as.matrix(fit.dc$trivar.est)

# Visualize the estimated HD coefficient functions
HD.alpha.est <- HD.bivar.est[,1:3]
HD.beta.est <- HD.bivar.est[,4:6]
plot.fun(uu, vv, HD.alpha.est, HD.beta.est, hs, ind.inside)

# Calculate Mean Integrated Squared Error (MISE) for the HD method
HD.trivar.all <- rowMeans((HD.trivar.est - trivar.true)^2, na.rm = TRUE)
HD.res <- round(c(mean(HD.trivar.all), apply((HD.alpha.est-alpha.true)^2, 2, mean), apply((HD.beta.est-beta.true)^2, 2, mean)), 3)
names(HD.res)[1] <- 'alpha0'
print(HD.res)

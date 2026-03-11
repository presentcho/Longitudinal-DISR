rm(list=ls())
library(BPST)
library(TPST)
library(dplyr)
library(grpnet)
library(Triangulation)
library(pracma)
library(MASS)
library(mgcv)
library(colorRamps)
library(ggplot2)
library(splines2)
require(MGLM)
source('funcs/basis.tensor.R')
source('funcs/energy.tensor.R')
source('funcs/simul2.data.gen.R')
source('funcs/trivar.fit.R')
source('funcs/trivar.fit.ho.R')
source('funcs/trivar.fit.cv.R')
source('funcs/gc.fit.cv.R')
source('funcs/gc.fit.ho.R')

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
Y <- as.matrix(read.csv('data/simul2/Y.csv'))
X <- as.matrix(read.csv('data/simul2/X.csv'))
V <- as.matrix(read.csv('data/simul2/V.csv'))
M <- as.vector(read.csv('data/simul2/M.csv')$x)
tij <- as.vector(read.csv('data/simul2/tij.csv')$x)

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

# Find optimal tuning parameters for the Baseline Model (Section 5) via Cross-Validation
X.tilde <- cbind(1, X[,1:3])
tuning.param <- trivar.fit.cv(nfold=5, X.tilde, as.matrix(Y), V, tij, U0, B0, Q2, P11, P12)
lambda11.optimal <- tuning.param$lambda11
lambda12.optimal <- tuning.param$lambda12

# Evaluate Baseline Model (Section 5) performance via Train-Test Split
n <- nrow(Y); npix <- ncol(Y) 
train.idx <- sample(1:n, size = floor(0.7*n))

# Prepare training dataset for the Baseline Model
set.seed(123)
X.train <- as.matrix(X.tilde[train.idx,])
Y.train <- as.matrix(Y[train.idx,])
U0.train <- U0[train.idx,]
test.idx <- setdiff(c(1:n), train.idx)

# Prepare testing dataset for the Baseline Model
X.test <- as.matrix(X.tilde[test.idx,])
Y.test <- as.matrix(Y[test.idx,])
U0.test <- U0[test.idx,]

# Fit the Baseline Model and calculate SSE/MSE
cv.res <- trivar.fit.ho(X.train = X.train, Y.train = Y.train, U0.train = U0.train,
                        X.test = X.test, Y.test = Y.test, U0.test = U0.test, 
                        B0 = B0, Q2 = Q2, P11 = P11, P12 = P12,
                        lambda11 = lambda11.optimal, 
                        lambda12 = lambda12.optimal)
cv.res.Yhat <- as.matrix(cv.res$Yhat)
cv.res$sse
cv.res$mse

# Find optimal tuning parameters for the Proposed Model via Cross-Validation
cv.fit <- gc.fit.cv(nfold = 5, X=X, Y=Y, V=V, M=M, tij=tij)
lambda11.optimal <- cv.fit$lambda11
lambda12.optimal <- cv.fit$lambda12
lambda2.optimal <- cv.fit$lambda2

# Prepare training dataset for the Proposed Model
X.train <- X[train.idx,]
Y.train <- Y[train.idx,]
U0.train <- U0[train.idx,]
M.train <- M[train.idx]
test.idx <- setdiff(c(1:n), train.idx)

# Prepare testing dataset for the Proposed Model
X.test <- X[test.idx,]
Y.test <- Y[test.idx,]
U0.test <- U0[test.idx,]
M.test <- M[test.idx]

# Fit the Proposed Model and calculate SSE/MSE
gc.fit.cv <- gc.fit.ho(X.train, Y.train, M.train, U0.train,
                       X.test, Y.test, M.test, U0.test, 
                       V, B0, Q2, K, P11, P12, P2, lambda11.optimal, lambda12.optimal, lambda2.optimal)
gc.fit.cv$sse
gc.fit.cv$mse

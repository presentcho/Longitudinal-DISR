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
source('funcs/simul1_data_gen.R')
source('funcs/bivar.gc.fit.R')
source('funcs/ring.dc.R')
source('funcs/sampling.HC.R')
source('funcs/local.fit.bivar.R')

################################################################################
simulation <- function(n.sim){
  MISE <- list()
  bias <- list()
  for (i in 1:n.sim){
    set.seed(i)
    # Set boundary 
    xm <- seq(-1, 3.5, length = 101)
    yn <- seq(-1, 1, length = 101)
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
    dat <- simul.data.gen(n, Tr, Ver, V)
    Y <- dat$Y 
    X <- dat$X 
    M <- dat$M 
    V <- dat$V
    ind.inside <- dat$ind.inside
    alpha.true <- dat$bivar.alpha
    beta.true <- dat$bivar.beta
    bivar.true <- cbind(alpha.true, beta.true)
    
    Y.est <- Y[,ind.inside]
    V <- V[ind.inside, ]
    
    # basis function
    if(d > 1){
      Ball <- basis(Ver, Tr, d, r, V)
      K <- Ball$K
      Q2 <- Ball$Q2
      B <- Ball$B 
    }
    
    X_tilde <- M * X
    Y_tilde <- M * Y.est
    
    # Hilbert Curve sampling 
    result <- sampling.HC(n.samp, Ver, Tr, n.layer = n.layer)
    gamma.all <- matrix(0, ncol = ncol(X_tilde), nrow = nrow(Tr) * (d+1) * (d+2)/2)
    count.tri <- rep(0, nrow(Tr))
    nbasis.tri <- (d+1) * (d+2)/2
    TV <- tdata(Ver, Tr)$TV
    idx.sample.tri <- result$sample.tri
    
    # distributed learning
    fit.all <- mclapply(idx.sample.tri, FUN = local.fit.bivar, mc.cores = core,
                        V0 = Ver, Tr0 = Tr, TV0 = TV, n.layer = n.layer, 
                        Y.all = Y_tilde, X.all = X_tilde,
                        Z.all = V, d = 5, r = 1, lambda = 10^(seq(-6, 0, by = 1)))
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
  
    MISE[[i]] <- apply((bivar.est - bivar.true[ind.inside,])^2, 2, mean, na.rm = TRUE)
    bias[[i]] <- apply((bivar.est - bivar.true[ind.inside,]), 2, mean, na.rm = TRUE)
  }
  MISE.res <- round(colMeans(do.call(rbind, MISE)), 4)
  bias.res <- round(colMeans(do.call(rbind, bias)), 4)
  return(list(bias = bias.res, MISE = MISE.res))
}


nT <- 6
n.samp <- 20
n.layer <- 2
core <- 1
n <- 100; d <- 5; r <- 1
simulation(100)


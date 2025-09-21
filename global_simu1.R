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
source('funcs/simul1_data_gen.R')
source('funcs/bivar.gc.fit.R')

################################################################################
simulation <- function(iter){
  MISE <- list()
  bias <- list()
  for (i in 1:iter){
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
    V <- dat$V
    M <- dat$M
    ind.inside <- dat$ind.inside
    alpha.true <- dat$bivar.alpha
    beta.true <- dat$bivar.beta
    bivar.true <- cbind(alpha.true, beta.true)
    
    lambda <- 10^(seq(-6, 0, by = 1))
    
    Y.est <- Y[,ind.inside]
    V <- V[ind.inside, ]
    
    # image on regression 
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
    bias.res <- apply((bivar.est - bivar.true[ind.inside,]), 2, mean, na.rm = TRUE)
    mse.res <- apply((bivar.est - bivar.true[ind.inside,])^2, 2, mean, na.rm = TRUE)
    MISE[[i]] <- mse.res
    bias[[i]] <- bias.res
  }
  MISE.res <- round(colMeans(do.call(rbind, MISE)), 4)
  bias.res <- round(colMeans(do.call(rbind, bias)), 4)
  return(list(bias = bias.res, MISE = MISE.res))
}


nT <- 6
n <- 100; d <- 5; r <- 1
simulation(100)


# gc.fit.dc: fits estimators through a parallel algorithm based on domain decomposition

gc.fit.dc <- function(n.samp, n.layer, Ver, Tr, V, X, Y, 
                      d, r, L, rho, tij, M,
                      n.core, lambda11 = exp(seq(log(0.001), log(1000), length.out = 5)), 
                      lambda12 = exp(seq(log(0.001), log(1000), length.out = 5)), 
                      lambda2 = exp(seq(log(0.001), log(1000), length.out = 5))){
  n.pix <- nrow(V)
  time.bound <- c(min(tij), max(tij))
  probs <- seq(0, 1, length.out = L + 2)
  time.knots <- quantile(tij, probs = probs)[-c(1, (L+2))]
  Basis <- basis.tensor(ss = V, tt = tij, V = Ver, Tri = Tr, 
                        d = d, r = r, time.knots = time.knots, time.bound = time.bound, rho = rho)
  B0 <- Basis$B0
  U0 <- Basis$U0
  Q2 <- Basis$Q2
  K <- Basis$K
  P11 <- Basis$D1
  P12 <- Basis$D2
  P2 <- kronecker(diag(ncol(X)), as.matrix(crossprod(Q2, K) %*% Q2))
  ind.inside <- Basis$ind.inside
  V <- V[ind.inside, ]
  Y.est <- Y[, ind.inside]
  B0 <- B0[ind.inside,]
  
  sampling.res <- sampling.HC(n.samp, Ver, Tr, n.layer)
  count.tri <- rep(0, nrow(Tr))
  TV <- tdata(Ver, Tr)$TV
  idx.sample.tri <- sampling.res$sample.tri
  
  fit.all <- mclapply(idx.sample.tri, function(iter) {
    local.fit(iter, Ver0 = Ver, Tr0 = Tr, TV0 = TV, n.layer = n.layer, 
              X = X, Y = Y.est, V = V, M = M, tij = tij, 
              d = d, r = r, L = L, rho = rho, 
              time.knots = time.knots, time.bound = time.bound,
              lambda11 = exp(seq(log(0.001), log(1000), length.out = 5)), 
              lambda12 = exp(seq(log(0.001), log(1000), length.out = 5)), 
              lambda2 = exp(seq(log(0.001), log(1000), length.out = 5)))
  }, mc.cores = n.core, mc.preschedule = TRUE)
  
  gamma.all <- matrix(0, ncol = L+rho, nrow = nrow(Tr)*(d+1)*(d+2)/2)
  theta.all <- matrix(0, ncol = ncol(X), nrow = nrow(Tr)*(d+1)*(d+2)/2)
  
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
  
  p = dim(X)[2]
  Basis <- basis.tensor(ss = V, tt = tij, V = Ver, Tri = Tr,
                        d = d, r = r, time.knots = time.knots, time.bound = time.bound, rho = rho)
  
  est.gamma = solve(crossprod(Basis$Q2), crossprod(Basis$Q2, mean.gamma))
  est.theta = solve(crossprod(Basis$Q2), crossprod(Basis$Q2, mean.theta))
  
  alpha0.hat <- tcrossprod(Basis$U0, (Basis$B0 %*% Basis$Q2 %*% est.gamma)) 
  bivar.est <- Basis$B0 %*% Basis$Q2 %*% est.theta 
  
  beta.hat <- matrix(NA, n.pix, dim(X)[2])
  beta.hat[ind.inside,] <- as.matrix(bivar.est) 
  
  trivar.hat <- matrix(NA, n.pix, dim(X)[1])
  trivar.hat[ind.inside,] <- t(as.matrix(alpha0.hat))
  
  Yhat <- t(trivar.hat) + tcrossprod(X, beta.hat) 
  Rhat <- Y - Yhat
  
  return(list(trivar.hat = trivar.hat, beta.hat=beta.hat, gamma = est.gamma, theta=est.theta, Yhat=Yhat, Rhat=Rhat))
}

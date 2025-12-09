# gc.fit.dc.k: fits estimators through a parallel algorithm based on domain decomposition

gc.fit.dc.k <- function(idx.sample.tri, Tr, Ver, TV, n.layer, 
                        X, Y, V, M, tij, d, r, L, rho, 
                        lambda11, lambda12, lambda2, Basis){
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
  
  fit.all <- mclapply(idx.sample.tri, function(iter) {
    local.fit.k(iter, Ver0 = Ver, Tr0 = Tr, TV0 = TV, n.layer = n.layer, 
                X = X, Y = Y, V = V, M = M, tij = tij,
                d = d, r = r, L = L, rho = rho, 
                time.knots = time.knots, time.bound = time.bound, 
                lambda11=lambda11, lambda12=lambda12, lambda2=lambda2)
  }, mc.cores = n.core, mc.preschedule = TRUE)
  
  count.tri <- rep(0, nrow(Tr))
  
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
  est.gamma = solve(crossprod(Q2), crossprod(Q2, mean.gamma))
  est.theta = solve(crossprod(Q2), crossprod(Q2, mean.theta))
  
  trivar.est <- tcrossprod(U0, (B0 %*% Q2 %*% est.gamma)) 
  bivar.est <- B0 %*% Q2 %*% est.theta 
  
  return(list(trivar.est = trivar.est, bivar.est=bivar.est))
}

trivar.fit <- function(X, Y, U0, B0, Q2, P11, P12, L, rho,
                       lambda11, lambda12){
  X <- as.matrix(X)
  U0 <- as.matrix(U0)
  B.tilde <- as.matrix(B0 %*% Q2)
  XU <- row.kronecker(X, U0)
  XUtXU <- crossprod(XU)
  BtB <- crossprod(as(B.tilde, 'dgCMatrix'))
  DD = kronecker(XUtXU, BtB)
  P.lambda <- Diagonal(ncol(X)) %x% bdiag(lambda11*P11 + lambda12*P12)
  lhs <- DD + P.lambda
  rhs <- rowSums(sapply(1:nrow(Y), function(iter)
    as.matrix(kronecker(XU[iter,], crossprod(B.tilde, Y[iter,])))))
  gamma.hat <- solve(lhs, rhs)
  Q2.all <- kronecker(diag(L + rho), Q2)
  gamma.mat <- matrix(gamma.hat, nrow = ncol(Q2.all), ncol = ncol(X))
  trivar.est <- array(NA_real_, dim=c(nrow(B.tilde), nrow(U0), ncol(X)))
  for (k in seq_len(ncol(X))){
    #theta_k <- matrix(gamma.mat[,k], ncol = ncol(B.tilde), nrow = ncol(U0))
    theta_k <- gamma.mat[,k]
    trivar.est[,,k] <- sapply(1:nrow(U0), function(iter) kronecker(as.matrix(t(U0[iter,])), B.tilde) %*% theta_k)
    #trivar.est[,,k] <- as.matrix(B.tilde %*% t(theta_k) %*% t(U0))
  }
  return(list(trivar.est = trivar.est, gamma.hat = gamma.mat))
}

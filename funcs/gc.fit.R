# gc.fit: fits estimators through global methods 

gc.fit <- function(X, Y, V, M, U0, B0, Q2, K, P11, P12, P2, L, rho,
                    lambda11, lambda12, lambda2){
  B.tilde <- as.matrix(B0 %*% Q2)
  X.tilde <- M * X
  Y.tilde <- M * Y
  U.tilde <- M * U0
  D <- cbind(U.tilde, X.tilde)
  DD <- kronecker(crossprod(D), crossprod(B.tilde))
  P.lambda <- bdiag(lambda11*P11 + lambda12*P12, lambda2*P2)
  lhs <- DD + P.lambda
  rhs1 <- rowSums(sapply(1:nrow(Y.tilde), function(iter)
    as.matrix(kronecker(U.tilde[iter,], crossprod(B.tilde, Y.tilde[iter,])))))
  rhs2 <- rowSums(sapply(1:nrow(Y.tilde), function(iter)
    as.matrix(kronecker(X.tilde[iter,],crossprod(B.tilde, Y.tilde[iter,])))))
  rhs <- c(rhs1, rhs2)
  phi.est <- solve(lhs, rhs)
  gamma.hat <- phi.est[1:length(rhs1)]
  Q2.all <- kronecker(diag(L + rho), Q2)
  gamma.mat <- matrix(as.vector(Q2.all %*% gamma.hat), nrow = nrow(Q2), ncol = L+rho) 
  theta.hat <- phi.est[(length(rhs1)+1):length(rhs)]
  theta.mat <- Q2 %*% matrix(theta.hat, nrow = ncol(Q2), ncol = ncol(X))
  trivar.est <- t(sapply(1:nrow(U0), function(iter) kronecker(as.matrix(t(U0[iter,])), B.tilde) %*% gamma.hat))

  bivar.est <- B.tilde %*% matrix(theta.hat, nrow = ncol(Q2), ncol = ncol(X))
  return(list(trivar.est = trivar.est, bivar.est = bivar.est,
              gamma.hat = gamma.mat, theta.hat = theta.mat))
}

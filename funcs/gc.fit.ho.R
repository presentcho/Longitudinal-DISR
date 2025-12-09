# gc.fit.ho: obtain sse/mse using train and test dataset

gc.fit.ho <- function(X.train, Y.train, M.train, U0.train,
                      X.test, Y.test, M.test, U0.test, 
                      V, B0, Q2, K, P11, P12, P2, lambda11, lambda12, lambda2){
  # train fit 
  B.tilde <- as.matrix(B0 %*% Q2)
  X.tilde <- M.train * X.train
  Y.tilde <- M.train * Y.train
  U.tilde <- M.train * U0.train
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
  r.hat <- phi.est[1:length(rhs1)]
  theta.hat <- phi.est[(length(rhs1)+1):length(rhs)]
  bivar.est <- B.tilde %*% matrix(theta.hat, nrow = ncol(Q2), ncol = ncol(X.train))
  
  # test fit
  alpha0.test <- kronecker(U0.test, B.tilde) %*% r.hat
  trivar.test <- t(matrix(alpha0.test, ncol = nrow(U0.test), nrow = nrow(B0)))
  Yhat.test <- trivar.test + X.test %*% t(bivar.est)
  sse.test <- sum(apply((Y.test - Yhat.test)^2, 1, sum))
  mse.test <- mean(apply((Y.test - Yhat.test)^2, 1, sum))
  return(list(sse = sse.test, mse = mse.test, Yhat = Yhat.test))
}

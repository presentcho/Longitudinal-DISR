trivar.fit.ho <- function(X.train, Y.train, U0.train,
                      X.test, Y.test, U0.test, 
                      B0, Q2, P11, P12, lambda11, lambda12, lambda2){
  # train fit
  B.tilde <- B0 %*% Q2
  fit0 <- trivar.fit(X.train, Y.train, U0.train, B0, Q2, P11, P12, L, rho, lambda11, lambda12)
  trivar.train <- fit0$trivar.est
  gamma.train <- fit0$gamma.hat
  
  # test fit
  trivar.test <- array(NA_real_, dim=c(nrow(B.tilde), nrow(U0.test), ncol(X.test)))
  for (k in seq_len(ncol(X.test))){
    theta_k <- as.vector(gamma.train[,k])
    trivar.test[, , k] <- sapply(seq_len(nrow(U0.test)), function(iter) {
      Ki <-  kronecker(as.matrix(t(U0[iter,])), B.tilde)
      as.vector(Ki %*% theta_k)
    })  }
  
  Yhat.test <- Reduce('+', lapply(1:ncol(X.test), function(k){
    sweep(t(trivar.test[,,k]), 1, X.test[,k], `*`)
  }))
  sse.test <- sum(apply((Y.test - Yhat.test)^2, 1, sum))
  mse.test <- mean(apply((Y.test - Yhat.test)^2, 1, sum))
  return(list(sse = sse.test, mse = mse.test, Yhat = Yhat.test, trivar = trivar.train))
}

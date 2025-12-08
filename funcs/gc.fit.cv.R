# gc.fit.gcv: finds optimal tuning parameters using cross-valdation

gc.fit.cv <- function(nfold, X, Y, V, M, tij, Basis, lambda11 = exp(seq(log(0.001), log(1000), length.out = 5)),
                      lambda12 = exp(seq(log(0.001), log(1000), length.out = 5)),
                      lambda2 = 10^seq(-3,3,by=0.5), iter = 2024){
  
  U0 <- Basis$U0
  B0 <- Basis$B0
  Q2 <- Basis$Q2
  K <- Basis$K
  P11 <- Basis$D1
  P12 <- Basis$D2
  P2 <- kronecker(diag(Basis$dimU), as.matrix(crossprod(Q2, K) %*% Q2))
  
  cv.function <- function(find.optimal, fix.lambda11, fix.lambda12, fix.lambda2){
    set.seed(iter)
    n <- nrow(Y)
    npix <- ncol(Y)
    sfold <- round(n/nfold)
    idx.sample <- sample(1:n)
    cv.error <- list()
    for(ii in 1:nfold){
      if(ii < nfold){
        test.idx <- sort(idx.sample[((ii-1)*sfold+1):(ii*sfold)])
      } else {
        test.idx <- sort(idx.sample[((ii-1)*sfold+1):n])
      }
      train.idx <- setdiff(1:n, test.idx)
      X.train <- X[train.idx,]
      Y.train <- Y[train.idx,]
      M.train <- M[train.idx]
      U0.train <- U0[train.idx,]
      
      # Test set
      X.test <- X[test.idx,]
      Y.test <- Y[test.idx,]
      M.test <- M[test.idx]
      U0.test <- U0[test.idx,]
      
      sse.lambda <- sapply(find.optimal, function(lambda) {
      res <- gc.fit.ho(X.train = X.train, Y.train = Y.train, M.train = M.train, U0.train = U0.train,
                       X.test = X.test, Y.test = Y.test, M.test = M.test, U0.test = U0.test, 
                       V = V, B0 = B0, Q2 = Q2, K = K, P11 = P11, P12 = P12, P2 = P2,
                       lambda11 = ifelse(is.null(fix.lambda11), lambda, fix.lambda11), 
                       lambda12 = ifelse(is.null(fix.lambda12), lambda, fix.lambda12), 
                       lambda2 = ifelse(is.null(fix.lambda2), lambda, fix.lambda2))
        return(res$sse)
      })
      cv.error[[ii]] <- unlist(sse.lambda)
    }
    lambda.cv.error <- do.call(rbind, cv.error)
    lambda.min <- which.min(apply(lambda.cv.error, 2, sum))
    optimal.lambda <- find.optimal[lambda.min]
    lambda.mspe <- apply(lambda.cv.error, 2, sum) / n / npix
    return(list(optimal.lambda = optimal.lambda, lambda.mspe = lambda.mspe))
  }
  
  lambda11.cv <- cv.function(lambda11, NULL, 1e-3, 1e-3)
  optimal.lambda11 <- lambda11.cv$optimal.lambda
  lambda12.cv <- cv.function(lambda12, optimal.lambda11, NULL, 1e-3)
  optimal.lambda12 <- lambda12.cv$optimal.lambda
  lambda2.cv <- cv.function(lambda2, optimal.lambda11, optimal.lambda12, NULL)
  optimal.lambda2 <- lambda2.cv$optimal.lambda
  
  return(list(lambda11 = optimal.lambda11, lambda11.mspe = lambda11.cv$lambda.mspe,
              lambda12 = optimal.lambda12, lambda12.mpse = lambda12.cv$lambda.mspe,
              lambda2 = optimal.lambda2, lambda2.mspe = lambda2.cv$lambda.mspe))
}

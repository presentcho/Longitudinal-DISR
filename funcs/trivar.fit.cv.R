trivar.fit.cv <- function(nfold, X, Y, V, tij, U0, B0, Q2, P11, P12,
                          lambda11 = exp(seq(log(0.001), log(1000), length.out = 5)),
                       lambda12 = exp(seq(log(0.001), log(1000), length.out = 5)), iter=2024){
  
    cv.function <- function(find.optimal, fix.lambda11, fix.lambda12){
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
      U0.train <- U0[train.idx,]
      
      # Test set
      X.test <- X[test.idx,]
      Y.test <- Y[test.idx,]
      U0.test <- U0[test.idx,]
      
      sse.lambda <- sapply(find.optimal, function(lambda) {
        res <- trivar.fit.ho(X.train = X.train, Y.train = Y.train, U0.train = U0.train,
                         X.test = X.test, Y.test = Y.test, U0.test = U0.test, 
                         B0 = B0, Q2 = Q2, P11 = P11, P12 = P12,
                         lambda11 = ifelse(is.null(fix.lambda11), lambda, fix.lambda11), 
                         lambda12 = ifelse(is.null(fix.lambda12), lambda, fix.lambda12))
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
  
  lambda11.cv <- cv.function(lambda11, NULL, 1e-3)
  optimal.lambda11 <- lambda11.cv$optimal.lambda
  lambda12.cv <- cv.function(lambda12, optimal.lambda11, NULL)
  optimal.lambda12 <- lambda12.cv$optimal.lambda
  
  return(list(lambda11 = optimal.lambda11, lambda11.mspe = lambda11.cv$lambda.mspe,
              lambda12 = optimal.lambda12, lambda12.mpse = lambda12.cv$lambda.mspe))
}

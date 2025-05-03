# gc.fit.gcv: finds optimal tuning parameters using GCV

gc.fit.gcv <- function(X, Y, V, M, U0, B0, Q2, K, P11, P12, P2, 
                       lambda11, lambda12, lambda2){
  B.tilde <- as.matrix(B0 %*% Q2)
  X.tilde <- M * X
  Y.tilde <- M * Y
  U.tilde <- M * U0
  D <- cbind(U.tilde, X.tilde)
  DD <- kronecker(crossprod(D), crossprod(B.tilde))
  n <- nrow(Y)
  npix <- ncol(Y)
  
  rhs1 <- rowSums(sapply(1:nrow(Y.tilde), function(iter)
    as.matrix(kronecker(U.tilde[iter,], crossprod(B.tilde, Y.tilde[iter,])))))
  rhs2 <- rowSums(sapply(1:nrow(Y.tilde), function(iter)
    as.matrix(kronecker(X.tilde[iter,],crossprod(B.tilde, Y.tilde[iter,])))))
  eps <- 1e-6
  DD_reg <- DD + diag(eps, nrow(DD))
  Ainv <- chol(DD_reg, pivot = TRUE) 
  A <- solve(t(Ainv))
  D <- bdiag(P11+ P12, P2) 
  ADA <- A %*% D %*% t(A)
  eigs <- eigen(ADA)
  C <- eigs$values 
  
  optimize_lambda <- function(lambda_vals, fixed_lambda11, fixed_lambda12, fixed_lambda2) {
    sse_all <- c()
    gcv_all <- c()
    
    for (il in 1:length(lambda_vals)) {
      lam11 <- if (is.null(fixed_lambda11)) lambda_vals[il] else fixed_lambda11
      lam12 <- if (is.null(fixed_lambda12)) lambda_vals[il] else fixed_lambda12
      lam2 <- if (is.null(fixed_lambda2)) lambda_vals[il] else fixed_lambda2
      
      lambda1_mat <- Matrix(lam11 * P11 + lam12 * P12, sparse = TRUE)
      lambda2_mat <- Matrix(lam2 * P2, sparse = TRUE)
      P.lambda <- Matrix::bdiag(lambda1_mat, lambda2_mat)
      lhs <- DD + P.lambda
      rhs <- c(rhs1, rhs2)
      phi.est <- solve(lhs, rhs)
      gamma.hat <- phi.est[1:length(rhs1)]
      theta.hat <- phi.est[(length(rhs1) + 1):length(rhs)]
      trivar.est <- t(sapply(1:nrow(U0), function(iter) kronecker(as.matrix(t(U0[iter,])), B.tilde) %*% gamma.hat))
      bivar.est <- B.tilde %*% matrix(theta.hat, nrow = ncol(Q2), ncol = ncol(X))
      Yhat <- trivar.est + X %*% t(bivar.est)
      
      sseij <- apply((Y - Yhat)^2, 1, sum)
      sse <- sum(sseij)
      sse_all[il] <- sse
      df <- sum(1/(1+C*lambda_vals[il]))
      gcv <- n*npix*sse/(n*npix-df)^2
      gcv_all[il] <- gcv
    }
    optimal_lambda <- lambda_vals[which.min(gcv_all)]
    return(list(optimal_lambda = optimal_lambda, gcv_all = gcv_all, sse_all = sse_all))
  }
  
  lambda11_result <- optimize_lambda(lambda11, fixed_lambda11 = NULL, fixed_lambda12 = 1e-3,
                                     fixed_lambda2 = 1e-3)
  optimal.lambda11 <- lambda11_result$optimal_lambda
  lambda12_result <- optimize_lambda(lambda12, fixed_lambda11 = optimal.lambda11, fixed_lambda12 = NULL,
                                     fixed_lambda2 = 1e-3)
  optimal.lambda12 <- lambda12_result$optimal_lambda
  lambda2_result <- optimize_lambda(lambda2, fixed_lambda11 = optimal.lambda11, 
                                    fixed_lambda12 = optimal.lambda12, fixed_lambda2 = NULL)
  optimal.lambda2 <- lambda2_result$optimal_lambda
  
  return(list(lambda11 = optimal.lambda11, lambda12 = optimal.lambda12, lambda2 = optimal.lambda2))
}

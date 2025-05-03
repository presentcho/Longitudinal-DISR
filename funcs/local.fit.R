# local.fit: fits a local estimators for each subregions

local.fit <- function(iter, Ver0, Tr0, TV0, n.layer, X, Y, V, M, tij, d, r, L, rho, time.knots, time.bound,
                      lambda11 = exp(seq(log(0.001), log(1000), length.out = 5)), 
                      lambda12 = exp(seq(log(0.001), log(1000), length.out = 5)), 
                      lambda2 = exp(seq(log(0.001), log(1000), length.out = 5))) {
  # Generate basis function for subregions
  nbasis.tri <- (d+1)*(d+2)/2
  result <- ring.dc(iter, Ver0, Tr0, TV0, n.layer)
  Basis.local <- basis.tensor.local(ss = V, tt = tij, V = result$V1, Tri = result$Tr1, 
                                    d = d, r = r, time.knots = time.knots, time.bound = time.bound, rho = rho)
  U0.s <- Basis.local$U0
  B0.s <- Basis.local$B0
  Q2.s <- Basis.local$Q2
  K.s <- Basis.local$K
  P11.s <- Basis.local$D1
  P12.s <- Basis.local$D2
  P2.s <- kronecker(diag(ncol(X)), as.matrix(crossprod(Q2.s, K.s) %*% Q2.s))
  ind.s <- Basis.local$ind.inside
  X <- X
  Y <- Y[,ind.s]
  V <- V[ind.s,]
  M <- M
  tij <- tij
  # Finds tuning parameters in each subregion
  gcv.lambda <- gc.fit.gcv(X=X, Y=Y, V=V, M=M, U0=U0.s, B0=B0.s, Q2=Q2.s, K=K.s, P11=P11.s, P12=P12.s, P2=P2.s, 
                           lambda11 = exp(seq(log(0.001), log(1000), length.out = 5)), 
                           lambda12 = exp(seq(log(0.001), log(1000), length.out = 5)), 
                           lambda2 = exp(seq(log(0.001), log(1000), length.out = 5)))
  lambda11.optimal <- gcv.lambda$lambda11
  lambda12.optimal <- gcv.lambda$lambda12
  lambda2.optimal <- gcv.lambda$lambda2
  # Fits the estimators 
  mfit0 <- gc.fit(X, Y, V, M, U0.s, B0.s, Q2.s, K.s, P11.s, P12.s, P2.s, L, rho,
                   lambda11=lambda11.optimal, lambda12=lambda12.optimal, lambda2=lambda2.optimal)
  idx.tr <- prodlim::row.match(data.frame(result$V.tr), data.frame(Tr0))
  idx.psi <- rep((idx.tr - 1) * nbasis.tri, each = nbasis.tri) + rep(1:nbasis.tri, times = length(idx.tr))  
  lambdac <- c(lambda11.optimal, lambda12.optimal, lambda2.optimal)
  list(idx.tr = idx.tr, idx.psi = idx.psi, gamma.hat = mfit0$gamma.hat, theta.hat = mfit0$theta.hat, lambdac = lambdac)
}

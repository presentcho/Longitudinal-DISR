local.fit.k <- function(iter, Ver0, Tr0, TV0, n.layer, X, Y, V, M, tij, 
                        d, r, L, rho, time.knots, time.bound, lambda11, lambda12, lambda2
) {
  # generate basis function
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
  mfit0 <- gc.fit(X, Y, V, M, U0.s, B0.s, Q2.s, K.s, P11.s, P12.s, P2.s, L, rho,
                   lambda11=lambda11, lambda12=lambda12, lambda2=lambda2)
  idx.tr <- prodlim::row.match(data.frame(result$V.tr), data.frame(Tr0))
  idx.psi <- rep((idx.tr - 1) * nbasis.tri, each = nbasis.tri) + rep(1:nbasis.tri, times = length(idx.tr))  
  lambdac <- c(lambda11, lambda12, lambda2)
  list(idx.tr = idx.tr, idx.psi = idx.psi, gamma.hat = mfit0$gamma.hat, theta.hat = mfit0$theta.hat, lambdac = lambdac)
}

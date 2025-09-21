local.fit.bivar <- function(iter, V0, Tr0, TV0, n.layer = 2, Y.all, X.all, Z.all, 
                      d = 5, r = 1, lambda = lambda) {
  nbasis.tri <- (d + 1) * (d + 2) / 2
  result <- ring.dc(iter, V0, Tr0, TV0, n.layer = n.layer)
  
  B.local <- basis(result$V1,result$Tr1, d, r, Z.all)
  K.local <- B.local$K
  Q2.local <- B.local$Q2
  B.mtr.local <- B.local$B
  ind.local <- B.local$Ind.inside
  Y.sub <- Y.all[,ind.local]
  mfit0 <- bivar.gc.fit(B.mtr.local, Q2.local, K.local, lambda,
                        X.all, Y.sub)
  idx.tr <- prodlim::row.match(data.frame(result$V.tr), data.frame(Tr0))
  idx.gamma <- rep((idx.tr - 1) * nbasis.tri, each = nbasis.tri) + 
    rep(1:nbasis.tri, times = length(idx.tr))
  
  list(idx.tr = idx.tr, idx.gamma = idx.gamma, gamma.local = mfit0$gamma)
}
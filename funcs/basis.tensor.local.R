#basis.tensor.local: generates bivariate and trivariate spline basis matrix in each subregions

basis.tensor.local <- function (ss, tt, V, Tri, d = 2, r = 1, time.knots, time.bound, 
                                rho = 3){ 
  ss <- as.matrix(ss)
  ss.uni <- unique(ss)
  index.ss <- match(data.frame(t(ss)), data.frame(t(ss.uni)))
  B.uni <- BPST::basis(V, Tri, d, r, ss.uni)
  B0 <- B.uni$B
  Q2 <- B.uni$Q2
  K <- B.uni$K
  BQ2.uni <- B0 %*% Q2
  BQ2 <- BQ2.uni
  tt.uni <- unique(tt)
  index.tt <- match(tt, tt.uni)
  U0.uni <- bSpline(tt.uni, knots = time.knots, intercept = FALSE, 
                    degree = rho, Boundary.knots = time.bound)
  U0 <- U0.uni[index.tt, ]
  U0 <- Matrix(U0, sparse = TRUE)
  B0 <- Matrix(B0, sparse = TRUE)
  dimU <- ncol(U0)
  Q2.all <- kronecker(diag(dimU), Q2)
  H <- B.uni$H
  H.all <- kronecker(diag(dimU), H)
  Energ <- energy.tensor(V, Tri, d, time.knots, degr = rho, time.bound)
  P1 <- Energ$Eng1
  P2 <- Energ$Eng2
  D1 <- crossprod(Q2.all, as.matrix(P1)) %*% Q2.all
  D2 <- crossprod(Q2.all, as.matrix(P2)) %*% Q2.all
  list(Q2 = Q2, H = H, K = K, Q2.all = Q2.all, U0 = U0, B0 = B0, 
       H.all = H.all, dimB = ncol(B0), dimU = dimU, P1 = P1, 
       P2 = P2, D1 = D1, D2 = D2, ind.inside = B.uni$Ind.inside)
}

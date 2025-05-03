# energy.tensor: generates bivariate and trivariate energy functions

energy.tensor <- function (V, Tr, d, time.knots, degr, tt) 
{
  K <- TPST:::energy(V, Tr, d)
  UU <- bSpline(seq(range(tt)[1], range(tt)[2], 0.01), knots = time.knots, 
                intercept = FALSE, degree = degr)
  energy.U <- 0.01 * t(UU) %*% UU
  UU2 <- dbs(seq(range(tt)[1], range(tt)[2], 0.01), derivs = 2, 
             knots = time.knots, degree = degr, intercept = FALSE)
  energy.U2 <- 0.01 * t(UU2) %*% UU2
  list(Eng1 = kronecker(energy.U, K$K2deri), Eng2 = kronecker(energy.U2, 
                                                              K$K2))
}

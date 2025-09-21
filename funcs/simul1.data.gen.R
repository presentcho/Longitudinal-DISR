simul.data.gen <- function(n, Tr, Ver, V) {
  ind <- inVT(Ver, Tr, V[, 1], V[, 2])
  ind.inside <- ind$ind.inside
  v1 <- V[, 1]
  v2 <- V[, 2]
  n1 <- length(unique(v1))
  n2 <- length(unique(v2))
  
  # Randomly generate time points 
  mi <- sample(2:7, n, replace = TRUE)
  tij <- list()
  for(i in 1:n){
    ti1 <- runif(1, min = 0.5, max = 1)
    deltaij <- sapply(mi[i] - 1, function(m) runif(m, min = 0.5, max = 1))
    tij[[i]] <- cumsum(c(ti1, deltaij))
  }
  tij <- unlist(tij)
  
  # Generate covariate matrix X 
  mu <- rep(0, 2)
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  X <- scale(mvrnorm(sum(mi), mu, Sigma))
  X.base <- cbind(1, X)
  X.slope <- tij*X.base
  # Define alpha coefficient functions 
  alpha1 <- 3 * cos(0.06 * pi * ((-1.5 * v1^2) + (-1.5 * v2^2)))
  alpha2 <- 2 * sin(pi * v2) - 1
  alpha3 <- 4 * cos(pi * v2)
  bivar.alpha <- cbind(alpha1, alpha2, alpha3)
  ind.inside.alpha <- which(!is.na(alpha1))
  # Define beta coefficient functions 
  beta1 <- mgcv::fs.test(v1,v2)
  beta2 <- 2 * cos(pi * v2) - 1
  beta3 <- (0.5 * (v1 - 0.5)^2 + (v2 - 0.5)^2) + 2 * cos(pi* v2) - 1 
  bivar.beta <- cbind(beta1, beta2, beta3)
  ind.inside.beta <- which(!is.na(beta1))
  
  idx <- intersect(ind.inside, ind.inside.alpha)
  idx <- intersect(idx, ind.inside.beta)
  ind.outside <- setdiff(1:(n1 * n2), idx)
  bivar.alpha[ind.outside, ] <- NA
  bivar.beta[ind.outside, ] <- NA
  
  # Set alpha coefficient functions error terms
  lambda1 <- 0.5413 * sin(pi * v1)
  lambda2 <- 0.7805 * sin(pi * v2)
  vartheta1 <- rnorm(n, mean = 0, sd = 1)
  vartheta2 <- rnorm(n, mean = 0, sd = 1)
  Xi <- matrix(0, nrow = n, ncol = n1*n2)
  for(i in 1:n){
    Xi[i,] <- vartheta1[i]*lambda1 + vartheta2[i]*lambda2 
  }
  Xij <- do.call(rbind, lapply(1:n, function(i) Xi[i, , drop = FALSE][rep(1, mi[i]), ]))
  Xij[,ind.outside] <- NA
  
  # Set beta coefficient functions error terms
  omega1 <- 0.6472 * cos(pi * v1)
  omega2 <- 0.7237 * cos(pi * v2)
  rho1 <- rnorm(n, mean = 0, sd = 0.5)
  rho2 <- rnorm(n, mean = 0, sd = 0.5)
  Zetai <- matrix(0, nrow = n, ncol = n1*n2)
  for(i in 1:n){
    Zetai[i,] <- rho1[i]*omega1 + rho2[i]*omega2 
  }
  Zetaij <- do.call(rbind, lapply(1:n, function(i) Zetai[i, , drop = FALSE][rep(1, mi[i]), ]))
  Zetaij[,ind.outside] <- NA
  
  # Generate baseline function $\mu_i(v)$ and slope function $\eta_i(v)$
  Muij <- X.base %*% t(bivar.alpha) + Xij
  Etaij <- X.slope %*% t(bivar.beta) + tij * Zetaij 
  error.term <- matrix(0, nrow = sum(mi), ncol = n1*n2)
  for(ij in 1:sum(mi)){
    error.term[ij,] <- rnorm(n1*n2, mean = 0, sd = 0.5)
  }
  error.term[,ind.outside] <- NA
  
  Yij <- Muij + Etaij + error.term
  X.mat <- cbind(X.base, X.slope)
  M.vec <- sqrt(1/rep(mi, times = mi))
  Y.mat <- Yij
  X.mat <- cbind(X.base, X.slope)
  V.mat <- V
  list(Y = Y.mat, X = X.mat, M.vec = M.vec, tij = tij,
       V = V.mat, bivar.alpha = bivar.alpha,
       bivar.beta = bivar.beta, ind.inside = idx)
}

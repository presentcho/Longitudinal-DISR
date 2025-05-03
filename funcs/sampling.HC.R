# sampling.HC: identifies the index of a triangle within a triangulation using Hilbert space filling-curve.

sampling.HC <- function(n, V, Tr, n.layer = 2){
  points.tri <- c()
  for(iter in 1:nrow(Tr)){
    points.tri <- rbind(points.tri, colMeans(V[Tr[iter,], ]))
  }
  covered.idx = FALSE
  covered.idx.all = c()
  time.sample = 0
  idx.all <- c()
  mean.cover.all <- c()
  sd.cover.all <- c()
  
  while(time.sample < 5) {
    time.sample = time.sample + 1
    i <- 4
    df <- tibble(val = 1:4^i,
                 size = runif(4^i, 1, 5),
                 color = rep(c(1, 2, 3, 4), 4^(i - 1)))
    values <- df$val
    
    ggh <- as.matrix(hilbertd2xy(n = 4^i, values)) + 1
    grid.x <- seq(min(points.tri[,1]), max(points.tri[,1]), length.out = 2^i)
    grid.y <- seq(min(points.tri[,2]), max(points.tri[,2]), length.out = 2^i)
    grid.all <- cbind(grid.x[ggh[, 2]], grid.y[ggh[, 1]])
    
    idx <- apply(points.tri, 1, function(x) which.min((x[1] - grid.all[ , 1])^2 + 
                                                        (x[2] - grid.all[ , 2])^2))
    cut.region <- cut(idx/4^i, seq(0, 1, 1/n))
    s.all <- c()
    for(iter in levels(cut.region)){
      s1 <- which(cut.region == iter)
      if (length(s1) >= 1){
        s1 <- sample(which(cut.region == iter), 1)
        s.all <- c(s.all, s1)
      }
    }
    idx <- s.all 
    idx
    
    TV <- tdata(V, Tr)$TV
    covered <- c()
    TriPlot(V, Tr)
    points(points.tri[idx, ])
    for(iter in idx) {
      result <- ring.dc(iter, V, Tr, TV, n.layer = n.layer)
      covered <- c(covered, prodlim::row.match(data.frame(result$V.tr), data.frame(Tr)))
    }
    covered.idx = sum(1:nrow(Tr) %in% unique(covered)) == nrow(Tr)
    idx.all = cbind(idx.all, idx)
    covered.idx.all = c(covered.idx.all, covered.idx)
    mean.cover.all = c(mean.cover.all, mean(table(covered)))
    sd.cover.all = c(sd.cover.all, sd(table(covered)))
  }
  
  if(!any(covered.idx.all == TRUE)) {
    cat("Triangles are not fully covered; Try to increase sample size or n.layer.")
  } else {
    best.sample = which.min(sd.cover.all[covered.idx.all])
    best.sample = which(covered.idx.all)[best.sample]
    cat("selected set:", best.sample, "\n")
    return(list(sample.tri = idx.all[, best.sample ], 
                mean.cover = mean.cover.all[best.sample],
                sd.cover = sd.cover.all[best.sample]))
  }
}

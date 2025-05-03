# ring.dc: identifies the neighborhood of a triangle within a triangulation.

ring.dc <- function(ind.T, V, Tr, TV, n.layer){
  Tr0 = Tr[ind.T, ]
  I.tr = c(Tr0)
  V1 = c(); Tr1 = c();
  for(i in 1:n.layer){
    k = length(I.tr)
    J = c()
    for(j in 1:k){
      J1 = which(TV[, I.tr[j]] != 0)
      J = c(J, J1)
    }
    V2 = Tr[J, ]
    I.tr = c(V2)
    V1 = c(V1, I.tr)
    Tr1 = c(Tr1, J)
  }
  V1 = sort(unique(V1))
  I.tr = V1
  V1 = V[I.tr, ]
  V.tr = Tr[sort(unique(Tr1)), ]
  Tr1 = Tr[sort(unique(Tr1)), ]
  nt = nrow(Tr1)
  for(i in 1:nt){
    j = Tr1[i, 1]; nj = which(I.tr == j); Tr1[i, 1] = nj
    j = Tr1[i, 2]; nj = which(I.tr == j); Tr1[i, 2] = nj
    j = Tr1[i, 3]; nj = which(I.tr == j); Tr1[i, 3] = nj
  }
  j = Tr0[1]; nj = which(I.tr == j); Tr0[1] = nj
  j = Tr0[2]; nj = which(I.tr == j); Tr0[2] = nj
  j = Tr0[3]; nj = which(I.tr == j); Tr0[3] = nj
  
  ring.tri = list(V1 = V1, Tr1 = Tr1, Tr0 = Tr0, I.tr = I.tr, V.tr = V.tr)
}

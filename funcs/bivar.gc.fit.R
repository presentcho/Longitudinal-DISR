bivar.gc.fit =function(B,Q2,K,lambda,X,Y){
  
  n=nrow(Y)
  nx=ncol(X)
  J=ncol(Q2)
  npix=nrow(B)
  BQ2=B%*%Q2
  
  WW=kronecker(crossprod(X),crossprod(BQ2))
  
  rhs=rowSums(sapply(1:n, function(iter) as.matrix(kronecker(X[iter,],
                                                             crossprod(BQ2,Y[iter,])))))
  
  P=as.matrix(crossprod(Q2,K)%*%Q2)
  
  flag=(rankMatrix(WW)<(J*nx))
  if(!flag){
    Ainv=chol(WW,pivot=TRUE)
  }
  if(flag){
    warning("The solution may not be correct.")
    VV=WW+diag(rep(1e-6),nx*J)
    Ainv=chol(VV,pivot=TRUE)
  }
  A=solve(t(Ainv))
  D=as.matrix(kronecker(diag(rep(1,nx)),P))
  ADA=A%*%D%*%t(A)
  eigs=eigen(ADA)
  C=eigs$values
  
  lambda=as.matrix(lambda)
  nlam=nrow(lambda)
  if(ncol(lambda)==1){
    lambda=matrix(rep(lambda,times=nx),nlam,nx)
  }
  
  sse_all=c(); df_all=c(); gcv_all=c(); bic_all=c()
  for(il in 1:nlam){
    if(nx==1){
      Lam=as.matrix(lambda[il])
    }
    if(nx>1){
      Lam=diag(lambda[il,])
    }
    Dlam=as.matrix(kronecker(Lam,P))
    lhs=WW+Dlam
    theta=solve(lhs,rhs)
    theta.mtx=matrix(theta,J,nx)
    beta=BQ2%*%theta.mtx
    Yhat=tcrossprod(X,beta) 
    
    ssei=apply((Y-Yhat)^2,1,sum)
    sse=sum(ssei)
    sse_all=c(sse_all,sse)
    
    df=sum(1/(1+C*mean(lambda[il,])))
    df_all=c(df_all,df)
    gcv=n*sse/(n-df)^2
    gcv_all=c(gcv_all,gcv)
    bic=log(sse/n)+df*log(n)/n
    bic_all=c(bic_all,bic)
    
  }
  j=which.min(gcv_all)
  lambdac=lambda[j,]
  df=df_all[j]
  sse=sse_all[j]
  gcv=gcv_all[j]
  bic=bic_all[j]
  
  if(nx==1){
    Lam=as.matrix(lambdac)
  }
  if(nx>1){
    Lam=diag(lambdac)
  }
  Dlam=as.matrix(kronecker(Lam,P))
  lhs=WW+Dlam
  theta=solve(lhs,rhs)
  theta.mtx=matrix(theta,J,nx)
  gamma=as.matrix(Q2%*%theta.mtx)
  beta=as.matrix(BQ2%*%theta.mtx)
  Yhat=as.matrix(X%*%t(beta))
  
  list(beta=beta,gamma=gamma,Yhat=Yhat,sse=sse,df=df,gcv=gcv,
       bic=bic,lambdac=lambdac)
}